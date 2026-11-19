"""Fit Treasury curves, solve the continuous-time HJB, and replay ETF closes."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import make_smoothing_spline

from download import DATA, SOURCES

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
TENORS = dict(zip(
    ["1 Mo", "2 Mo", "3 Mo", "6 Mo", "1 Yr", "2 Yr", "3 Yr", "5 Yr", "7 Yr",
     "10 Yr", "20 Yr", "30 Yr"],
    [1/12, 2/12, 3/12, 0.5, 1, 2, 3, 5, 7, 10, 20, 30],
))
DT = 1 / 252
A_RATE = 0.05
PROXY_MATURITIES = (8.5, 25)
HORIZON_DAYS = 63
FEE = 0.0002
WARMUP = 61
BASE_DAYS = 5
PENALTY_DAYS = [2, 5, 10, 20]


def read_treasury():
    frames = [pd.read_csv(DATA / f"treasury_{year}.csv", parse_dates=["Date"])
              for year in range(2022, 2026)]
    rates = pd.concat(frames).set_index("Date").sort_index()[list(TENORS)] / 100
    if rates.index.has_duplicates or rates.isna().any().any():
        raise ValueError("Expected unique, complete Treasury observations at selected tenors")
    return rates


def read_etf(ticker):
    chart = json.loads((DATA / f"{ticker}.json").read_text(encoding="utf-8"))["chart"]
    raw, = chart["result"]
    if raw["meta"]["symbol"] != ticker or raw.get("events", {}).get("splits"):
        raise ValueError(f"{ticker}: unexpected symbol or split requiring share adjustment")

    def dates(stamps):
        return pd.to_datetime(stamps, unit="s", utc=True).tz_convert(
            "America/New_York"
        ).tz_localize(None).normalize()

    frame = pd.DataFrame({
        "yahoo_close": raw["indicators"]["quote"][0]["close"],
        "adjusted": raw["indicators"]["adjclose"][0]["adjclose"],
    }, index=dates(raw["timestamp"])).sort_index()
    exchange = json.loads((DATA / f"{ticker}-nasdaq.json").read_text(encoding="utf-8"))["data"]
    closes = pd.DataFrame(exchange["tradesTable"]["rows"])
    closes.index = pd.to_datetime(closes["date"], format="%m/%d/%Y")
    if exchange["symbol"] != ticker or not closes.index.sort_values().equals(frame.index):
        raise ValueError(f"{ticker}: inconsistent close histories")
    frame["close"] = closes["close"].str.replace(",", "").astype(float)
    events = json.loads((DATA / f"{ticker}-dividends.json").read_text(
        encoding="utf-8"
    ))["data"]["dividends"]["rows"]
    cash = pd.DataFrame(events)
    cash.index = pd.to_datetime(cash["exOrEffDate"], format="%m/%d/%Y")
    cash = cash.loc[(cash.index >= frame.index.min()) & (cash.index <= frame.index.max())]
    dividends = pd.Series(
        cash["amount"].str.replace("$", "", regex=False).astype(float),
    ).groupby(level=0).sum()
    yahoo_events = raw.get("events", {}).get("dividends", {}).values()
    yahoo_dates = dates([event["date"] for event in yahoo_events])
    if (frame.index.has_duplicates or frame.isna().any().any()
            or (frame <= 0).any().any() or dividends.empty
            or not dividends.index.isin(frame.index).all()
            or not dividends.index.equals(yahoo_dates.sort_values().unique())):
        raise ValueError(f"{ticker}: invalid prices or distribution date")
    frame["cash"] = dividends.reindex(frame.index, fill_value=0)
    frame["pnl"] = frame["close"].diff() + frame["cash"]
    frame["adjustment_gap"] = (
        frame["pnl"] / frame["close"].shift() - frame["adjusted"].pct_change()
    )
    return frame


def curve(tenors, par, multiplier=None):
    # The short end uses an annualized simple-interest approximation.
    # Semiannual coupon cash flows fill the discount grid from the smoothed par rates.
    half_years = np.arange(0.5, 30.5, 0.5)
    par_spline = make_smoothing_spline(tenors, par)
    discounts = []
    for coupon in par_spline(half_years) / 2:
        discounts.append((1 - coupon * sum(discounts)) / (1 + coupon))
    zero_tenors = np.r_[tenors[tenors < 0.5], half_years]
    zero_discounts = np.r_[1 / (1 + tenors[tenors < 0.5] * par[tenors < 0.5]),
                           discounts]
    if min(zero_discounts) <= 0:
        raise ValueError("Nonpositive bootstrapped discount factor")
    zero_nodes = -np.log(zero_discounts) / zero_tenors
    zero_spline = make_smoothing_spline(
        zero_tenors, zero_nodes,
        lam=None if multiplier is None else multiplier * gcv_lambda(zero_tenors, zero_nodes),
    )
    return par_spline, zero_spline, zero_tenors, zero_nodes


def gcv_lambda(x, y):
    # Recover the GCV-selected lambda from the fitted spline's normal equations.
    spline = make_smoothing_spline(x, y)
    second = spline(x, 2)
    # Test the normal equations against g minus a linear function for conditioning.
    fitted = spline(x)
    departure = fitted - fitted[0]
    departure -= (fitted[-1] - fitted[0]) * (x - x[0]) / (x[-1] - x[0])
    roughness = np.sum(np.diff(x) * (second[:-1]**2 + second[:-1]*second[1:]
                                    + second[1:]**2) / 3)
    return max(0.0, np.dot(y - fitted, departure) / roughness) if roughness else 0.0


def zero_value(spline, t, order=0):
    first = spline.t[0]
    t = np.asarray(t)
    value = spline(np.maximum(t, first), order)
    if order == 0:
        return value + np.minimum(t - first, 0) * spline(first, 1)
    if order == 2:
        return np.where(t < first, 0, value)
    return value


def forward(zero, t):
    return zero_value(zero, t) + t * zero_value(zero, t, 1)


def loading(t):
    return -np.expm1(-A_RATE * t) / A_RATE


def theta(zero, t, sigma):
    return (2 * zero_value(zero, t, 1) + t * zero_value(zero, t, 2) + A_RATE * forward(zero, t)
            + sigma**2 * loading(2*t) / 2)


def reconstruct_discount(zero, t, sigma):
    # Conditional mean and variance of the integrated Gaussian short rate.
    knots = np.unique(zero.t)
    interior = knots[(knots > 0) & (knots < t)]
    mean = loading(t) * forward(zero, 0) + quad(
        lambda s: loading(t - s) * theta(zero, s, sigma),
        0, t, points=interior, limit=150,
    )[0]
    variance = sigma**2 * quad(lambda s: loading(s)**2, 0, t)[0]
    return np.exp(-mean + variance / 2)


def riccati(beta, sigma, penalty, horizon, terminal):
    left, hedge = beta
    k = sigma**2

    def rhs(_, coefficients):
        a, b, _ = coefficients
        return [k * hedge**2 - a*a/penalty,
                2*k * left*hedge - a*b/penalty,
                k * left**2 - b*b/(4*penalty)]

    solution = solve_ivp(rhs, (0, horizon), [terminal, 0, 0],
                         rtol=1e-8, atol=1e-10)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[:, -1]


def replay(market, rates, diagnostics, *, penalty_reference=None):
    dates = market.index
    signal_dates = rates.index.to_series().reindex(dates, method="ffill")
    signal_rates = rates.reindex(dates, method="ffill")
    sigma = signal_rates["3 Mo"].diff().rolling(60, min_periods=60).std() / np.sqrt(DT)
    prices = market.xs("close", axis=1, level=1).to_numpy()
    asset_pnl = market.xs("pnl", axis=1, level=1).to_numpy()
    scale = prices[WARMUP - 1, 0]
    fitted = diagnostics.reindex(dates, method="ffill")
    hw_beta = -prices * fitted[["IEF_hw", "TLT_hw"]].to_numpy() / scale
    duration_beta = -prices * fitted[["IEF_duration", "TLT_duration"]].to_numpy() / scale
    # A comparison can hold both penalty schedules at their baseline values.
    penalty_beta_h = hw_beta[:, 1]
    if penalty_reference is not None:
        reference = penalty_reference["TLT_hw"].reindex(dates, method="ffill")
        penalty_beta_h = -prices[:, 1] * reference.to_numpy() / scale
    sigma0 = sigma.iloc[WARMUP - 1]
    penalties = {days: (sigma0 * penalty_beta_h[WARMUP - 1] * days * DT)**2
                 for days in PENALTY_DAYS}
    positions = {"No hedge": 0., "Static duration": 0., "Dynamic duration": 0.,
                 **{f"HJB {days}d": 0. for days in PENALTY_DAYS},
                 "HJB matched duration": 0.}
    rows = []
    for i in range(WARMUP, len(dates) - 1):
        j = i - 1
        duration_ratio = -duration_beta[j, 0] / duration_beta[j, 1]
        static = duration_ratio if i == WARMUP else positions["Static duration"]
        targets = {"No hedge": 0., "Static duration": static, "Dynamic duration": duration_ratio}
        remaining_days = len(dates) - i - 1
        horizon_days = min(HORIZON_DAYS, remaining_days)
        energy = {}
        for name in positions:
            if not name.startswith("HJB"):
                continue
            days = BASE_DAYS if name == "HJB matched duration" else int(name.split()[1][:-1])
            penalty = penalties[days]
            beta = hw_beta[j].copy()
            if name == "HJB matched duration":
                beta[0] = -duration_ratio * beta[1]
            # A rolling horizon freezes the latest loadings and volatility.
            # Penalize terminal inventory only when the actual end is in view.
            terminal = (np.sqrt(penalty) * sigma.iloc[j] * abs(penalty_beta_h[j])
                        if horizon_days == remaining_days else 0.)
            a, b, _ = riccati(beta, sigma.iloc[j], penalty, horizon_days * DT, terminal)
            rate = -(2*a * positions[name] + b) / (2*penalty)
            gain = a / penalty
            # Hold the beginning-of-step feedback coefficients over the session.
            if gain:
                targets[name] = positions[name] - rate * np.expm1(-gain * DT) / gain
                energy[name] = -rate**2 * np.expm1(-2 * gain * DT) / (2 * gain)
            else:
                targets[name] = positions[name] + rate * DT
                energy[name] = rate**2 * DT
        for name, q in targets.items():
            dq = q - positions[name]
            closing = i == len(dates) - 2
            turnover = abs(dq) * prices[i, 1] / scale
            if closing:
                turnover += abs(q) * prices[i + 1, 1] / scale
            gross = (asset_pnl[i + 1, 0] + q * asset_pnl[i + 1, 1]) / scale
            rows.append({
                "date": dates[i+1], "execution_date": dates[i], "signal_date": dates[j],
                "treasury_date": signal_dates.iloc[j], "strategy": name, "position": q,
                "gross": gross, "net": gross - FEE * turnover,
                "turnover": turnover, "shares_traded": abs(dq) + closing * abs(q),
                "rate_energy": energy.get(name, dq*dq / DT), "sigma": sigma.iloc[j],
                "hw_variance": (sigma.iloc[j] * (hw_beta[j, 0] + q * hw_beta[j, 1]))**2 * DT,
            })
        positions.update(targets)
    return pd.DataFrame(rows), penalties


def summarize(daily, penalties):
    groups = daily.groupby("strategy", sort=False)
    summary = groups.agg(
        gross_sd=("gross", "std"), net_sd=("net", "std"), gross_total=("gross", "sum"),
        net_total=("net", "sum"), turnover=("turnover", "sum"),
        shares_traded=("shares_traded", "sum"), rate_energy=("rate_energy", "sum"),
        hw_variance=("hw_variance", "sum"),
    )
    summary["gross_rmse"] = np.sqrt(groups["gross"].apply(lambda x: np.mean(x*x)))
    summary["net_rmse"] = np.sqrt(groups["net"].apply(lambda x: np.mean(x*x)))
    summary["variance_reduction"] = 1 - (summary["gross_sd"] / summary.loc["No hedge", "gross_sd"])**2
    summary["cost"] = FEE * summary["turnover"]
    summary["rate_penalty_common"] = penalties[BASE_DAYS] * summary["rate_energy"]
    summary["lambda"] = [
        penalties[BASE_DAYS] if name == "HJB matched duration" else
        penalties[int(name.split()[1][:-1])] if name.startswith("HJB") else np.nan
        for name in summary.index
    ]
    summary["rate_penalty_own"] = summary["lambda"] * summary["rate_energy"]
    return summary


def paired_risk_intervals(daily):
    """Resample aligned strategy P&Ls in circular, consecutive session blocks."""
    panel = daily.pivot(index="date", columns="strategy", values="gross")
    names = ["Dynamic duration", "HJB 5d", "HJB matched duration"]
    values = panel[names].to_numpy()
    rng = np.random.default_rng(2025)
    rows = []
    for length in (10, 20, 40):
        n = len(panel)
        blocks = (n + length - 1) // length
        starts = rng.integers(n, size=(4000, blocks))
        samples = ((starts[..., None] + np.arange(length)) % n).reshape(4000, -1)[:, :n]
        vols = 1e4 * values[samples].std(axis=1, ddof=1)
        observed = 1e4 * values.std(axis=0, ddof=1)
        for j, name in enumerate(names[1:], 1):
            low, high = np.quantile(vols[:, j] - vols[:, 0], (0.025, 0.975))
            rows.append({"strategy": name, "block_days": length,
                         "point_bp": observed[j] - observed[0], "low_bp": low, "high_bp": high})
    return pd.DataFrame(rows)


def fit_curves(rates):
    grid = np.linspace(min(TENORS.values()), 30, 601)
    tenors = np.array(list(TENORS.values()))
    example = rates.index[-1]
    rows = []
    for date, row in rates.iterrows():
        par, zero, nodes, zero_nodes = curve(tenors, row.to_numpy())
        fwd = forward(zero, grid)
        fitted = {
            "date": date,
            "par_rmse_bp": 1e4 * np.sqrt(np.mean((par(tenors) - row)**2)),
            "zero_rmse_bp": 1e4 * np.sqrt(np.mean((zero(nodes) - zero_nodes)**2)),
            "min_forward": min(fwd),
            "max_forward": max(fwd),
            "forward_turns": np.count_nonzero(np.diff(np.sign(np.diff(fwd)))),
            "min_discount": min(np.exp(-grid * zero(grid))),
            "theta_0_without_volatility": theta(zero, 0., 0.),
        }
        for ticker, maturity in zip(("IEF", "TLT"), PROXY_MATURITIES):
            payments = np.arange(0.5, maturity + 0.5, 0.5)
            cashflows = np.full_like(payments, par(maturity) / 2)
            cashflows[-1] += 1
            pv = cashflows * np.exp(-payments * zero(payments))
            fitted[f"{ticker}_duration"] = np.average(payments, weights=pv)
            fitted[f"{ticker}_hw"] = np.average(loading(payments), weights=pv)
        rows.append(fitted)
        if date == example:
            example_curves = par, zero
    return pd.DataFrame(rows).set_index("date"), example_curves


def save_plot(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_curves(rates, example_curves):
    tenors = np.array(list(TENORS.values()))
    maturity = np.linspace(tenors[0], tenors[-1], 601)
    observed = rates.iloc[-1].to_numpy()
    par, zero = example_curves
    fig, ax = plt.subplots()
    ax.plot(maturity, 100 * par(maturity), label="Smoothed par yield")
    ax.scatter(tenors, 100 * observed, label="Observed par yield", s=17)
    ax.plot(maturity, 100 * zero(maturity), "--", label="Smoothed zero yield")
    ax.set(xlabel="Maturity (years)", ylabel="Yield (%)")
    ax.legend()
    save_plot(fig, "yield")

    fig, ax = plt.subplots()
    sensitivity = {}
    for factor in (0.25, 1, 4):
        if factor != 1:
            _, spline, _, _ = curve(tenors, observed, factor)
        else:
            spline = zero
        fitted = forward(spline, maturity)
        sensitivity[str(factor)] = float(np.max(np.abs(fitted - forward(zero, maturity))) * 1e4)
        ax.plot(maturity, 100 * fitted, label=f"{factor:g} x GCV penalty")
    ax.set(xlabel="Maturity (years)", ylabel="Instantaneous forward rate (%)")
    ax.legend()
    save_plot(fig, "forward")
    return sensitivity


def plot_strategies(daily, summary):
    main = ["No hedge", "Static duration", "Dynamic duration", "HJB 5d"]
    for field, name, ylabel in (
        ("position", "position", "TLT shares per IEF share"),
        ("net", "pnl", "Cumulative net P&L (% of initial IEF value)"),
    ):
        fig, ax = plt.subplots()
        for strategy in main:
            frame = daily[daily["strategy"] == strategy]
            if field == "position" and strategy == "No hedge":
                continue
            values = frame[field] if field == "position" else 100 * frame[field].cumsum()
            dates = frame["execution_date"] if field == "position" else frame["date"]
            ax.plot(dates, values, label=strategy)
        ax.set(ylabel=ylabel)
        ax.legend()
        save_plot(fig, name)

    fig, (risk_ax, energy_ax) = plt.subplots(2, 1, sharex=True)
    hjb = summary.loc[[f"HJB {days}d" for days in PENALTY_DAYS]]
    x = hjb["lambda"]
    risk_ax.semilogx(x, 1e4 * hjb["gross_sd"], "o-", label="Hull–White HJB")
    risk_ax.axhline(1e4 * summary.loc["Dynamic duration", "gross_sd"],
                    linestyle="--", color="gray", label="Dynamic duration")
    risk_ax.set(ylabel="Daily P&L SD (bp)")
    risk_ax.legend()
    energy_ax.semilogx(x, hjb["rate_penalty_own"], "o-", label=r"Own $\lambda$")
    energy_ax.semilogx(x, hjb["rate_penalty_common"], "s--", label=r"Common $\lambda_{5d}$")
    energy_ax.set(xlabel=r"Trading-rate penalty $\lambda$", ylabel="Realized penalty")
    energy_ax.legend()
    save_plot(fig, "sensitivity")


def write_report_numbers(summary, daily, audit, intervals):
    names = {"No hedge": "Unhedged", "Static duration": "Static duration",
             "Dynamic duration": "Dynamic duration", "HJB 5d": "Hull--White HJB (5d)",
             "HJB matched duration": "Duration-target HJB (5d)"}
    br = r" \\"
    rows = []
    costs = []
    for name, label in names.items():
        row = summary.loc[name]
        rows.append(f"{label} & {1e4*row.gross_sd:.2f} & {1e4*row.net_sd:.2f} & "
                    f"{1e4*row.net_rmse:.2f} & {100*row.variance_reduction:.2f} & "
                    f"{100*row.gross_total:.2f} & {100*row.net_total:.2f}" + br)
        costs.append(f"{label} & {100*row.turnover:.2f} & {1e4*row.cost:.2f} & "
                     f"{row.rate_energy:.3f} & {row.rate_penalty_common:.6f}" + br)
    sensitivity = []
    for days in PENALTY_DAYS:
        row = summary.loc[f"HJB {days}d"]
        sensitivity.append(f"{days} & {row['lambda']:.7f} & {1e4*row.gross_sd:.2f} & "
                           f"{100*row.turnover:.2f} & {row.rate_energy:.3f} & "
                           f"{row.rate_penalty_own:.6f} & {100*row.net_total:.2f}" + br)
    for name, content in (("risk", rows), ("costs", costs), ("penalties", sensitivity)):
        (OUT / f"{name}.tex").write_text("\n".join(content), encoding="utf-8")
    macros = {
        "FirstDate": daily["date"].min().date(),
        "LastDate": daily["date"].max().date(),
        "SampleSize": daily["date"].nunique(),
        "CurveCount": audit["curve_count"],
        "ExampleDate": audit["example_date"],
        "MeanParRMSE": f"{audit['median_par_rmse_bp']:.2f}",
        "MeanZeroRMSE": f"{audit['median_zero_rmse_bp']:.2f}",
        "CurveMinDiscount": f"{audit['min_discount']:.3f}",
        "NegativeForwardDates": audit["negative_forward_dates"],
        "MedianTurns": f"{audit['median_forward_turns']:.0f}",
        "WeakForwardShift": f"{audit['forward_sensitivity_bp']['0.25']:.2f}",
        "StrongForwardShift": f"{audit['forward_sensitivity_bp']['4']:.2f}",
        "FitError": f"{audit['hw_max_discount_error']:.2e}",
        "RiccatiError": f"{audit['riccati_a_error']:.2e}",
        "InitialSigma": f"{100*daily.iloc[0]['sigma']:.2f}",
        "IefDistributionCount": audit["distribution_counts"]["IEF"],
        "TltDistributionCount": audit["distribution_counts"]["TLT"],
    }
    for prefix, strategy in (("Model", "HJB 5d"), ("Matched", "HJB matched duration")):
        row = intervals[(intervals["strategy"] == strategy) & (intervals["block_days"] == 20)].iloc[0]
        for suffix, field in (("Delta", "point_bp"), ("Low", "low_bp"), ("High", "high_bp")):
            macros[prefix + suffix] = f"{row[field]:.2f}"
    (OUT / "numbers.tex").write_text(
        "\n".join("\\newcommand{\\" + key + "}{" + str(value) + "}" for key, value in macros.items()),
        encoding="utf-8",
    )


def main():
    OUT.mkdir(exist_ok=True)
    rates = read_treasury()
    funds = {ticker: read_etf(ticker) for ticker in ("IEF", "TLT")}
    if not funds["IEF"].index.equals(funds["TLT"].index):
        raise ValueError("Unexpected mismatch in ETF close dates")
    market = pd.concat(funds, axis=1)
    diagnostics, example_curves = fit_curves(rates)
    daily, penalties = replay(market, rates, diagnostics)
    summary = summarize(daily, penalties)
    intervals = paired_risk_intervals(daily)
    assert (daily["treasury_date"] <= daily["signal_date"]).all()
    assert (daily["signal_date"] < daily["execution_date"]).all()
    assert (daily["execution_date"] < daily["date"]).all()
    assert np.isfinite(daily.select_dtypes("number")).all().all()
    market.to_csv(OUT / "market.csv")
    daily.to_csv(OUT / "daily.csv", index=False)
    summary.to_csv(OUT / "summary.csv")
    intervals.to_csv(OUT / "paired_risk_intervals.csv", index=False)

    _, zero = example_curves
    sigma = daily.iloc[-1]["sigma"]
    theta_tenors = np.linspace(0, 30, 301)
    discount_error = max(abs(reconstruct_discount(zero, t, sigma) - np.exp(-t * zero(t)))
                         for t in (0.25, 1, 5, 10, 20, 30))
    beta = np.array([-7.0, -12.0])
    test_lambda, test_sigma, test_horizon = 0.0001, 0.02, 0.25
    a_num, _, _ = riccati(beta, test_sigma, test_lambda, test_horizon, 0.)
    alpha = test_sigma * abs(beta[1])
    a_exact = alpha * np.sqrt(test_lambda) * np.tanh(alpha * test_horizon / np.sqrt(test_lambda))
    assert discount_error < 1e-8 and abs(a_num - a_exact) < 1e-8
    diagnostics.to_csv(OUT / "curve_diagnostics.csv")
    pd.DataFrame({
        "maturity": theta_tenors, "zero": zero_value(zero, theta_tenors),
        "discount": np.exp(-theta_tenors * zero_value(zero, theta_tenors)),
        "forward": forward(zero, theta_tenors), "theta": theta(zero, theta_tenors, sigma),
    }).to_csv(OUT / "example_curve.csv", index=False)
    audit = {
        "source_files": {
            name: {"url": url, "sha256": hashlib.sha256((DATA / name).read_bytes()).hexdigest()}
            for name, url in SOURCES.items()
        },
        "curve_count": len(diagnostics),
        "example_date": str(rates.index[-1].date()),
        "median_par_rmse_bp": diagnostics["par_rmse_bp"].median(),
        "median_zero_rmse_bp": diagnostics["zero_rmse_bp"].median(),
        "min_discount": diagnostics["min_discount"].min(),
        "negative_forward_dates": int((diagnostics["min_forward"] < 0).sum()),
        "median_forward_turns": diagnostics["forward_turns"].median(),
        "forward_sensitivity_bp": plot_curves(rates, example_curves),
        "hw_max_discount_error": discount_error,
        "riccati_a_error": abs(a_num - a_exact),
        "distribution_counts": {ticker: int((frame["cash"] > 0).sum())
                                for ticker, frame in funds.items()},
        "max_adjusted_return_gap_bp": {ticker: 1e4 * frame["adjustment_gap"].abs().max()
                                       for ticker, frame in funds.items()},
        "max_close_gap_dollars": {ticker: (frame["close"] - frame["yahoo_close"]).abs().max()
                                  for ticker, frame in funds.items()},
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    plot_strategies(daily, summary)
    write_report_numbers(summary, daily, audit, intervals)
    print(summary[["gross_sd", "net_total", "turnover", "variance_reduction"]].to_string())


if __name__ == "__main__":
    main()
