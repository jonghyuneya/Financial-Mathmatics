"""Compare the archived spline replay with a direct NSS fit to Treasury par quotes."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import least_squares

from compare_objective import STRATEGIES
from run import (
    A_RATE, BASE_DAYS, OUT, PROXY_MATURITIES, TENORS, curve,
    forward, loading, plt, read_treasury, replay, summarize, theta, zero_value,
)

OUT_NSS = Path(__file__).resolve().parent / "results_nss"
OBSERVED_TENORS = np.array(list(TENORS.values()))
COUPON_GRID = np.arange(0.5, 30.5, 0.5)
PAR_GRID = np.r_[OBSERVED_TENORS[OBSERVED_TENORS < 0.5], COUPON_GRID]
COUPON_INDICES = np.rint(2 * OBSERVED_TENORS[OBSERVED_TENORS >= 0.5]).astype(int) - 1
SHORT_COUNT = len(PAR_GRID) - len(COUPON_GRID)
DECAY_STARTS = [(0.3, 2), (0.7, 7), (1.5, 0.5), (2, 8), (4, 16), (8, 2)]


def implied_par(z, maturities):
    """Simple bill yields followed by a full semiannual coupon schedule."""
    log_discounts = -maturities * z(maturities)
    short = maturities < 0.5
    discounts = np.exp(log_discounts[~short])
    return np.r_[np.expm1(-log_discounts[short]) / maturities[short],
                 -2 * np.expm1(log_discounts[~short]) / np.cumsum(discounts)]


def quoted_par(z):
    rates = implied_par(z, PAR_GRID)
    return rates[np.r_[np.arange(SHORT_COUNT), SHORT_COUNT + COUPON_INDICES]]


def nss_zero(t, parameters):
    b0, b1, b2, b3, log_tau1, log_tau2 = parameters
    x = np.asarray(t) / np.exp(log_tau1)
    y = np.asarray(t) / np.exp(log_tau2)
    l1 = np.divide(-np.expm1(-x), x, out=np.ones_like(x), where=x != 0)
    l2 = np.divide(-np.expm1(-y), y, out=np.ones_like(y), where=y != 0)
    return b0 + b1 * l1 + b2 * (l1 - np.exp(-x)) + b3 * (l2 - np.exp(-y))


def nss_forward_and_slope(t, parameters):
    b0, b1, b2, b3, log_tau1, log_tau2 = parameters
    tau1, tau2 = np.exp([log_tau1, log_tau2])
    x, y = np.asarray(t) / tau1, np.asarray(t) / tau2
    e1, e2 = np.exp(-x), np.exp(-y)
    f = b0 + (b1 + b2 * x) * e1 + b3 * y * e2
    slope = (b2 * (1 - x) - b1) * e1 / tau1 + b3 * (1 - y) * e2 / tau2
    return f, slope


def nss_theta(t, parameters, sigma):
    f, slope = nss_forward_and_slope(t, parameters)
    return slope + A_RATE * f + sigma**2 * loading(2 * t) / 2


def fit_nss(yields):
    def residual(parameters):
        return 1e4 * (quoted_par(lambda t: nss_zero(t, parameters)) - yields)

    lower = np.r_[np.full(4, -0.30), np.log([0.1, 0.1])]
    upper = np.r_[np.full(4, 0.30), np.log([30, 30])]
    fits = []
    for taus in DECAY_STARTS:
        start = np.r_[yields[-1], yields[0] - yields[-1], 0., 0., np.log(taus)]
        fit = least_squares(residual, start, bounds=(lower, upper), max_nfev=160,
                            ftol=1e-7, xtol=1e-7, gtol=1e-7)
        if fit.success:
            fits.append(fit)
    if not fits:
        raise RuntimeError("All NSS optimization starts failed")
    best = min(fits, key=lambda f: f.cost)
    at_bound = np.any(np.isclose(best.x, lower, atol=1e-5) |
                      np.isclose(best.x, upper, atol=1e-5))
    return best.x, len(fits), bool(at_bound)


def fit_curves(rates, daily):
    signals = daily[daily["strategy"] == "No hedge"]
    signals_by_curve = {date: rows for date, rows in signals.groupby("treasury_date")}
    grid = np.linspace(0, 30, 601)
    drift_tenors = np.array([0., .25, 1., 5., 10., 20., 30.])
    fit_rows, diagnostics, drift_rows = [], [], []
    sigma_example = signals["sigma"].iloc[-1]
    for date, row in rates.iterrows():
        observed = row.to_numpy()
        parameters, starts, at_bound = fit_nss(observed)
        _, spline, nodes, targets = curve(OBSERVED_TENORS, observed)

        def z_nss(t):
            return nss_zero(t, parameters)

        def z_spline(t):
            return zero_value(spline, t)

        nss_fwd, _ = nss_forward_and_slope(grid, parameters)
        diagnostics_row = {"date": date}
        par_nss = implied_par(z_nss, COUPON_GRID)
        for ticker, maturity in zip(("IEF", "TLT"), PROXY_MATURITIES):
            count = int(2 * maturity)
            payments = COUPON_GRID[:count]
            cashflows = np.full(count, par_nss[count - 1] / 2)
            cashflows[-1] += 1
            pv = cashflows * np.exp(-payments * z_nss(payments))
            diagnostics_row[f"{ticker}_duration"] = np.average(payments, weights=pv)
            diagnostics_row[f"{ticker}_hw"] = np.average(loading(payments), weights=pv)
        diagnostics.append(diagnostics_row)
        fit_row = {"date": date, "successful_starts": starts, "at_bound": at_bound}
        for label, z, fwd in (("spline", z_spline, forward(spline, grid)),
                              ("nss", z_nss, nss_fwd)):
            fit_row[f"{label}_par_rmse_bp"] = 1e4 * np.sqrt(np.mean((quoted_par(z) - observed)**2))
            fit_row[f"{label}_zero_rmse_bp"] = 1e4 * np.sqrt(np.mean((z(nodes) - targets)**2))
            fit_row[f"{label}_min_forward"] = np.min(fwd)
        fit_row.update(zip(("beta0", "beta1", "beta2", "beta3", "tau1", "tau2"),
                           np.r_[parameters[:4], np.exp(parameters[4:])]))
        fit_rows.append(fit_row)
        for signal in signals_by_curve.get(date, pd.DataFrame()).itertuples(index=False):
            sigma = signal.sigma
            for t, nss, old in zip(drift_tenors, nss_theta(drift_tenors, parameters, sigma),
                                   theta(spline, drift_tenors, sigma)):
                drift_rows.append({"signal_date": signal.signal_date, "treasury_date": date,
                                   "maturity": t, "sigma": sigma,
                                   "theta_nss": nss, "theta_spline": old})
        if date == rates.index[-1]:
            example = pd.DataFrame({
                "maturity": grid, "nss_zero": z_nss(grid), "spline_zero": z_spline(grid),
                "nss_discount": np.exp(-grid * z_nss(grid)),
                "spline_discount": np.exp(-grid * z_spline(grid)),
                "nss_forward": nss_fwd, "spline_forward": forward(spline, grid),
                "nss_theta": nss_theta(grid, parameters, sigma_example),
                "spline_theta": theta(spline, grid, sigma_example),
            })
            example_fit = quoted_par(z_nss), quoted_par(z_spline)
            t = np.array([.25, 1., 5., 10., 20., 30.])
            h = 1e-4
            estimated = ((t + h) * z_nss(t + h) - (t - h) * z_nss(t - h)) / (2*h)
            np.testing.assert_allclose(estimated, nss_forward_and_slope(t, parameters)[0], atol=1e-9)
            errors = []
            for maturity in t:
                mean = loading(maturity) * nss_forward_and_slope(0., parameters)[0]
                mean += quad(lambda s: loading(maturity - s) * nss_theta(s, parameters, sigma_example),
                             0, maturity)[0]
                variance = sigma_example**2 * quad(lambda s: loading(s)**2, 0, maturity)[0]
                errors.append(abs(np.exp(-mean + variance/2) - np.exp(-maturity*z_nss(maturity))))
            assert max(errors) < 1e-8
    return (pd.DataFrame(diagnostics).set_index("date"), pd.DataFrame(fit_rows),
            pd.DataFrame(drift_rows), example, example_fit, max(errors))


def compare(old, new, baseline_objective):
    metrics = []
    for method, daily in (("Spline", old), ("NSS", new)):
        for key, name in STRATEGIES.items():
            path = daily[daily["strategy"] == key].sort_values("date")
            row = {
                "method": method, "strategy": name, "days": len(path),
                "gross_variance": path["gross"].var(), "net_variance": path["net"].var(),
                "turnover": path["turnover"].sum(), "net_pnl": path["net"].sum(),
                "trading_cost": (path["gross"] - path["net"]).sum(),
            }
            if key == f"HJB {BASE_DAYS}d":
                baseline = baseline_objective.loc[name]
                # The ledger stores variance contributions under its own curve.
                risk = baseline["gamma"] * path["hw_variance"].sum()
                rate = baseline["lambda"] * path["rate_energy"].sum()
                terminal = baseline["eta"] * path["position"].iloc[-1]**2
                row.update(risk_loss_native=risk, trading_rate_penalty=rate,
                           terminal_penalty=terminal, total_objective_native=risk + rate + terminal)
                row.update(baseline[["gamma", "lambda", "eta"]].to_dict())
                if method == "Spline":
                    np.testing.assert_allclose(
                        [risk, rate, terminal],
                        baseline[["risk_loss", "trading_rate_penalty", "terminal_penalty"]]
                        .to_numpy(dtype=float), atol=1e-13,
                    )
            metrics.append(row)
    return pd.DataFrame(metrics)


def save_figure(fig, name):
    fig.tight_layout()
    fig.savefig(OUT_NSS / name, bbox_inches="tight")
    plt.close(fig)


def plot_curves(example, example_fit, observed):
    fig, axes = plt.subplots(2, 2, figsize=(9, 5))
    par_ax, fwd_ax, discount_ax, theta_ax = axes.ravel()
    par_ax.scatter(OBSERVED_TENORS, 100 * observed, label="Observed", color="black", s=12)
    for method, prefix, fitted in zip(("NSS", "Spline"), ("nss", "spline"), example_fit):
        par_ax.plot(OBSERVED_TENORS, 100 * fitted, label=method)
        for ax, field, scale in ((fwd_ax, "forward", 100), (discount_ax, "discount", 1),
                                 (theta_ax, "theta", 100)):
            ax.plot(example["maturity"], scale * example[f"{prefix}_{field}"], label=method)
    for ax, label in zip(axes.ravel(), ("Implied par yield (%)", "Forward rate (%)",
                                        "Discount factor", "Hull–White drift (%/year)")):
        ax.set(xlabel="Maturity (years)", ylabel=label)
    par_ax.legend()
    save_figure(fig, "curves.pdf")


def plot_performance(comparison):
    fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    x = np.arange(len(STRATEGIES))
    fields = (("gross_variance", "Gross variance (bp²)", 1e8),
              ("turnover", "Turnover (% of initial IEF value)", 100),
              ("trading_cost", "Trading cost (bp)", 1e4))
    for ax, (field, title, scale) in zip(axes.ravel(), fields):
        for offset, method in ((-.19, "Spline"), (.19, "NSS")):
            values = comparison[comparison["method"] == method][field]
            ax.bar(x + offset, values * scale, .38, label=method)
        ax.set_title(title, fontsize=10)
    axes[0].legend()
    axes[-1].set_xticks(x, ["No hedge", "Static", "Dynamic", "Original\nHJB", "Duration-\ntarget HJB"])
    save_figure(fig, "metrics.pdf")


def write_tables(comparison, fits):
    fields = (("gross_variance", 1e8), ("turnover", 100), ("trading_cost", 1e4))
    rows = []
    for strategy, group in comparison.groupby("strategy", sort=False):
        methods = group.set_index("method")
        cells = [f"{methods.loc[method, field] * scale:.3f}"
                 for field, scale in fields for method in ("Spline", "NSS")]
        rows.append(" & ".join([strategy, *cells]) + r" \\")
    (OUT_NSS / "performance.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")
    hjb = comparison[comparison["strategy"] == "Original HJB"].set_index("method")
    fields = ("risk_loss_native", "trading_rate_penalty", "terminal_penalty", "total_objective_native")
    rows = []
    for method, label in (("Spline", "Smoothing Splines HJB"), ("NSS", "NSS HJB")):
        cells = [f"{hjb.loc[method, field] * 1e6:.3f}" for field in fields]
        rows.append(" & ".join([label, *cells]) + r" \\")
    (OUT_NSS / "objective.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")
    stats = {
        "NssSplineParMedian": fits["spline_par_rmse_bp"].median(),
        "NssParMedian": fits["nss_par_rmse_bp"].median(),
        "NssSplineZeroMedian": fits["spline_zero_rmse_bp"].median(),
        "NssZeroMedian": fits["nss_zero_rmse_bp"].median(),
        "NssBoundDates": int(fits["at_bound"].sum()),
    }
    for name, method in (("NssHjbSpline", "Spline"), ("NssHjbNss", "NSS")):
        row = hjb.loc[method]
        stats[f"{name}Sd"] = 1e4 * np.sqrt(row["gross_variance"])
        stats[f"{name}Turnover"] = 100 * row["turnover"]
    macros = [
        fr"\newcommand{{\{name}}}{{{value if isinstance(value, int) else f'{value:.3f}'}}}"
        for name, value in stats.items()
    ]
    (OUT_NSS / "numbers.tex").write_text("\n".join(macros) + "\n", encoding="utf-8")


def update_comparison(old, new, objective, fits):
    comparison = compare(old, new, objective)
    comparison.to_csv(OUT_NSS / "comparison.csv", index=False)
    write_tables(comparison, fits)
    plot_performance(comparison)
    return comparison


def main():
    OUT_NSS.mkdir(exist_ok=True)
    rates = read_treasury()
    market = pd.read_csv(OUT / "market.csv", index_col=0, header=[0, 1], parse_dates=True)
    market.columns = market.columns.set_names([None, None])
    date_columns = ["date", "execution_date", "signal_date", "treasury_date"]
    old = pd.read_csv(OUT / "daily.csv", parse_dates=date_columns)
    baseline_diagnostics = pd.read_csv(OUT / "curve_diagnostics.csv", index_col=0, parse_dates=True)
    objective = pd.read_csv(OUT / "objective_comparison.csv", index_col="strategy")
    diagnostics, fits, drift, example, example_fit, discount_error = fit_curves(rates, old)
    new, penalties = replay(market, rates, diagnostics, penalty_reference=baseline_diagnostics)
    assert np.isclose(penalties[BASE_DAYS], objective["lambda"].iloc[0], rtol=1e-12, atol=0)
    pd.testing.assert_frame_equal(old[date_columns + ["strategy", "sigma"]],
                                 new[date_columns + ["strategy", "sigma"]])
    np.testing.assert_allclose(old.loc[old["strategy"] == "No hedge", "gross"],
                               new.loc[new["strategy"] == "No hedge", "gross"], atol=1e-14)
    nss_summary = summarize(new, penalties)
    for name, frame in (("daily", new), ("curve_fit", fits),
                        ("signal_theta", drift), ("example_curve", example)):
        frame.to_csv(OUT_NSS / f"{name}.csv", index=False)
    diagnostics.to_csv(OUT_NSS / "curve_diagnostics.csv")
    nss_summary.to_csv(OUT_NSS / "summary.csv")
    audit = {
        "example_discount_reconstruction_max_error": discount_error,
        "curve_dates": len(rates),
        "evaluation_dates": int(new["date"].nunique()),
        "baseline_files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in sorted(OUT.glob("*.csv"))},
    }
    (OUT_NSS / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    comparison = update_comparison(old, new, objective, fits)
    plot_curves(example, example_fit, rates.iloc[-1].to_numpy())
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-saved", action="store_true",
                        help="Update scores, tables and charts from saved ledgers.")
    if parser.parse_args().from_saved:
        update_comparison(
            pd.read_csv(OUT / "daily.csv"), pd.read_csv(OUT_NSS / "daily.csv"),
            pd.read_csv(OUT / "objective_comparison.csv", index_col="strategy"),
            pd.read_csv(OUT_NSS / "curve_fit.csv"),
        )
    else:
        main()
