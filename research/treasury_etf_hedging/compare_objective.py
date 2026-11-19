"""Score the saved backtest paths under the original Hull–White objective."""
import csv
from math import fsum, isclose, sqrt
from pathlib import Path

OUT = Path(__file__).resolve().parent / "results"
STRATEGIES = {
    "No hedge": "No Hedge",
    "Static duration": "Static Duration",
    "Dynamic duration": "Dynamic Duration",
    "HJB 5d": "Original HJB",
    "HJB matched duration": "Duration-Target HJB",
}


def read_rows(name):
    with (OUT / name).open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def main():
    summary = {row["strategy"]: row for row in read_rows("summary.csv")}
    penalty = float(summary["HJB 5d"]["lambda"])
    daily = read_rows("daily.csv")
    unhedged = sorted((row for row in daily if row["strategy"] == "No hedge"),
                      key=lambda row: row["date"])
    last = unhedged[-1]
    curves = {row["date"]: row for row in read_rows("curve_diagnostics.csv")}
    with (OUT / "market.csv").open(newline="", encoding="utf-8") as source:
        reader = csv.reader(source)
        headers = list(zip(next(reader), next(reader)))
        ief = headers.index(("IEF", "close"))
        tlt = headers.index(("TLT", "close"))
        prices = {row[0]: (float(row[ief]), float(row[tlt])) for row in reader}
    reference = prices[unhedged[0]["signal_date"]][0]
    final_loading = float(curves[last["treasury_date"]]["TLT_hw"])
    beta_h = prices[last["signal_date"]][1] * final_loading / reference
    eta = sqrt(penalty) * float(last["sigma"]) * beta_h

    rows = []
    table = []
    for strategy, label in STRATEGIES.items():
        path = sorted((row for row in daily if row["strategy"] == strategy),
                      key=lambda row: row["date"])
        if [row["date"] for row in path] != [row["date"] for row in unhedged]:
            raise ValueError(f"Unaligned evaluation dates for {strategy}")
        risk = fsum(float(row["hw_variance"]) for row in path)
        trading = penalty * fsum(float(row["rate_energy"]) for row in path)
        assert isclose(risk, float(summary[strategy]["hw_variance"]))
        assert isclose(trading, float(summary[strategy]["rate_penalty_common"]))
        q_end = float(path[-1]["position"])
        terminal = eta * q_end**2
        total = fsum((risk, trading, terminal))
        rows.append([label, risk, trading, terminal, total, q_end, 1, penalty, eta])
        table.append(" & ".join([label, *(f"{1e6 * x:.3f}" for x in
                                           (risk, trading, terminal, total))]) + r" \\")

    with (OUT / "objective_comparison.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["strategy", "risk_loss", "trading_rate_penalty", "terminal_penalty",
                         "total_objective", "terminal_position", "gamma", "lambda", "eta"])
        writer.writerows(rows)
    (OUT / "objective_rows.tex").write_text("\n".join(table) + "\n", encoding="utf-8")
    parameters = []
    for name, value in (("ObjectiveLambda", penalty), ("ObjectiveEta", eta)):
        mantissa, exponent = f"{value:.6e}".split("e")
        parameters.append(fr"\newcommand{{\{name}}}{{{mantissa}\times10^{{{int(exponent)}}}}}")
    (OUT / "objective_parameters.tex").write_text("\n".join(parameters) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
