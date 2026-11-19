# Continuous-time optimal hedging of Treasury ETFs

The report studies whether trading-rate control improves the historical hedge of one IEF share with a short TLT position. The project includes raw Treasury and ETF source snapshots and reproduces the LaTeX tables, five vector figures, and daily strategy ledger.

From this directory, with Python 3.11 or newer, XeLaTeX and the Malgun Gothic font installed, run

```sh
python -m pip install -r requirements.txt
python run.py
python compare_objective.py
python nss_experiment.py
xelatex -interaction=nonstopmode -halt-on-error report.tex
xelatex -interaction=nonstopmode -halt-on-error report.tex
```

Use the included snapshots for an exact data replay. If fresh snapshots are needed, run `python download.py` before the analysis. The analysis records URL and SHA-256 provenance in `results/audit.json`. The Nasdaq sources provide the realized closes and cash distributions. The Yahoo Finance snapshots cross-check dates, closes and adjusted returns.

Run `python compare_objective.py` by itself to regenerate the common-objective table from existing result CSVs.

Run `python nss_experiment.py` with the existing spline outputs and source snapshots to reproduce the NSS appendix. This command writes its outputs to `results_nss/`. The shared controller uses the baseline penalty schedule, and each Original HJB score uses its own model's risk loadings. Use `python nss_experiment.py --from-saved` to refresh the objective and market-performance comparisons from the saved ledgers.

| Output | Contents |
| --- | --- |
| `report.pdf` and `report.tex` | Research report and source |
| `results/summary.csv` | Strategy comparisons and penalty sweep |
| `results/objective_comparison.csv` | Five saved paths scored under the original HJB objective |
| `results/daily.csv` | Signal, execution and P&L dates, positions, costs and variance estimates |
| `results/market.csv` | Aligned closes, cash distributions and adjusted-price checks |
| `results/curve_diagnostics.csv` | Daily fit errors, forward diagnostics and exposure estimates |
| `results/example_curve.csv` | Zero yields, forwards, discounts and Hull–White drift for the report date |
| `results_nss/comparison.csv` | Five strategies under each curve estimator with empirical and objective comparisons |
| `results_nss/curve_fit.csv` | Paired daily fit diagnostics and fitted NSS coefficients |
| `results_nss/signal_theta.csv` | Paired Hull–White drift on the lagged signal dates at selected maturities |
