"""Download public source snapshots for the 2022–2025 sample."""
from pathlib import Path
from urllib.request import Request, urlopen

DATA = Path(__file__).resolve().parent / "data"
TREASURY = "https://home.treasury.gov/resource-center/data-chart-center/interest-rates"
SOURCES = {
    **{
        f"treasury_{year}.csv": (
            f"{TREASURY}/daily-treasury-rates.csv/{year}/all"
            f"?type=daily_treasury_yield_curve&field_tdr_date_value={year}&page&_format=csv"
        )
        for year in range(2022, 2026)
    },
    **{
        f"{ticker}.json": (
            f"https://query2.finance.yahoo.com/v8/finance/chart/{ticker}"
            "?period1=1640995200&period2=1767225600&interval=1d&events=div%2Csplits"
        )
        for ticker in ("IEF", "TLT")
    },
    **{
        f"{ticker}-{kind}.json": f"https://api.nasdaq.com/api/quote/{ticker}/{endpoint}"
        for ticker in ("IEF", "TLT")
        for kind, endpoint in (
            ("nasdaq", "historical?assetclass=etf&fromdate=2022-01-01&todate=2025-12-31&limit=1500"),
            ("dividends", "dividends?assetclass=etf"),
        )
    },
}


if __name__ == "__main__":
    DATA.mkdir(exist_ok=True)
    for name, url in SOURCES.items():
        agent = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                 "AppleWebKit/537.36 Chrome/130.0.0.0 Safari/537.36")
        with urlopen(Request(url, headers={"User-Agent": agent}), timeout=30) as response:
            content = response.read()
        temporary = DATA / f"{name}.tmp"
        temporary.write_bytes(content)
        temporary.replace(DATA / name)
        print(name, len(content))
