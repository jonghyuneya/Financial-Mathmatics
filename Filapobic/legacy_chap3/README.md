# Financial Mathematics — Chapter 3 C++

Filipović, *Term-Structure Models*, Chapter 3 연습문제 3.1·3.2·3.5의 C++17 풀이입니다.
외부 수치 라이브러리나 Python 없이 실행합니다. 실제 SNB CSV 데이터도 포함했습니다.

## VS Code에서 실행

1. 이 GitHub 저장소를 Clone한 뒤 VS Code에서 저장소 폴더를 엽니다.
2. Microsoft의 **C/C++** 확장을 설치합니다.
3. C++17 컴파일러를 설치합니다. 확장 설치만으로 컴파일러가 설치되지는 않습니다.
   Windows에서는 MSYS2 UCRT64의 `mingw-w64-ucrt-x86_64-gcc`를 설치하고
   `C:\msys64\ucrt64\bin`을 PATH에 추가한 후 VS Code를 다시 실행합니다.
   `g++ --version`으로 설치를 확인합니다.
4. `main.cpp`를 열고 **Ctrl+Shift+B**로 빌드합니다.
5. 메뉴 **Terminal → Run Task → Run C++ Chapter 3**로 실행합니다.

작업 설정은 저장소 루트 폴더를 VS Code로 열었을 때 적용됩니다.
현재 기본 작업은 GCC/g++용입니다. macOS의 clang++은 아래 터미널 명령을 쓰거나 작업의 command를 변경하세요.

### 터미널로 실행

Windows PowerShell (MinGW GCC):

```powershell
g++ -std=c++17 -O2 -Wall -Wextra main.cpp -o ch3.exe
.\ch3.exe
```

Linux:

```bash
g++ -std=c++17 -O2 -Wall -Wextra main.cpp -o ch3
./ch3
```

macOS (Xcode Command Line Tools):

```bash
clang++ -std=c++17 -O2 -Wall -Wextra main.cpp -o ch3
./ch3
```

MSVC 개발자 터미널:

```powershell
cl /std:c++17 /EHsc /utf-8 main.cpp /Fe:ch3.exe
.\ch3.exe
```

선택적으로 CMake도 지원합니다: `cmake -S . -B build`, `cmake --build build`.
실행 시 작업 디렉터리는 저장소 루트로 두세요. 실행파일 뒤 인수는 입력 폴더, 출력 폴더, PCA 마지막 월입니다.
예: `./ch3 data results 2026-09` (Windows는 `.\ch3.exe data results 2026-09`).

## 문제와 구현

### 3.1: 현금흐름 행렬

문제 요약: §3.2.2의 1997-10-06 미국 단기금융시장 자료로 현금흐름 행렬 처음 14열을 구하세요.

`data/money_market_3_1.csv`에서 표 3.3의 19개 상품을 읽습니다.
기준일 1997-10-08, 액면 1, LIBOR Actual/360, 선물 기간 91/360,
연간 스왑 δ=1을 적용합니다. LIBOR는 만기 1+δL, 선물은 시작 −1·종료 1+δF,
스왑은 매년 K·마지막 1+K입니다. 선물가격은 100(1−F)입니다.
휴일 조정은 추가하지 않습니다.

답: `results/ex3_1_first14.csv`. 전체 행렬은 19×43이며 `ex3_1_cashflows_full.csv`에 저장합니다.
가격 벡터는 `ex3_1_prices.csv`입니다. 책 37쪽의 19×47 표기는 같은 쪽 N=43과 불일치합니다.

### 3.2: 스무딩 스플라인

문제 요약: 만기 2,3,4,5,7,10년의 주어진 수익률에 대해 α=0.01,0.1,1로
식 (3.5)의 선도곡선 계수 f(0), a1,…,a6를 구하고 선도·수익률 곡선을 그리세요.

`data/yields_3_2.csv`를 읽습니다. 금리 입력과 출력은 % 단위입니다.
m=min(u,Tk), hk(u)=Tk+Tk·m−m²/2를 사용하여
f(u)=f0+Σak hk(u), R(u)=∫₀ᵘf(s)ds/u를 계산합니다.
식 (3.7)–(3.8)의 블록 연립방정식에 부분 피벗 가우스 소거법을 적용합니다.

답: `results/ex3_2_coefficients.csv`. α별 곡선 CSV, 관측 만기 적합값 CSV,
RMSE·거칠기·목적함수 파일도 생성합니다.
`results/ex3_2_curves_alpha_*.svg`를 브라우저로 열면 선도·수익률 그래프를 볼 수 있습니다.
검은 점은 관측 수익률입니다. α별 비교에는 동일한 만기 축을 쓰되 세로축 범위는 각 파일에 표시된 값을 확인하세요.

### 3.5: 수익률 월별 변화의 PCA

문제 요약: 2,3,4,5,7,10,20,30년 스위스 연방채 현물금리의 최근 10년 월별 변화에 대해
공분산, 내림차순 고유값·고유벡터, 설명분산을 구하세요.

`data/snb_yields.csv`를 읽으며 형식은 `date,2,3,4,5,7,10,20,30`입니다.
입력 월은 YYYY-MM, 금리는 %입니다. 최근 구간의 누락 월·누락 값·중복 월은 오류로 처리합니다.
코드는 정규화한 CSV 입력을 지원합니다. Excel을 사용할 경우 먼저 위 형식의 CSV로 저장하세요.

교재 옛 주소 `https://www.snb.ch/ext/stats/statmon/xls/en/statmon_E3_M1_M.xls`는
2026-10-31 확인 시 404였으므로 현행 [SNB 공식 데이터](https://data.snb.ch/en/topics/ziredev/cube/rendeiduebm)를 사용했습니다.
[원본 CSV API](https://data.snb.ch/api/cube/rendeiduebm/data/csv/en)에서 D0=CHF,
D1=2J,3J,4J,5J,7J,10J,20J,30J를 추출한 자료입니다.
이 파일은 실행 때 다시 다운로드하지 않는 고정 데이터입니다. 옛 Excel의 당시 빈티지와 동일함을 보장하지 않습니다.
원본 게시일 2026-10-01, 다운로드 확인일 2026-10-31, 마지막 관측월 2026-09입니다.

2016-09~2026-09의 121개 수준으로 2016-10~2026-09의 120개 차분을 계산합니다.
책 §3.4.2에 맞춰 공분산 분모는 N=120입니다. 상관행렬 표준화는 하지 않습니다.
고유값·공분산 단위는 (%포인트)²입니다.
8×8 대칭행렬은 Jacobi 회전법으로 분해합니다. 고유벡터는 열별이고 가장 큰 절댓값 성분의 부호를 양수로 통일합니다.

답 파일:

- `results/ex3_5_covariance.csv`: 공분산
- `results/ex3_5_eigenvectors.csv`: 고유벡터
- `results/ex3_5_variance.csv`: 내림차순 고유값·설명분산
- `results/ex3_5_monthly_changes.csv`, `ex3_5_mean_changes.csv`, `ex3_5_scores.csv`: 중간 계산

PC1=88.4477%, PC2=9.3831%, PC3=1.7414%, 첫 세 성분 합계=99.5721%입니다.

3.3·3.4는 증명 문제이므로 이 프로그램의 계산 대상에서 제외했습니다.

## 검증

GCC에서 `-std=c++17 -O2 -Wall -Wextra -pedantic`으로 경고 없이 빌드하고 실행했습니다.
3.1 전체 행렬, 3.2 계수·오차, 3.5 공분산·고유값·고유벡터를
기존 NumPy 풀이와 대조하여 허용 오차 내 일치를 확인했습니다.
Windows·macOS에서는 직접 실행 검증하지 않았습니다.
