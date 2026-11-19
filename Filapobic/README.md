# Filapobic

## Chapter 3

`chap3/`에는 연습문제 3.1, 3.2, 3.5의 C++ 소스와 연습문제 3.2의 결과 CSV가 있습니다. C++17과 Eigen이 필요합니다. 아래 예시는 Eigen 헤더가 `/usr/include/eigen3`에 있는 Linux 또는 WSL 환경에서 실행합니다.

```bash
cd Filapobic/chap3
g++ -std=c++17 -O2 -I /usr/include/eigen3 ex31.cpp -o ex31
g++ -std=c++17 -O2 -I /usr/include/eigen3 ex32.cpp -o ex32
g++ -std=c++17 -O2 -I /usr/include/eigen3 ex35.cpp -o ex35
./ex35
```

연습문제 3.5의 기본 입력은 `../legacy_chap3/data/snb_yields.csv`이며, 다른 CSV 경로는 첫 번째 인수로 지정할 수 있습니다. 포함된 rapidcsv 헤더는 [BSD 3-Clause 라이선스](chap3/LICENSE.rapidcsv)를 따릅니다.

## 기존 Chapter 3 구현

기존 루트의 C++ 코드, 입력 CSV, CMake, VS Code 설정은 [legacy_chap3/](legacy_chap3/)에 모았습니다. 이전 실행 절차는 [기존 README](legacy_chap3/README.md)에 있습니다.
