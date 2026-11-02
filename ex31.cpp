#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
using namespace Eigen;
// input: 19 prices of derivatives (3 LIBOR, 7 futures, 9 swaps) and output the first 14 columns of the cash-flow matrix in Table 3.3.
// quote: "rate for LIBOR, price for futures, rate for swaps"
// date: "date of giving the interest"
int main() {
    // Exercise 3.1: the first 14 columns of the cash-flow matrix in Table 3.3.
    // The dates/row positions are identified in the textbook on p. 37.
    MatrixXd C = MatrixXd::Zero(19, 14);
    double libor[] = {5.59375, 5.625, 5.71875}; // rate of LIBOR in percent
    int liborDays[] = {1, 33, 92}, liborCols[] = {0, 2, 5}; // day to maturity in days, column of cash-flow matrix
    for(int i=0;i<3;i++) C(i,liborCols[i]) = 1 + liborDays[i]/360.0*libor[i]/100.0; // cash-flow of LIBOR

    double futures[] = {94.27,94.26,94.24,94.23,94.18,94.12,94.00};
    int starts[] = {1,3,4,8,9,10,12}; // stating date of the FRA contract: money is taken
    int ends[]   = {6,7,8,9,10,12,13}; // terminal date of the FRA contract: money is given
    for(int i=0;i<7;i++) {
        C(i+3,starts[i]) = -1;
        C(i+3,ends[i]) += 1 + (91.0/360.0)*(1-futures[i]/100.0);
    }
    double swaps[] = {6.01253,6.10823,6.16,6.22,6.32,6.42,6.56,6.56,6.56};
    for(int i=0;i<9;i++) C(i+10,11) = swaps[i]/100.0;
    std::cout << std::fixed << std::setprecision(6) << C << '\n';
}