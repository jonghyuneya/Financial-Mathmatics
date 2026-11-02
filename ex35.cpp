#include <Eigen/Dense>
#include "rapidcsv.h"
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace Eigen;
int main() {
    // CSV: date,2,3,4,5,7,10,20,30; oldest -> newest; one row/month.
    rapidcsv::Document csv("data/snb_yields.csv");
    std::vector<std::string> tenors={"2","3","4","5","7","10","20","30"};
    const int n=8;
    int total=(int)csv.GetRowCount();
    if(total<121) throw std::runtime_error("Need at least 121 monthly yield observations");
    MatrixXd yields(121,n);

    for(int j=0;j<n;j++) {
    auto col = csv.GetColumn<std::string>(tenors[j]);

    for(int i=0;i<121;i++) {
        int row = total-121+i;

        if(col[row].empty()) {
            std::cerr << "Missing value: "
                      << "row=" << row+2
                      << ", tenor=" << tenors[j] << '\n';
            return 1;
        }

        yields(i,j) = std::stod(col[row]);
    }
    }
    MatrixXd changes=yields.bottomRows(120)-yields.topRows(120);
    RowVectorXd mean=changes.colwise().mean();
    MatrixXd centered=changes.rowwise()-mean;
    MatrixXd Q=centered.transpose()*centered/120.0; // empirical definition in Sect. 3.4
    SelfAdjointEigenSolver<MatrixXd> es(Q);
    if(es.info()!=Success) throw std::runtime_error("Eigen decomposition failed");
    std::cout << "Covariance Q:\n" << Q << "\n\n";
    double totalVariance=es.eigenvalues().sum(), cumulative=0;
    std::cout << "PC  eigenvalue  explained(%)  cumulative(%)\n";
    for(int j=n-1;j>=0;j--) { // Eigen returns ascending order
        double p=100*es.eigenvalues()(j)/totalVariance;
        cumulative+=p;
        std::cout << n-j << "  " << es.eigenvalues()(j) << "  " << p << "  " << cumulative << '\n';
    }
    std::cout << "Eigenvectors (columns: PC1 ... PC8):\n"
              << es.eigenvectors().rowwise().reverse() << '\n';
}