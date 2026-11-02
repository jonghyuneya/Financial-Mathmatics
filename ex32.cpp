#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace Eigen;
int main() {
    const int n=6;
    VectorXd T(n),Y(n);
    T << 2,3,4,5,7,10;
    Y << .948,1.129,1.300,1.454,1.704,1.955; // percent
    for(double alpha : {0.01,0.1,1.0}) {
        MatrixXd A=MatrixXd::Zero(n+1,n+1);
        VectorXd b=VectorXd::Zero(n+1);
        for(int i=0;i<n;i++) {
            A(0,i+1)=T(i); // for sum(a_i*T_i)
            A(i+1,0)=alpha*T(i);
            b(i+1)=alpha*T(i)*Y(i); // for alpha*Y_i*T_i=alpha*(sum(a_l<h_l,h_i>)+f(0)*T_i)
            for(int j=0;j<n;j++) {
                double m=std::min(T(i),T(j));
                double inner=T(i)*T(j)*(1+m) -(T(i)+T(j))*m*m/2 + m*m*m/3;
                A(i+1,j+1)=alpha*inner+(i==j ? 1.0 : 0.0);
            }
        }
        VectorXd x=A.fullPivLu().solve(b); // f(0), a_1,...,a_6
        std::cout << "alpha=" << alpha << "\n" << x.transpose() << "\n";
        std::string label=(alpha==.01?"001":alpha==.1?"01":"1");
        std::ofstream out("ex32_alpha_"+label+".csv");
        out << "maturity,forward,yield\n" << std::setprecision(12);
        for(int step=0;step<=500;step++) {
            double u=10.0*step/500.0, f=x(0), integrated=x(0)*u;
            for(int i=0;i<n;i++) {
                double t=T(i),m=std::min(u,t);
                // h_i(0)=T_i, h_i'(u)=(T_i-u)_+.
                f += x(i+1)*(t+t*m-m*m/2); // = a_i*h_i(u)
                // Integral from 0 to u of h_i(s) ds, cumulated sum for forward rate
                integrated += x(i+1)*(t*u+t*m*m/2-m*m*m/6
                                     +(u-m)*(t*m-m*m/2));
                //calculating yield as a integral of forward rate divided by time to maturity
            }
            out << u << ',' << f << ',' << (u==0?f:integrated/u) << '\n';
        }
    }
}