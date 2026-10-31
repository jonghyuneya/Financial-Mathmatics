// C++17, standard library only. Rates are in percent; covariance in percentage-points squared.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using Vec = std::vector<double>;
using Mat = std::vector<Vec>;
namespace fs = std::filesystem;
using Rows = std::vector<std::vector<std::string>>;
Rows readCSV(const fs::path& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot read: " + path.string());
    Rows rows; std::string line;
    while (std::getline(f,line)) {
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> row; std::string cell; std::istringstream s(line);
        while (std::getline(s,cell,',')) row.push_back(cell);
        rows.push_back(row);
    }
    if(rows.empty()) throw std::runtime_error("Empty CSV");
    return rows;
}
double number(const std::string& s) {
    size_t used=0; double v=std::stod(s,&used);
    if(used!=s.size() || !std::isfinite(v)) throw std::runtime_error("Invalid number: "+s);
    return v;
}
std::ofstream output(const fs::path& p) {
    std::ofstream f(p); if(!f) throw std::runtime_error("Cannot write: "+p.string());
    f << std::setprecision(16); return f;
}
void writeMat(const fs::path& p,const Mat& a,const std::vector<std::string>& labels) {
    auto f=output(p); f<<"maturity"; for(auto& s:labels) f<<','<<s; f<<'\n';
    const Vec tenors={2,3,4,5,7,10,20,30};
    for(size_t i=0;i<a.size();++i){ f<<tenors[i]; for(double x:a[i]) f<<','<<x; f<<'\n'; }
}
Vec solve(Mat a,Vec b) {
    const size_t n=b.size();
    for(size_t k=0;k<n;++k){
        size_t pivot=k; for(size_t i=k+1;i<n;++i) if(std::abs(a[i][k])>std::abs(a[pivot][k])) pivot=i;
        if(std::abs(a[pivot][k])<1e-14) throw std::runtime_error("Singular system");
        std::swap(a[k],a[pivot]); std::swap(b[k],b[pivot]);
        for(size_t i=k+1;i<n;++i){ double q=a[i][k]/a[k][k];
            for(size_t j=k;j<n;++j) { a[i][j]-=q*a[k][j]; }
            b[i]-=q*b[k]; }
    }
    Vec x(n); for(int i=int(n)-1;i>=0;--i){ double v=b[i]; for(size_t j=i+1;j<n;++j) v-=a[i][j]*x[j]; x[i]=v/a[i][i]; }
    return x;
}
bool leap(int y){return y%4==0 && (y%100!=0 || y%400==0);}
int ordinal(const std::string& s){
    int y=std::stoi(s.substr(0,4)),m=std::stoi(s.substr(5,2)),d=std::stoi(s.substr(8,2));
    int months[]={31,28,31,30,31,30,31,31,30,31,30,31};
    if(y<1900 || y>2200 || m<1 || m>12 || d<1 || d>months[m-1]+(m==2 && leap(y))) throw std::runtime_error("Invalid date");
    int n=0; for(int k=1900;k<y;++k)n+=365+leap(k);
    for(int k=1;k<m;++k) { n+=months[k-1]+(k==2 && leap(y)); }
    return n+d-1;
}
std::string date(int n){
    int y=1900,m=1; while(n>=365+leap(y)){n-=365+leap(y);++y;}
    int months[]={31,28,31,30,31,30,31,31,30,31,30,31};
    while(n>=months[m-1]+(m==2 && leap(y))){n-=months[m-1]+(m==2 && leap(y));++m;}
    std::ostringstream s; s<<y<<'-'<<std::setfill('0')<<std::setw(2)<<m<<'-'<<std::setw(2)<<n+1; return s.str();
}
void exercise31(const fs::path& data,const fs::path& out){
    auto rows=readCSV(data/"money_market_3_1.csv"); int spot=ordinal("1997-10-08");
    std::vector<std::map<int,double>> flows; std::vector<std::string> names; Vec prices; std::map<int,int> dates;
    for(size_t i=1;i<rows.size();++i){
        auto& r=rows[i]; if(r.size()!=4) throw std::runtime_error("Invalid money market row");
        double q=number(r[3]); std::map<int,double> flow; double price=1;
        if(r[0]=="libor"){int t=ordinal(r[2]);flow[t]=1+(t-spot)/360.0*q/100;}
        else if(r[0]=="future"){int t=ordinal(r[2]);flow[t]=-1;flow[t+91]=1+91.0/360*(1-q/100);price=0;}
        else if(r[0]=="swap"){int years=std::stoi(r[1]);for(int k=1;k<=years;++k)flow[ordinal(std::to_string(1997+k)+"-10-08")]=q/100;flow.rbegin()->second+=1;}
        else throw std::runtime_error("Unknown instrument");
        for(auto& x:flow)dates[x.first]=0;
        flows.push_back(flow);prices.push_back(price);names.push_back(r[0]+"_"+r[1]);
    }
    for(bool first14:{false,true}){
        auto f=output(out/(first14?"ex3_1_first14.csv":"ex3_1_cashflows_full.csv"));
        std::vector<int> selected; for(auto& d:dates){if(first14 && selected.size()==14)break;selected.push_back(d.first);}
        f<<"instrument";for(int d:selected)f<<','<<date(d);f<<'\n';
        for(size_t i=0;i<flows.size();++i){f<<names[i];for(int d:selected)f<<','<<(flows[i].count(d)?flows[i].at(d):0);f<<'\n';}
    }
    auto f=output(out/"ex3_1_prices.csv");f<<"instrument,price\n";for(size_t i=0;i<prices.size();++i)f<<names[i]<<','<<prices[i]<<'\n';
    std::cout<<"3.1: "<<flows.size()<<" x "<<dates.size()<<" cash flow matrix\n";
}
std::pair<double,double> curves(double u,const Vec& T,const Vec& x){
    double forward=x[0],integral=x[0]*u;
    for(size_t i=0;i<T.size();++i){double t=T[i],m=std::min(u,t);
        forward+=x[i+1]*(t+t*m-m*m/2);
        integral+=x[i+1]*(t*u+t*m*m/2-m*m*m/6+(u-m)*(t*m-m*m/2));
    }
    return {forward,u==0?forward:integral/u};
}
void exercise32(const fs::path& data,const fs::path& out){
    auto rows=readCSV(data/"yields_3_2.csv");Vec T,Y;
    for(size_t i=1;i<rows.size();++i){T.push_back(number(rows[i].at(0)));Y.push_back(number(rows[i].at(1)));}
    for(size_t i=0;i<T.size();++i)if(T[i]<=0 || (i && T[i]<=T[i-1]))throw std::runtime_error("Maturities must increase");
    size_t n=T.size();if(!n)throw std::runtime_error("No yields");
    auto coef=output(out/"ex3_2_coefficients.csv");coef<<"alpha,f0";for(size_t i=1;i<=n;++i)coef<<",a"<<i;coef<<'\n';
    auto metrics=output(out/"ex3_2_metrics.csv");metrics<<"alpha,yield_rmse,roughness,objective\n";
    for(double alpha:{.01,.1,1.}){
        Mat a(n+1,Vec(n+1)),H(n,Vec(n));Vec b(n+1);
        for(size_t i=0;i<n;++i){a[0][i+1]=a[i+1][0]=T[i];b[i+1]=T[i]*Y[i];
            for(size_t j=0;j<n;++j){double m=std::min(T[i],T[j]);H[i][j]=T[i]*T[j]*(1+m)-(T[i]+T[j])*m*m/2+m*m*m/3;a[i+1][j+1]=H[i][j]+(i==j?1/alpha:0);}}
        auto x=solve(a,b);coef<<alpha;for(double v:x)coef<<','<<v;coef<<'\n';
        std::ostringstream name;name<<alpha;
        auto f=output(out/("ex3_2_curves_alpha_"+name.str()+".csv"));f<<"maturity,forward_percent,yield_percent\n";
        for(int k=0;k<=500;++k){double u=T.back()*k/500.;auto c=curves(u,T,x);f<<u<<','<<c.first<<','<<c.second<<'\n';}
        // Portable SVG charts: open directly in any web browser.
        auto svg=output(out/("ex3_2_curves_alpha_"+name.str()+".svg"));
        double lo=0,hi=0;
        for(int k=0;k<=500;++k){auto c=curves(T.back()*k/500.,T,x);lo=std::min({lo,c.first,c.second});hi=std::max({hi,c.first,c.second});}
        for(double y:Y){lo=std::min(lo,y);hi=std::max(hi,y);}
        double pad=std::max(.1,(hi-lo)*.1);lo-=pad;hi+=pad;
        svg<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"400\" viewBox=\"0 0 1000 400\"><rect width=\"1000\" height=\"400\" fill=\"white\"/>";
        for(int panel=0;panel<2;++panel){double left=60+panel*500;
            svg<<"<text x=\""<<left<<"\" y=\"25\" font-size=\"18\">"<<(panel==0?"Forward":"Yield")<<" curve, alpha="<<alpha<<"</text>";
            svg<<"<path d=\"M "<<left<<" 50 V 330 H "<<left+400<<"\" fill=\"none\" stroke=\"black\"/>";
            for(int tick=0;tick<=5;++tick){double v=lo+(hi-lo)*tick/5.;svg<<"<text x=\""<<left-50<<"\" y=\""<<330-280.*tick/5<<"\" font-size=\"11\">"<<std::setprecision(3)<<v<<"%</text>";}
            svg<<"<text x=\""<<left<<"\" y=\"355\">0</text><text x=\""<<left+330<<"\" y=\"355\">"<<T.back()<<" years</text><polyline fill=\"none\" stroke=\"#1768ac\" stroke-width=\"2\" points=\""<<std::setprecision(16);
            for(int k=0;k<=500;++k){auto c=curves(T.back()*k/500.,T,x);double v=panel==0?c.first:c.second;svg<<left+400.*k/500<<','<<330-280*(v-lo)/(hi-lo)<<' ';}
            svg<<"\"/>";
            if(panel==1)for(size_t i=0;i<n;++i)svg<<"<circle cx=\""<<left+400*T[i]/T.back()<<"\" cy=\""<<330-280*(Y[i]-lo)/(hi-lo)<<"\" r=\"4\" fill=\"black\"/>";
        }
        svg<<"</svg>";
        auto fit=output(out/("ex3_2_fit_alpha_"+name.str()+".csv"));fit<<"maturity,observed_yield,fitted_yield\n";
        double mse=0,sse=0,rough=0;
        for(size_t i=0;i<n;++i){double y=curves(T[i],T,x).second,e=Y[i]-y;mse+=e*e;sse+=T[i]*T[i]*e*e;fit<<T[i]<<','<<Y[i]<<','<<y<<'\n';
            for(size_t j=0;j<n;++j)rough+=x[i+1]*(H[i][j]-T[i]*T[j])*x[j+1];}
        metrics<<alpha<<','<<std::sqrt(mse/n)<<','<<rough<<','<<rough+alpha*sse<<'\n';
    }
    std::cout<<"3.2: coefficients and curves saved\n";
}
// Jacobi diagonalization for a small real symmetric matrix, eigenvectors in columns.
void eigen(Mat a,Vec& values,Mat& vectors){
    size_t n=a.size();vectors=Mat(n,Vec(n));for(size_t i=0;i<n;++i)vectors[i][i]=1;
    bool done=false;
    for(int iter=0;iter<10000;++iter){
        size_t p=0,q=1;double largest=0;
        for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j)if(std::abs(a[i][j])>largest){largest=std::abs(a[i][j]);p=i;q=j;}
        double scale=0;for(size_t i=0;i<n;++i)scale=std::max(scale,std::abs(a[i][i]));
        if(largest<=1e-14*std::max(scale,1e-300)){done=true;break;}
        double theta=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),c=std::cos(theta),s=std::sin(theta);
        double pp=a[p][p],qq=a[q][q],pq=a[p][q];
        for(size_t k=0;k<n;++k)if(k!=p && k!=q){double kp=a[k][p],kq=a[k][q];a[k][p]=a[p][k]=c*kp-s*kq;a[k][q]=a[q][k]=s*kp+c*kq;}
        a[p][p]=c*c*pp-2*s*c*pq+s*s*qq;a[q][q]=s*s*pp+2*s*c*pq+c*c*qq;a[p][q]=a[q][p]=0;
        for(size_t k=0;k<n;++k){double vp=vectors[k][p],vq=vectors[k][q];vectors[k][p]=c*vp-s*vq;vectors[k][q]=s*vp+c*vq;}
    }
    if(!done)throw std::runtime_error("Eigensolver did not converge");
    std::vector<size_t> order(n);std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](size_t i,size_t j){return a[i][i]>a[j][j];});
    Mat sorted(n,Vec(n));values.resize(n);
    for(size_t j=0;j<n;++j){values[j]=a[order[j]][order[j]];size_t pivot=0;for(size_t i=1;i<n;++i)if(std::abs(vectors[i][order[j]])>std::abs(vectors[pivot][order[j]]))pivot=i;
        double sign=vectors[pivot][order[j]]<0?-1:1;for(size_t i=0;i<n;++i)sorted[i][j]=sign*vectors[i][order[j]];}
    vectors=sorted;
}
int monthID(const std::string& s){
    if(s.size()!=7 || s[4]!='-')throw std::runtime_error("Expected YYYY-MM");
    int y=std::stoi(s.substr(0,4)),m=std::stoi(s.substr(5,2));if(m<1 || m>12)throw std::runtime_error("Invalid month");return y*12+m-1;
}
std::string monthName(int id){std::ostringstream s;s<<id/12<<'-'<<std::setfill('0')<<std::setw(2)<<id%12+1;return s.str();}
void exercise35(const fs::path& input,const fs::path& out,const std::string& end){
    auto rows=readCSV(input);std::vector<std::string> tenors={"2","3","4","5","7","10","20","30"};
    std::vector<size_t> cols;for(auto& t:tenors){auto it=std::find(rows[0].begin(),rows[0].end(),t);if(it==rows[0].end())throw std::runtime_error("Missing tenor column "+t);cols.push_back(it-rows[0].begin());}
    auto dit=std::find(rows[0].begin(),rows[0].end(),"date");if(dit==rows[0].end())throw std::runtime_error("Missing date column");size_t dc=dit-rows[0].begin();
    std::map<int,std::vector<std::string>> levels;
    for(size_t r=1;r<rows.size();++r){int m=monthID(rows[r].at(dc));if(levels.count(m))throw std::runtime_error("Duplicate month");levels[m]=rows[r];}
    if(levels.empty())throw std::runtime_error("No observations");
    int last=end.empty()?levels.rbegin()->first:monthID(end);
    Mat x;Vec previous;
    for(int m=last-120;m<=last;++m){if(!levels.count(m))throw std::runtime_error("Missing month: "+monthName(m));Vec v;
        for(auto c:cols)v.push_back(number(levels.at(m).at(c)));
        if(!previous.empty()){Vec change(8);for(int j=0;j<8;++j)change[j]=v[j]-previous[j];x.push_back(change);}previous=v;}
    Vec mean(8);for(auto& row:x)for(int j=0;j<8;++j)mean[j]+=row[j]/120.;
    Mat Q(8,Vec(8));for(auto& row:x)for(int i=0;i<8;++i)for(int j=0;j<8;++j)Q[i][j]+=(row[i]-mean[i])*(row[j]-mean[j])/120.;
    Vec vals;Mat A;eigen(Q,vals,A);double total=std::accumulate(vals.begin(),vals.end(),0.);
    if(total<=0)throw std::runtime_error("Zero total variance");
    std::vector<std::string> pcs;for(int j=1;j<=8;++j)pcs.push_back("PC"+std::to_string(j));
    writeMat(out/"ex3_5_covariance.csv",Q,tenors);writeMat(out/"ex3_5_eigenvectors.csv",A,pcs);
    auto f=output(out/"ex3_5_variance.csv");f<<"component,eigenvalue,explained_percent,cumulative_percent\n";double cumulative=0;
    for(int j=0;j<8;++j){double pct=100*vals[j]/total;cumulative+=pct;f<<pcs[j]<<','<<vals[j]<<','<<pct<<','<<cumulative<<'\n';}
    auto changes=output(out/"ex3_5_monthly_changes.csv"),scores=output(out/"ex3_5_scores.csv"),means=output(out/"ex3_5_mean_changes.csv");
    changes<<"month";scores<<"month";for(int j=0;j<8;++j){changes<<','<<tenors[j];scores<<','<<pcs[j];}changes<<'\n';scores<<'\n';means<<"maturity,mean_change\n";
    for(int j=0;j<8;++j)means<<tenors[j]<<','<<mean[j]<<'\n';
    for(size_t k=0;k<x.size();++k){changes<<monthName(last-119+int(k));scores<<monthName(last-119+int(k));for(int j=0;j<8;++j){changes<<','<<x[k][j];double score=0;for(int i=0;i<8;++i)score+=(x[k][i]-mean[i])*A[i][j];scores<<','<<score;}changes<<'\n';scores<<'\n';}
    std::cout<<"3.5: "<<monthName(last-119)<<" to "<<monthName(last)<<", first 3 PCs = "<<100*(vals[0]+vals[1]+vals[2])/total<<"%\n";
}
int main(int argc,char** argv){
    try {
        fs::path data=argc>1?argv[1]:"data",out=argc>2?argv[2]:"results";
        std::string end=argc>3?argv[3]:"";fs::create_directories(out);
        exercise31(data,out);exercise32(data,out);exercise35(data/"snb_yields.csv",out,end);
        std::cout<<"Results: "<<fs::absolute(out)<<'\n';return 0;
    } catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
