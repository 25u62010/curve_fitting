#include <iostream>
#include <curve_fitting.h>
#include <opencv2/core/core.hpp>

using namespace std; 
using namespace curve_fitting;
int main( int argc, char** argv )
{
	sin_fitting sinFitting;
    double a=2.0, b=3.0, c=0.5,d=0.2;         // 真实参数值
    int N=500;                          // 数据点
    double w_sigma=0.05;                 // 噪声Sigma值
    cv::RNG rng;                        // OpenCV随机数产生器


    vector<double> x_data, y_data;      // 数据
    
    //cout<<"generating data: "<<endl;
    for ( int i=0; i<N; i++ ){
        double x = i/100.0;
        double y=a*sin ( b*x  + c ) +d + rng.gaussian ( w_sigma );
        sinFitting.addData(x,y);
    }
    
    sinFitting.Optimzer();
    double aE,bE,cE,dE;
    sinFitting.readParameters(aE,bE,cE,dE);
    printf("fx=%lfsin(%lfx+%lf)+%lf\r\n",aE,bE,cE,dE);
    double magnitude=sinFitting.readMagnitude();
    double frequence=sinFitting.readFrequence();
    cout<<"fre="<<frequence<<" magnitude="<<magnitude<<endl;
    return 0;
}
