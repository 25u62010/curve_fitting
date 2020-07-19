#include <iostream>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <Eigen/Core>
#include <cmath>
#include <chrono>
using namespace std;
namespace curve_fitting {
class CurveFittingVertex: public g2o::BaseVertex<4, Eigen::Vector4d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    virtual void setToOriginImpl() // 重置
    {
        _estimate << 1.,1.,0,0;
    }
    
    virtual void oplusImpl( const double* update ) // 更新
    {
        _estimate += Eigen::Vector4d(update);
    }
    // 存盘和读盘：留空
    virtual bool read( istream& in ) {}
    virtual bool write( ostream& out ) const {}
};

// 误差模型 模板参数：观测值维度，类型，连接顶点类型
class CurveFittingEdge: public g2o::BaseUnaryEdge<1,double,CurveFittingVertex>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    CurveFittingEdge( double x ): BaseUnaryEdge(), _x(x) {}
    // 计算曲线模型误差
    void computeError()
    {
        const CurveFittingVertex* v = static_cast<const CurveFittingVertex*> (_vertices[0]);
        const Eigen::Vector4d abc = v->estimate();
        _error(0,0) = _measurement - abc(0,0)*std::sin( abc(1,0)*_x + abc(2,0)) - abc(3,0);

    }
    virtual void linearizeOplus()
    {
	  CurveFittingVertex* vi = static_cast<CurveFittingVertex*>(_vertices[0]);
	  const Eigen::Vector4d xyz = vi->estimate();

	  double a = xyz[0];
	  double b = xyz[1];
	  double c = xyz[2];
	  double d = xyz[3];

	  _jacobianOplusXi(0,0) = -std::sin(b*_x+c);
	  _jacobianOplusXi(0,1) = -a*_x*std::cos(b*_x+c);
	  _jacobianOplusXi(0,2) = -a*std::cos(b*_x+c);
	  _jacobianOplusXi(0,3) = -1.;
	}
	virtual bool read( istream& in ) {}
	virtual bool write( ostream& out ) const {}
public:
    double _x;  // x 值， y 值为 _measurement
};

class sin_fitting
{
	public:
		sin_fitting();
		void addData(double x,double y);
		void Optimzer();
		bool readParameters(double& a,double& b,double& c,double& d);
		double readFrequence();
		double readMagnitude();
	private:
		g2o::SparseOptimizer optimizer;     // 图模型
		Eigen::Vector4d abcd_estimate;
		vector<double> x_data, y_data;	  
		int dataNum=0; 
		CurveFittingVertex* vertex;  
};
}
