#include <curve_fitting.h>

using namespace curve_fitting;

sin_fitting::sin_fitting()
{
	abcd_estimate[0]=1;abcd_estimate[1]=1;abcd_estimate[2]=0;abcd_estimate[3]=0;
	typedef g2o::BlockSolver< g2o::BlockSolverTraits<4,1> > Block;  // 每个误差项优化变量维度为4，误差值维度为1
    Block::LinearSolverType* linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>(); // 线性方程求解器
    Block* solver_ptr = new Block( linearSolver );      // 矩阵块求解器
    // 梯度下降方法，从GN, LM, DogLeg 中选
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg( solver_ptr );
    // g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton( solver_ptr );
    // g2o::OptimizationAlgorithmDogleg* solver = new g2o::OptimizationAlgorithmDogleg( solver_ptr );
    
    optimizer.setAlgorithm( solver );   // 设置求解器
    optimizer.setVerbose( true );       // 打开调试输出
    
    // 往图中增加顶点
    vertex = new CurveFittingVertex();
    vertex->setEstimate( abcd_estimate );
    vertex->setId(0);
    optimizer.addVertex( vertex );
}
void sin_fitting::addData(double x,double y)
{
	int i=dataNum;
	dataNum++;
	x_data.push_back ( x );
    y_data.push_back ( y );
	CurveFittingEdge* edge = new CurveFittingEdge( x );
	edge->setId(i);
	edge->setVertex( 0, vertex );                // 设置连接的顶点
	edge->setMeasurement( y_data[i] );      // 观测数值
	edge->setInformation( Eigen::Matrix<double,1,1>::Identity()); // 信息矩阵：协方差矩阵之逆
	optimizer.addEdge( edge );
}
void sin_fitting::Optimzer()
{
	optimizer.initializeOptimization();
    optimizer.optimize(100);
    abcd_estimate = vertex->estimate();
    while(abcd_estimate[2]<0||abcd_estimate[2]>3.14){   
    	if(abcd_estimate[2]>3.14){
    		abcd_estimate[2]-=3.14;
    		abcd_estimate[0]*=-1;
    	}
    	if(abcd_estimate[2]<0){
    		abcd_estimate[2]+=3.14;
    		abcd_estimate[0]*=-1;
    	}
    }
}
bool sin_fitting::readParameters(double& a,double& b,double& c,double& d)
{
	a=abcd_estimate[0];
	b=abcd_estimate[1];
	c=abcd_estimate[2];
	d=abcd_estimate[3];
}
double sin_fitting::readFrequence()
{
	return abcd_estimate[1]/2*3.14;
}
double sin_fitting::readMagnitude()
{
	return abcd_estimate[0];
}

