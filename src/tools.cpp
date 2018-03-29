#include <iostream>
#include "tools.h"
#include "measurement_package.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size() || estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}
	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
  	rmse += residual;
	}
	//calculate the mean
	rmse = rmse/estimations.size();
	//calculate the squared root
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}

VectorXd Tools::InitializeRadar(MeasurementPackage measurement_pack){
  VectorXd x_state(5); 
  double rho = measurement_pack.raw_measurements_[0];
  double phi = measurement_pack.raw_measurements_[1];
  double px = rho*cos(phi);
  double py = rho*sin(phi);
  double v = 0;
  double yaw_angle = 0.3;
  double yaw_rate  = 0.2;
  if(px<0.00001){px = 0.00001;}
  if(py<0.00001){py = 0.00001;}  
  x_state<<px,py,v,yaw_angle,yaw_rate;
  return x_state;
}