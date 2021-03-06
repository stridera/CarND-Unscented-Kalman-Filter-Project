#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
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

double Tools::NormalizePhi(double phi)
{
	while (phi >  M_PI) phi -= 2.0 * M_PI;
    while (phi < -M_PI) phi += 2.0 * M_PI;
    return phi;
}


// Pretty Print Functions.  Used to easily watch items while debugging.
std::string Tools::pp(VectorXd v1) {
    std::ostringstream oss;
    oss << v1;
    return oss.str();
}

std::string Tools::pp(MatrixXd v1) {
    std::ostringstream oss;
    oss << v1;
    return oss.str();
}