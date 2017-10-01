#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
	rmse << 0,0,0,0;

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
	
	if(estimations.size() != ground_truth.size()
		|| estimations.size() == 0){
	 	cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
		VectorXd c = (estimations[i] - ground_truth[i]);
		c = c.array() * c.array();
		rmse += c;
	}

	//calculate the mean
	// ... your code here
	rmse = rmse/estimations.size();

	//calculate the squared root
	// ... your code here
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE
	float sum_sq = pow(px, 2) + pow(py, 2);

	//check division by zero
  if (fabs(px) < 0.0001 and fabs(py) < 0.0001){
	  px = 0.0001;
	  py = 0.0001;
  }

  if(fabs(sum_sq) < 0.0000001){
		sum_sq = 0.0000001;
  }

	//compute the Jacobian matrix

	Hj << px/pow(sum_sq, 0.5), py/pow(sum_sq, 0.5), 0, 0,
	      -py/sum_sq         , px/sum_sq          , 0, 0, 
	      py*(vx*py - vy*px)/pow(sum_sq, 1.5), px*(vy*px-vx*py)/pow(sum_sq, 1.5), px/pow(sum_sq,0.5), py/pow(sum_sq, 0.5);

	return Hj;
}
