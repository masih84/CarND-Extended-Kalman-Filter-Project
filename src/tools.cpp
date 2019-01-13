#include "tools.h"
#include "FusionEKF.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size()
		|| estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	// accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		// coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// calculate the mean
	rmse = rmse / estimations.size();

	// calculate the squared root
	rmse = rmse.array().sqrt();

	// return the result
	return rmse;
}

VectorXd Tools::cartesian_to_polar(const VectorXd& x_state) {
	VectorXd z_hat(3, 1);

	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = atan2(py, px);
	if (c3 < 0) {
		float PI = atan(1) * 4;
		c3 += 2 * PI;
	}
	float c4 = (px * vx) + (py * vy);

	// check division by zero
	if (fabs(c2) < 0.001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		c2 = 0.001;
	}

	// compute the z_hat
	z_hat << c2,
		c3,
		c4 / c2;
	return z_hat;
}

VectorXd Tools::polar_to_cartesian(const VectorXd& z_state) {
	VectorXd x_state(4, 1);

	// recover Radar measured state parameters
	float ro = z_state(0);
	float phi = z_state(1);
	float dro = z_state(2);

	// compute the x_state
	x_state << ro * cos(phi),
		 ro * sin(phi),
		dro * cos(phi),
		 dro * sin(phi);
	return x_state;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3, 4);
	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	// check division by zero
	if (fabs(c1) < 0.001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		c1 = 0.001;
	}

	// compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	return Hj;
}

float Tools::ConstrainAngle(float x) {
	float PI = atan(1) * 4;
	/*x = fmod(x + PI, 2*PI);
	if (x < 0)
		x += 2 * PI;
	return x - PI;
		*/
	x = fmod(x, 2 * PI);
	if (x < 0)
		x += 2 * PI;
	return x;

}

