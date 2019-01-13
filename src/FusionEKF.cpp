#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;



/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;


  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */

  // set the acceleration noise components
  noise_ax = 5;
  noise_ay = 5;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * Convert radar from polar to cartesian coordinates.
     */
	 // initializing matrices
	  ekf_.P_ = MatrixXd(4, 4);
	  ekf_.F_ = MatrixXd(4, 4);
	  ekf_.Q_ = MatrixXd(4, 4);

	  ekf_.x_ = VectorXd(4);

	 // state covariance matrix P
	  ekf_.P_ << 1, 0, 0, 0,
		  0, 1, 0, 0,
		  0, 0, 1000, 0,
		  0, 0, 0, 1000;

	  // the initial transition matrix F_
	  ekf_.F_ << 1, 0, 1, 0,
		  0, 1, 0, 1,
		  0, 0, 1, 0,
		  0, 0, 0, 1;

	  //systematic covariance matrix 
	  ekf_.Q_ << 1.0 / 4.0 * noise_ax, 0, 1.0 / 2.0 * noise_ax, 0,
		  0, 1.0 / 4.0 * noise_ay, 0, 1. / 2.0 * noise_ay,
		  1. / 2.0 * noise_ax, 0, 1. * noise_ax, 0,
		  0, 1. / 2.0 * noise_ay, 0, 1. * noise_ay;

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ << 1, 1, 1, 1;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

		// Convert radar from polar to cartesian coordinates 
		//         and initialize state.
			
		ekf_.x_ = tools.polar_to_cartesian(measurement_pack.raw_measurements_);
		// first measurement

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

		// Initialize state.
		// set the state with the initial location and zero velocity
		ekf_.x_ << measurement_pack.raw_measurements_[0],
			measurement_pack.raw_measurements_[1],
			0,
			0;

	}

	ekf_.Predict();

    // done initializing, no need to predict or update
    is_initialized_ = true;
	cout << "done initializing" << endl;

    return;
  }



  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   // 1. Modify the F matrix so that the time is integrated
   // Modified F to include time 

  noise_ax = 5;
  noise_ay = 5;


  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  if (dt > 10) {
	  dt = 0;
  }
  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;

  // 2. Set the process covariance matrix Q

  ekf_.Q_ << dt4 / 4.0 * noise_ax, 0, dt3 / 2.0 * noise_ax, 0,
	  0, dt4 / 4.0 * noise_ay, 0, dt3 / 2.0 * noise_ay,
	  dt3 / 2.0 * noise_ax, 0, dt2 * noise_ax, 0,
	  0, dt3 / 2.0 * noise_ay, 0, dt2 * noise_ay;

  // Predict 
  ekf_.Predict();


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //  Radar updates

	  float phi = measurement_pack.raw_measurements_(1);

	// initializing matrices

	  ekf_.R_ = MatrixXd(3, 3);

	  // LASER measurement covariance
	  ekf_.R_ = R_radar_;

	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
		// initializing matrices
	  ekf_.R_ = MatrixXd(2, 2);
	  ekf_.H_ = MatrixXd(2, 4);

	  // LASER measurement matrix
	  ekf_.H_ = H_laser_;

	  // LASER measurement covariance
	  ekf_.R_ = R_laser_;

	  ekf_.Update(measurement_pack.raw_measurements_);

  }
  // print the output
  cout << "final x_ = " << ekf_.x_ << endl;
  cout << "final P_ = " << ekf_.P_ << endl;
}


