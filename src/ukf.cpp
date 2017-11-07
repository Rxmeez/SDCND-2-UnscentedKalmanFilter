#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //[px, py, v, yaw_angle, yaw_rate]
  n_x_ = x_.size();

  //augment state dimensions + [nu, psi]
  n_aug_ = n_x_ + 2;

  //initial corvariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  //sigma points
  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //sigma point spreading variable
  lambda_ = 3 - n_aug_;

  //weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  //measurement noise covariance matrices initialization
  R_rader_ = MatrixXd(3, 3);

  R_rader_ << std_radr_ * std_radr_,  0,  0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2, 2);

  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy;
}

UKF::~UKF() {}

void UKF::NormAng(double *ang) {
  //Angle normalization to [-Pi, Pi]
  while (*ang > M_PI) *ang -= 2. * M_PI;
  while (*ang < -M_PI) *ang += 2. * M_PI;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){
    x_ << 0, 0, 0, 0, 0;
    //initialize covariance matrix
    P << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0,
         0, 0, 1, 0, 0,
         0, 0, 0, 1, 0,
         0, 0, 0, 0, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
       //polar to cartesian coordinates
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rhodot = measurement_pack.raw_measurements_(2);

      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rhodot * cos(phi);
      float vy = rhodot * sin(phi);
      float v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER){
      //velocity == 0
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);

      x_ << px, py, 0, 0;

      //dealing with < 0.001 problems
      if (fabs(px) < 0.001 and fabs(py) < 0.001){
        px = 0.001;
        py = 0.001;
      }
    }

    //initialize weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < weights_.size(); i++){
      weights_(i) = 0.5 / (n_aug_ + lambda_);
    }

    //initial timestamp for dt
    time_us_ = measurement_pack.timestamp_;

    //complete intializing
    is_initialized_ = true;

    return;
  }

  // dt in seconds
  float dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  time_us_ = measurement_pack.timestamp_

  Prediction(dt);
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(measurement_pack);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(measurement_pack);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  float delta_t2 = delta_t * delta_t;

  //initial augmented state VectorXd
  VectorXd x_aug_ = VectorXd(n_aug_);

  //augment state covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_aug_);

  //fill in augmented matrices
  x_aug_.fill(0);
  x_aug_.head(n_x_) = x_;

  P.aug_.fill(0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  //square root of P matrix
  MatrixXd L = P_aug_.llt().matrixL();

  Xsig_aug_.col(0) = x_aug_;
  for (int = 0; i < n_aug_; i++){
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
