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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  is_initialized_ = false;

  //initialize timestamp
	time_us_ = 0;

  //[px, py, v, yaw_angle, yaw_rate]
  n_x_ = x_.size();

  //augment state dimensions + [nu, psi]
  n_aug_ = n_x_ + 2;

  //number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  //predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //sigma point spreading variable
  lambda_ = 3 - n_aug_;

  //weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  //measurement noise covariance matrices initialization
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_,  0,  0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;
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
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){
    //initialize covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //reading in the data
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rhodot = measurement_pack.raw_measurements_[2];

      //polar to cartesian coordinates
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rhodot * cos(phi);
      float vy = rhodot * sin(phi);
      float v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER){
      //velocity == 0
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;

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

  // timesteps between measurements, dt in seconds
  double dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  time_us_ = measurement_pack.timestamp_;

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
  //initial augmented state VectorXd
  VectorXd x_aug_ = VectorXd(n_aug_);

  //augment state covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_sig_);

  //sigma points
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //fill in augmented matrices
  x_aug_.fill(0);
  x_aug_.head(n_x_) = x_;

  P_aug_.fill(0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  //square root of P matrix
  MatrixXd L = P_aug_.llt().matrixL();

  Xsig_aug_.col(0) = x_aug_;
  for (int i= 0; i < n_aug_; i++){
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i < n_sig_; i++){
    //extract values for better readability
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    double delta_t2 = delta_t * delta_t;
    double yawdd_dt = yawd * delta_t;

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd) * (sin(yaw + yawdd_dt) - sin(yaw));
        py_p = p_y + (v/yawd) * (cos(yaw) - cos(yaw+yawdd_dt));
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t2*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t2*sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t2;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  // Predicted state mean
  x_ = Xsig_pred_ * weights_; // vectorised sum
  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization
    NormAng(&(x_diff(3)));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
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
  // Set measurement dimension
  int n_z = 2;
  // Create matrix for sigma points in measurement space
  // Transform sigma points into measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);
  UpdateUKF(meas_package, Zsig, n_z);
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

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // Measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);          //r
    Zsig(1,i) = atan2(p_y,p_x);                   //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
  }
  UpdateUKF(meas_package, Zsig, n_z);
}

// Universal update function
void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred  = Zsig * weights_;
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Angle normalization
    NormAng(&(z_diff(1)));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    R = R_lidar_;
  }
  S = S + R;

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
      // Angle normalization
      NormAng(&(z_diff(1)));
    }
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalization
    //--> CHANGE:: NormAng(&(x_diff(3)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Measurements
  VectorXd z = meas_package.raw_measurements_;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // Residual
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    // Angle normalization
    NormAng(&(z_diff(1)));
  }
  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate NIS
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
	   NIS_radar_ = z.transpose() * S.inverse() * z;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
	   NIS_laser_ = z.transpose() * S.inverse() * z;
  }
}
