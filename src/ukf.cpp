#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5,5);
  P_ << .15, 0, 0, 0, 0,
      0, .15, 0, 0 ,0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_aug_;
  previous_timestamp_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;    x_ = VectorXd(5);
    x_ << 1, 1, 1, 1, 1;
    Tools tools;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      x_ = tools.InitializeRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  } 

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ =  meas_package.timestamp_;

  //PREDICTION STEP
  Prediction(dt);

  //UPDATE STEP
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
  else {
    UpdateLidar(meas_package);
  }
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t);
  PredictMeanCovariance(Xsig_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  n_z = 2;
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred_ = VectorXd(n_z);
  S_ = MatrixXd(n_z,n_z);
  VectorXd z_ = VectorXd(n_z);
  z_ << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1];
  PredictMeasurementSigmaPoints(n_z);
  UpdateStates(z_,n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  n_z = 3;
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred_ = VectorXd(n_z);
  S_ = MatrixXd(n_z,n_z);
  //measurement vector
  VectorXd z_ = VectorXd(n_z);
  z_ << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];
  PredictMeasurementSigmaPoints(n_z);
  UpdateStates(z_,n_z);
  
}
MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t){
  MatrixXd Xsig_predict_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for(int i = 0; i < 2*n_aug_+1; i++){
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw_angle = Xsig_aug(3,i);
    double yaw_rate = Xsig_aug(4,i);
    double longacc_noise = Xsig_aug(5,i);
    double yaw_acc_noise = Xsig_aug(6,i);
    
    double px_pred;
    double py_pred;
    double v_pred;
    double yaw_pred;
    double yaw_rate_pred;
    if(yaw_rate >.001){
      px_pred = px + (v*(sin(yaw_angle+(yaw_rate*delta_t)) - sin(yaw_angle)))/yaw_rate;
      py_pred = py + (v*(-cos(yaw_angle+(yaw_rate*delta_t)) +cos(yaw_angle)))/yaw_rate;
      yaw_pred = yaw_angle + (yaw_rate*delta_t);
    }
    else{
      px_pred = px + (v*cos(yaw_angle)*delta_t);
      py_pred = py + (v*sin(yaw_angle)*delta_t);
      yaw_pred = yaw_angle + (yaw_rate*delta_t);   
    }
    //add noise componenet
    px_pred = px_pred + (delta_t*delta_t*cos(yaw_angle)*longacc_noise)/2;
    py_pred = py_pred + (delta_t*delta_t*sin(yaw_angle)*longacc_noise)/2;
    v_pred = v + (delta_t*longacc_noise);
    yaw_pred = yaw_pred + (delta_t*delta_t*yaw_acc_noise)/2;
    yaw_rate_pred = yaw_rate + (delta_t*yaw_acc_noise);
    
    //set sigma points to the matrix
    Xsig_predict_(0,i)=px_pred;
    Xsig_predict_(1,i)=py_pred;
    Xsig_predict_(2,i)=v_pred;
    Xsig_predict_(3,i)=yaw_pred;
    Xsig_predict_(4,i)=yaw_rate_pred;
  }
  return Xsig_predict_;
}

void UKF::PredictMeanCovariance(MatrixXd Xsig_pred){
  weights_ = VectorXd(2*n_aug_+1); 
  for(int i = 0; i<2*n_aug_+1; i++){
      double w = 0;
      if(i == 0){
          w = lambda_/(lambda_+n_aug_);
      }
      else{
          w = 1/(2*(lambda_+n_aug_));
      }
      weights_(i) = w;
  }
  x_.fill(0);
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++){
      x_ = x_ + weights_(i) * Xsig_pred.col(i);
  }
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    //angle normalization
    if (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    if (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::PredictMeasurementSigmaPoints(int n_z){
  for(int i = 0; i<2*n_aug_+1; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw_angle = Xsig_pred_(3,i);
    double yaw_rate = Xsig_pred_(4,i);
    
    if(n_z==3){
      double rho = sqrt(px*px+py*py);
      double phi = atan2(py,px);
      double rho_dot = ((px*cos(yaw_angle)+py*sin(yaw_angle))*v)/rho;
      
      Zsig_(0,i)=rho;
      Zsig_(1,i)=phi;
      Zsig_(2,i)=rho_dot;
    }
    else if(n_z == 2){
      Zsig_(0,i)=px;
      Zsig_(1,i)=py;
    }
  }

  z_pred_.fill(0);
  for(int i=0; i<2*n_aug_+1;i++){
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  S_.fill(0.0);

  R_ = MatrixXd(n_z,n_z);
  if(n_z == 3){
    R_ << std_radr_*std_radr_,0,0,
      0,std_radphi_*std_radphi_,0,
      0,0,std_radrd_*std_radrd_;
  }
  else if(n_z == 2){
    R_ << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;
  }

  for(int i =0; i<2*n_aug_+1; i++){
      VectorXd z_diff = Zsig_.col(i) - z_pred_;

      if(n_z == 3){
        if(z_diff(1)<-M_PI){z_diff(1)=z_diff(1)+2*M_PI;}
        if(z_diff(1)>M_PI){z_diff(1)=z_diff(1)-2*M_PI;}
      }

      S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_ = S_ + R_;
}

void UKF::UpdateStates(VectorXd z, int n_z){
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){
    VectorXd z_diff = Zsig_.col(i)-z_pred_;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    if(n_z == 3){
      //phi normalization
      if(z_diff(1)>M_PI){z_diff(1) -= 2*M_PI;}
      if(z_diff(1)<-M_PI){z_diff(1) += 2*M_PI;}
      //yaw angle normalization
      if(x_diff(3)>M_PI){x_diff(3) -= 2*M_PI;}
      if(x_diff(3)<-M_PI){x_diff(3) += 2*M_PI;}
    }
    
    //calculating cross covariance matrix
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = Tc * S_.inverse();
  VectorXd z_diff = z - z_pred_;

  if(n_z == 3){
    //angle normalization
    if (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    if (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  }
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  //NIS CALCULATION
  double nis = z_diff.transpose() * S_.inverse() * z_diff;
  if(n_z == 3)
    WriteNIS(nis, 'R');
  else if(n_z == 2){ WriteNIS(nis, 'L');}
}

void UKF::WriteNIS(double nis, char sensor)
{
    ofstream outputfile;
    if(sensor == 'L'){
      outputfile.open("D:/SDCND Projects/Project 7/CarND-Unscented-Kalman-Filter-Project-master/src/dataLidar.txt");
    }
    else{
      outputfile.open("D:/SDCND Projects/Project 7/CarND-Unscented-Kalman-Filter-Project-master/src/dataRadar.txt");
    }
    outputfile << nis << endl;
}

