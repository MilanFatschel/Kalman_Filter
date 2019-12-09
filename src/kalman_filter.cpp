#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

  // Prediction Step
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

  // Update measurements
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse(); 

  // New State
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Predict State
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Radar Measurements 
  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rho_dot;
  if(rho == 0) rho_dot = 0;
  else rho_dot = (px * vx + py * vy) / rho; 

  // Prediction to compare to measured
  VectorXd zPrediction(3);
  zPrediction << rho, phi, rho_dot;

  // Measurement update
  VectorXd y = z - zPrediction;  
  double offsetValue = y(1) + M_PI;  
  y(1) = (offsetValue - (floor(offsetValue / (2 * M_PI)) * (2 * M_PI))) - M_PI;
 
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse(); 

  // New State 
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;


}
