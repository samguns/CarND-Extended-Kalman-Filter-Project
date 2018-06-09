#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
    PI = atan(1.0)*4;
}

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
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     * update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
     * update the state by using Extended Kalman Filter equations
     */

    double px = x_[0];
    double py = x_[1];
    double vx = x_[2];
    double vy = x_[3];
    double c1 = px * px + py * py;
    double rho;
    double theta;
    double ro_dot;
    if ((fabs(c1) < 0.0001) ||
        ((px == 0) && (py == 0))) {
        rho = 0;
        theta = 0;
        ro_dot = 0;
    } else {
        rho = sqrt(px * px + py * py);
        theta = atan2(py, px);
        ro_dot = px * vx + py * vy / rho;
    }
    VectorXd z_pred = VectorXd(3);
    z_pred << rho, theta, ro_dot;

    VectorXd y = z - z_pred;
    while (y[1] < -PI) {
        y[1] += 2 * PI;
    }
    while (y[1] > PI) {
        y[1] -= 2 * PI;
    }
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
