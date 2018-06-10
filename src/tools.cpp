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
  rmse << 0, 0, 0, 0;
  int vector_size = estimations.size();

  if (vector_size == 0 ||
          (estimations.size() != ground_truth.size())) {
      return rmse;
  }

  for (int i = 0; i < vector_size; i++) {
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array() * residual.array();
      rmse += residual;
  }

  /*
   * rmse.mean();
   */
  rmse = rmse / vector_size;
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3, 4);

    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj.setZero();
    }

    cout << "px " << px << " py " << py << " vx " << vx << " vy " << vy << endl;
    cout << "c1 " << fabs(c1) << " c2 " << c2 << " c3 " << c3 << endl;

    if (isinf(c1) || isinf(c2) || isinf(c3) ||
            isinf(px / c2) || isinf(py / c2) ||
            isinf(-(py/c1)) || isinf(px / c1) ||
            isinf(py*(vx*py - vy*px)/c3) ||
            isinf(px*(px*vy - py*vx)/c3) ||
            isinf(px/c2) || isinf(py/c2)) {
        cout << __func__ << " - Error - Division by Zero" << endl;
        return Hj.setZero();
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;


//
//    float h00 = isinf(px / c2) ? 0 : (px / c2);
//    h00 = isinf(c2) ? 0 : h00;
//    float h01 = isinf(py / c2) ? 0 : (py / c2);
//    h01 = isinf(c2) ? 0 : h01;
//    float h10 = isinf(-(py/c1)) ? 0 : (-(py/c1));
//    h10 = isinf(c1) ? 0 : h10;
//    float h11 = isinf(px / c1) ? 0 : (px / c1);
//    h11 = isinf(c1) ? 0 : h11;
//    float h20 = isinf(py*(vx*py - vy*px)/c3) ? 0 : (py*(vx*py - vy*px)/c3);
//    h20 = isinf(c3) ? 0 : h20;
//    float h21 = isinf(px*(px*vy - py*vx)/c3) ? 0 : (px*(px*vy - py*vx)/c3);
//    h21 = isinf(c3) ? 0 : h21;
//    float h22 = isinf(px/c2) ? 0 : (px/c2);
//    h22 = isinf(c2) ? 0 : h22;
//    float h23 = isinf(py/c2) ? 0 : (py/c2);
//    h23 = isinf(c2) ? 0 : h23;
//
//    //compute the Jacobian matrix
//    Hj << h00, h01, 0, 0,
//            h10, h11, 0, 0,
//            h20, h21, h22, h23;

    return Hj;
}
