//
// Created by Gang Wang on 6/11/2018.
//

#include <iostream>
#include <vector>
#include <fstream>

#include "Eigen/Dense"
#include "measurement_package.h"
#include "FusionEKF.h"
#include "tools.h"


using namespace std;
using Eigen::VectorXd;

static void ekfProcess() {
    FusionEKF fusionEKF;
    Tools tools;
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;
    vector<MeasurementPackage> measurement_pack_list;

    string in_file_name_ = "data/obj_pose-laser-radar-synthetic-input.txt";
    ifstream in_file(in_file_name_.c_str(), std::ifstream::in);

    if (!in_file.is_open()) {
        cout << "Cannot open input file: " << in_file_name_ << endl;
    }

    string line;

    while (getline(in_file, line)) {
        MeasurementPackage meas_package;

        istringstream iss(line);
        string sensor_type;
        iss >> sensor_type;	//reads first element from the current line
        int64_t timestamp;

        if (sensor_type.compare("L") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float px;
            float py;
            iss >> px;
            iss >> py;
            meas_package.raw_measurements_ << px, py;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
        } else if (sensor_type.compare("R") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float ro;
            float theta;
            float ro_dot;
            iss >> ro;
            iss >> theta;
            iss >> ro_dot;
            meas_package.raw_measurements_ << ro, theta, ro_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
        }

        float x_gt;
        float y_gt;
        float vx_gt;
        float vy_gt;
        iss >> x_gt;
        iss >> y_gt;
        iss >> vx_gt;
        iss >> vy_gt;
        VectorXd gt_values(4);
        gt_values(0) = x_gt;
        gt_values(1) = y_gt;
        gt_values(2) = vx_gt;
        gt_values(3) = vy_gt;
        ground_truth.push_back(gt_values);

        fusionEKF.ProcessMeasurement(meas_package);

        VectorXd estimate(4);
        double p_x = fusionEKF.ekf_.x_(0);
        double p_y = fusionEKF.ekf_.x_(1);
        double v1 = fusionEKF.ekf_.x_(2);
        double v2 = fusionEKF.ekf_.x_(3);

        estimate(0) = p_x;
        estimate(1) = p_y;
        estimate(2) = v1;
        estimate(3) = v2;

        estimations.push_back(estimate);

        VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);

        cout << "estimation " << estimate << endl;
        cout << "RMSE " << RMSE << endl;
    }
}

int main() {
    ekfProcess();
    return 0;
}