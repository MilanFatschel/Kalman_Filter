#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    // Calculates rmse in order to evaluate how accurate
    // the EKF is 

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // Check to make sure inputs are of same size
    // or have no data at all
    if(estimations.size() == 0 || ground_truth.size() == 0
        || estimations.size() != ground_truth.size())
    {
        std::cout << "ERROR : Invalid data in for RMSE" << std::endl;
        return rmse;
    } 


    // Start calculating each data entry in estimations 
    // and ground truths

    for(unsigned int i = 0; i < estimations.size(); i++)
    {
        // Calculate residual
        VectorXd residual = estimations[i] - ground_truth[i];

        // Square each entry
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    // Calculate the mean 
    rmse = rmse / estimations.size();

    // Calculate the squared root 
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);

    // State Parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // Precompute a set of terms
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);

    // Check division by 0
    if(fabs(c1) < 0.0001)
    {
        std::cout << "Error: Division by Zero" << std::endl;
        return Hj;
    }

    // Compute Jacobian Matrix
    Hj << (px / c2), (py / c2), 0, 0,
            -(py / c1), (px / c1), 0, 0,
            py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

    return Hj;
}
