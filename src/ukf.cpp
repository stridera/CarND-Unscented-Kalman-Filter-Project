#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

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
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.2;

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

    // Initially set to fase, true after initialization with first measurement
    is_initalized_ = false;

    //set state dimension
    n_x_ = 5;

    //define spreading parameter
    lambda_ = 3 - n_x_;

    //set augmented dimension
    n_aug_ = 7;

    //vector for weights
    weights_ = VectorXd(2*n_aug_+1);

    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
}

UKF::~UKF() {}

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
    if ((!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::Radar) ||
        (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::Laser))
        return;

    if (!is_initalized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::Radar) {
            float ro = meas_package.raw_measurements_[0];
            float phi = tools::NormalizePhi(meas_package.raw_measurements_[1]);
            float ro_dot = meas_package.raw_measurements_[2];
            
            float px = ro * cos(phi);
            float py = ro * sin(phi);

            x_ = px, py, ro_dot, 0., 0.;
        } else {
            float px = meas_package.raw_measurements_[0];
            float py = meas_package.raw_measurements_[1];
            x_ << px, py, 0., 0., 0.;
        }

        previous_timestamp_ = meas_package.timestamp_;
        is_initalized_ = true;
        return;
    }

    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;

    prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::Radar) {
        UpdateRadar();
    } else {
        UpdateLidar;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // Generate sigma points:
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0.0);
    AugmentedSigmaPoints(&Xsig_aug);
    Xsig_pred_ = SigmaPointPrediction(Xsig_aug, delta_t);
    PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 **/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
    // Lidar is 2dim (x, y)
    int n_z = 2;

    float px = meas_package.raw_measurements_[0];
    float py = meas_package.raw_measurements_[1];
    VectorXd z = VectorXd(n_z);
    z << px, py;

    // Kalman Filter (From Project 1)
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

    // Calculate NIS: z_diff^T * S^-1 * z_diff
    nis_lidar_ = y.transpose() * Si * y
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 **/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // Radar is 3dim (ro, phi, ro_dot)
    int n_z = 3;

    // Break out and normalize measurements
    float ro = meas_package.raw_measurements_[0];
    float phi = tools::NormalizePhi(meas_package.raw_measurements_[1]);
    float ro_dot = meas_package.raw_measurements_[2];

    VectorXd(n_z) z;
    z << ro, phi, ro_dot;

    // Todo:  Figure out z_pred and S
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);
    MatrixXd Zsig = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    PredictRadarMeasurement(&z_pred, &S, &Zsig);

    //calculate cross correlation matrix
    MatrixXd(n_z, n_z) Tc;
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        z_diff(1) = tools::NormalizePhi(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        x_diff(3) = tools::NormalizePhi(x_diff(3));

        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    z_diff(1) = tools::NormalizePhi(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

    // Calculate NIS: z_diff^T * S^-1 * z_diff
    nis_radar_ = z_diff.transpose()*S.inverse()*z_diff ;
}

/**
 * From Lesson 7, Secion 17-18
 **/
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {


    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a*std_a;
    P_aug(6,6) = std_yawdd*std_yawdd;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
    }
 
    //write result
    *Xsig_out = Xsig_aug;
}
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

    //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //create example sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
         Xsig_aug <<
        5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
            1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
        2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
        0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
        0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
                 0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
                 0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

    double delta_t = 0.1; //time diff in sec

    //predict sigma points
    for (int i = 0; i< 2*n_aug+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
                px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
                py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
                px_p = p_x + v*delta_t*cos(yaw);
                py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }

    //write result
    *Xsig_out = Xsig_pred;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

    //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //define spreading parameter
    double lambda = 3 - n_aug;

    //create example matrix with predicted sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred <<
                 5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
                     1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
                    2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
                 0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
                    0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    //create vector for weights
    VectorXd weights = VectorXd(2*n_aug+1);
    
    //create vector for predicted state
    VectorXd x = VectorXd(n_x);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x, n_x);

    // set weights
    double weight_0 = lambda/(lambda+n_aug);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug+lambda);
        weights(i) = weight;
    }

    //predicted state mean
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        x = x+ weights(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P = P + weights(i) * x_diff * x_diff.transpose() ;
    }

    //write result
    *x_out = x;
    *P_out = P;
}

/**
 * Lesson 7: Section 26-27
 **/
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug+1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        z_diff(1) = tools::NormalizePhi(z_diff(1));
        S = S + weights(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std_radr*std_radr, 0, 0,
          0, std_radphi*std_radphi, 0,
          0, 0,std_radrd*std_radrd;
    S = S + R;

    //write result
    *z_out = z_pred;
    *S_out = S;
    *Zsig_out = Zsig;
}