#include "RegularizedStokeslet.h"

RegularizedStokeslet::RegularizedStokeslet(vector<elasticRod *> m_rod_vec, timeStepper &m_stepper, double m_viscosity,
                                           double m_epsilon) {
    rod_vec = m_rod_vec;
    stepper = &m_stepper;
    viscosity = m_viscosity;
    epsilon = m_epsilon;
    numRod = rod_vec.size();
    nv = rod_vec[0]->nv;

    // U = AF
    // matrix (nRods * 3nv, nRod * 3nv)

    A = MatrixXd::Zero(numRod * 3 * nv, numRod * 3 * nv);
    VelocityVec = VectorXd::Zero(numRod * 3 * nv); //U
    ViscousForce = VectorXd::Zero(numRod * 3 * nv); //F

    Id3 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;
}

RegularizedStokeslet::~RegularizedStokeslet() {
    ;
}


void RegularizedStokeslet::prepareForViscousForce() {
    A.setZero();
    VelocityVec.setZero();
    ViscousForce.setZero();

    // loop for velocity idx
    int idx = 0;
    for (int n = 0; n < rod_vec.size(); n++) {
        rod_1 = rod_vec[n]; // get rod
        for (int i = 0; i < rod_1->nv; i++) {
            uPos = rod_1->getVertexOld(i); //Pose
            uVelocity = rod_1->getVelocityOld(i); //Velocity

            VelocityVec(3 * idx) = uVelocity(0);
            VelocityVec(3 * idx + 1) = uVelocity(1);
            VelocityVec(3 * idx + 2) = uVelocity(2);
            // idx++;
            // determine the discretized matrix for each node
            int idx1 = 0;
            for (int n1 = 0; n1 < rod_vec.size(); n1++) {
                rod_2 = rod_vec[n1];
                for (int j = 0; j < rod_2->ne; j++) {
                    y_0 = rod_2->getVertexOld(j);
                    y_1 = rod_2->getVertexOld(j + 1);
                    x_0 = uPos - y_0;
                    x_1 = uPos - y_1;

                    vDirection = y_0 - y_1;
                    edgeLength = vDirection.norm();

                    R_0 = sqrt(x_0.norm() * x_0.norm() + epsilon * epsilon);
                    R_1 = sqrt(x_1.norm() * x_1.norm() + epsilon * epsilon);

                    computeTCoeffs();

                    M2 = ((T1_1 + epsilon * epsilon * T1_3) * Id3 + T1_3 * (x_0 * x_0.adjoint()) +
                          T2_3 * (x_0 * vDirection.adjoint() + vDirection * x_0.adjoint()) +
                          T3_3 * (vDirection * vDirection.adjoint()));
                    M1 = ((T0_1 + epsilon * epsilon * T0_3) * Id3 + T0_3 * (x_0 * x_0.adjoint()) +
                          T1_3 * (x_0 * vDirection.adjoint() + vDirection * x_0.adjoint()) +
                          T2_3 * (vDirection * vDirection.adjoint())) - M2;

                    A.block(3 * idx, 3 * idx1, 3, 3) = A.block(3 * idx, 3 * idx1, 3, 3) + M1;
                    A.block(3 * idx, 3 * (idx1 + 1), 3, 3) = A.block(3 * idx, 3 * (idx1 + 1), 3, 3) + M2;

                    idx1++;
                }
                idx1++;
            }
            idx++;
        }
    }
    ViscousForce = -A.llt().solve(VelocityVec) * 8 * M_PI * viscosity;
}


double RegularizedStokeslet::getViscousForceNorm() {
    return ViscousForce.norm();
}


void RegularizedStokeslet::computeFrs() {
    int forceIdx = 0;
    for (int n = 0; n < numRod; n++) {
        for (int i = 0; i < nv; i++) {
            for (int j = 0; j < 3; j++) {
                stepper->addForce(4 * i + j, -ViscousForce(forceIdx), n);
                forceIdx++;
            }
        }
    }
}


void RegularizedStokeslet::computeTCoeffs() {
    T0_1 = (log(edgeLength * R_1 + x_1.dot(vDirection)) - log(edgeLength * R_0 + x_0.dot(vDirection))) / edgeLength;
    T0_3 = -(1 / (R_1 * (edgeLength * R_1 + x_1.dot(vDirection))) -
             1 / (R_0 * (edgeLength * R_0 + x_0.dot(vDirection))));
    T1_1 = (R_1 / (edgeLength * edgeLength) - R_0 / (edgeLength * edgeLength)) -
           T0_1 * x_0.dot(vDirection) / (edgeLength * edgeLength);
    T1_3 = -(1 / (R_1 * edgeLength * edgeLength) - 1 / (R_0 * edgeLength * edgeLength)) -
           T0_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
    T2_3 = -(1 / (R_1 * edgeLength * edgeLength) - 0 / (R_0 * edgeLength * edgeLength)) +
           T0_1 / (edgeLength * edgeLength) - T1_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
    T3_3 = -(1 / (R_1 * edgeLength * edgeLength) - 0 / (R_0 * edgeLength * edgeLength)) +
           2 * T1_1 / (edgeLength * edgeLength) - T2_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
}
