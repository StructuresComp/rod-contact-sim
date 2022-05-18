#ifndef REGULARIZEDSTOKESLET_H
#define REGULARIZEDSTOKESLET_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class RegularizedStokeslet {
public:
    RegularizedStokeslet(vector<elasticRod *> m_rod_vec, timeStepper &m_stepper, double m_viscosity, double m_epsilon);

    ~RegularizedStokeslet();

    void prepareForViscousForce();

    void computeFrs();

    void computeJrs();

    double getViscousForceNorm();

    VectorXd ViscousForce;

private:

    vector<elasticRod *> rod_vec;
    elasticRod *rod_1;
    elasticRod *rod_2;

    int numRod;
    int nv;

    timeStepper *stepper;

    double viscosity;
    double epsilon;

    Matrix3d Id3;

    Vector3d uPos, y_0, y_1;
    Vector3d x_0, x_1;
    Vector3d vDirection;
    double edgeLength;
    double R_0, R_1;

    Vector3d uVelocity;

    double T0_1, T0_3, T1_1, T1_3, T2_3, T3_3;

    void computeTCoeffs();

    Matrix3d M1;
    Matrix3d M2;

    MatrixXd A;
    VectorXd VelocityVec;
};

#endif
