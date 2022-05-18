#ifndef CONTACTPOTENTIALIMC_H
#define CONTACTPOTENTIALIMC_H

#include "baseContactPotential.h"


enum FrictionType {
    ZeroVel = 0,
    Sliding = 1,
    Sticking = 2
};


class contactPotentialIMC : public baseContactPotential {
public:
    contactPotentialIMC(vector<elasticRod *> m_rod, timeStepper &m_stepper, collisionDetector &m_col_detector,
                        double m_delta, double m_mu, double m_nu);
    void computeFc();
    void computeFcJc();

private:
    bool friction;
    double mu;
    double nu;
    double K1;
    double K2;

    Vector<double, 39> friction_input;
    Vector<double, 12> friction_forces;
    Matrix<double, 12, 12> friction_partials_dfr_dx;
    Matrix<double, 12, 12> friction_partials_dfr_dfc;
    Matrix<double, 12, 12> friction_jacobian;
    Matrix<double, 3, 12> friction_zero_matrix;

    ContactPiecewise contact_type;
    FrictionType friction_type;

    void prepFrictionInput();

    void computeFriction();
};

#endif