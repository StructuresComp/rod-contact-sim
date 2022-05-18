#ifndef BASECONTACTPOTENTIAL_H
#define BASECONTACTPOTENTIAL_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"
#include "collisionDetector.h"
#include "symbolicEquations.h"


class baseContactPotential {
public:
    baseContactPotential(vector<elasticRod *> m_rod, timeStepper &m_stepper, collisionDetector &m_col_detector,
                         double m_delta);
    void updateContactStiffness();
    void computeFc();
    void computeFcJc();
    double contact_stiffness;
protected:
    vector<elasticRod *> rod_vec;
    timeStepper *stepper;
    collisionDetector *col_detector;
    symbolicEquations *sym_eqs;
    double h2;
    double delta;
    double scale;
    int nv;

    int idx1;
    int idx2;
    int idx3;
    int idx4;
    ConstraintType constraint_type;

    Vector<double, 14> contact_input;
    Vector<double, 12> contact_gradient;
    Matrix<double, 12, 12> contact_hessian;

    Vector<double, 8> p2p_input;
    Vector<double, 11> e2p_input;
    Vector<double, 14> e2e_input;

    Vector<double, 6> p2p_gradient;
    Vector<double, 9> e2p_gradient;
    Vector<double, 12> e2e_gradient;

    Matrix<double, 6, 6> p2p_hessian;
    Matrix<double, 9, 9> e2p_hessian;
    Matrix<double, 12, 12> e2e_hessian;

    void prepContactInput();
};

#endif