#ifndef CONTACTPOTENTIALIPC_H
#define CONTACTPOTENTIALIPC_H

#include "baseContactPotential.h"

class contactPotentialIPC : public baseContactPotential {
public:
    contactPotentialIPC(std::vector<elasticRod *> m_rod_vec, timeStepper &m_stepper, collisionDetector &m_col_detector,
                        double m_delta);
    ~contactPotentialIPC();
    void computeFc();
    void computeFcJc();
    double computeUpperBound(const VectorXd &p);

private:
    bool LQF_CCD3(Vector3d &a0s, Vector3d &a1s, Vector3d &b0s, Vector3d &b1s, Vector3d &p1, Vector3d &p2, Vector3d &p3,
                  Vector3d &p4, double &toi, double OrigLen = 1);
};

#endif
