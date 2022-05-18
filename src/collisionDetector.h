#ifndef COLLISIONDETECTOR_H
#define COLLISIONDETECTOR_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"


enum ConstraintType {
    PointToPoint=0,
    PointToEdge=1,
    EdgeToEdge=2
};


enum ContactPiecewise {
    Normal=0,
    Penetrated=1
};


class collisionDetector
{
public:
    collisionDetector(std::vector<elasticRod *> m_rod_vec, timeStepper &m_stepper, double m_delta,
                      double m_col_limit);

    void constructCandidateSet();
    void detectCollisions();
    void computeMinDistance(const Vector3d &v1s, const Vector3d &v1e, const Vector3d &v2s, const Vector3d &v2e, double& dist);

    MatrixXi edge_ids;
    MatrixXi contact_ids;
    vector<Vector2i> candidate_set;
    int num_collisions;
    double min_dist;

private:
    vector<elasticRod*> rod_vec;
    elasticRod* rod;
    timeStepper* stepper;
    double delta;
    double col_limit;
    int num_edge_combos;
    double scale;
    int rod_num;
    int nv;

    double contact_limit;
    double candidate_limit;
    double numerical_limit;

    void fixbound(double &x);
    void computeMinDistance(int &idx1, int &idx2, int&idx3, int &idx4, double &dist, ConstraintType& constraintType);
};
#endif
