#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "inertialForce.h"

// include external force
#include "RegularizedStokeslet.h"
#include "contactPotentialIPC.h"
#include "contactPotentialIMC.h"
#include "externalGravityForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

#include "collisionDetector.h"

//
// clamp function is copied from BASim code (Columbia Univ.)
/** Clamps scalar to the range [min,max]. */
template<typename T>
inline T clamp(const T &scalar, const T &min, const T &max) {
    if (scalar < min) return min;
    if (scalar > max) return max;
    return scalar;
}

class world {
public:
    world();

    world(setInput &m_inputData);

    ~world();

    void setRodStepper();

    void updateTimeStep();

    int simulationRunning();

    int numPoints();

    double getScaledCoordinate(int i, int j);

    double getCurrentTime();

    double getTotalTime();

    bool isRender();

    // file output
    void OpenFile(ofstream &simfile, ofstream &configfile);

    void CloseFile(ofstream &simfile, ofstream &configfile);

    void CoutData(ofstream &simfile, ofstream &configfile, double &time_taken);

    int numVertices;


private:

    // Physical parameters
    double RodLength;
    double helixRadius, helixPitch;
    double rodRadius;
    double youngM;
    double Poisson;
    double shearM;
    double deltaTime;
    double totalTime;
    double density;
    Vector3d gVector;
    double viscosity;
    double epsilon;
    int numRod;
    double omega;
    double distance;

    double tol, stol;
    int maxIter; // maximum number of iterations
    double characteristicForce;
    double forceTol;

    // Geometry
    MatrixXd vertices;
    double currentTime;
    double mu;
    double nu;
    double delta;
    double col_limit;

    // flag for ipc or imc
    int ipc;
    double Kappa;
    double dHat;

    std::vector<elasticRod *> rodsVector;

    // set up the time stepper
    timeStepper *stepper;
    double *dx;

    std::vector<timeStepper *> stepperVector;
    std::vector<double *> totalForceVector;

    std::vector<elasticStretchingForce *> v_stretchForce;
    std::vector<elasticBendingForce *> v_bendingForce;
    std::vector<elasticTwistingForce *> v_twistingForce;
    std::vector<inertialForce *> v_inertialForce;
    std::vector<externalGravityForce *> v_gravityForce;

    RegularizedStokeslet *m_RegularizedStokeslet;
    collisionDetector *m_collisionDetector;
    contactPotentialIPC *m_contactPotentialIPC;
    contactPotentialIMC *m_contactPotentialIMC;

    int Nstep;
    int timeStep;
    int iter;

    void rodGeometry(double offset_x, double offset_y, double offset_theta);

    void rodBoundaryCondition(int n);

    // Variables about angular velocity
    double deltaTheta;

    bool render; // should the OpenGL rendering be included?
    bool saveData; // should data be written to a file?

    void updateEachRod(int n);

    double axisLength;
    int ne, nv;

    int contactNum;

    double axisLengthInput;
    double deltaLengthInput;

    void newtonMethod(bool &solved);

    int friction;
    int line_search;

    void newtonDamper();

    double lineSearch();

    double IPCLineSearch();

    void updateBoundary();

    void printSimData();

    double alpha;

    double normf;
    double normf0;
};

#endif
