#include "inertialForce.h"

inertialForce::inertialForce(elasticRod &m_rod, timeStepper &m_stepper, int idx)
{
	rod = &m_rod;
	stepper = &m_stepper;

	ForceVec = VectorXd::Zero(rod->ndof);
	rod_idx = idx;
}

inertialForce::~inertialForce()
{
	;
}

void inertialForce::computeFi()
{
	ForceVec = VectorXd::Zero(rod->ndof);

	for (int i=0; i<rod->ndof; i++)
	{
		f = rod->massArray[i] * (rod->x[i] - rod->x0[i]) / ((rod->dt) *(rod->dt))
				- (rod->massArray[i] * rod->u[i])/(rod->dt);

		ForceVec(i) = f;

		stepper->addForce(i, f, rod_idx);
	}
}

void inertialForce::computeJi()
{
	for (int i=0; i<rod->ndof; i++)
    {
		jac = rod->massArray(i) / ((rod->dt) *(rod->dt));
		stepper->addJacobian(i, i, jac, rod_idx);
	}
}
