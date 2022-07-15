#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticRod.h"

#include "eigenIncludes.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif


//extern "C" void dgbsv_( int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info );

/* PARDISO prototype. */
// extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
// extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
//                   double *, int    *,    int *, int *,   int *, int *,
//                      int *, double *, double *, int *, double *);
// extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
// extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
// extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
//                            double *, int *);

class timeStepper
{
public:
	timeStepper(std::vector<elasticRod*> m_rod_vec);
	~timeStepper();
	double* getForce();
	double* getJacobian();
	void setZero();
	void addForce(int ind, double p, int idx);
	void addJacobian(int ind1, int ind2, double p, int idx);
  void addJacobian(int ind1, int ind2, double p, int idx1, int idx2);
	void integrator();
	double *dx;
  VectorXd DX;
	VectorXd Force;
  VectorXd force;
  MatrixXd Jacobian;



private:
	elasticRod *rod;
  elasticRod *rod1;

	std::vector<elasticRod*> rod_vec;
	int kl, ku, freeDOF;

	double *totalForce;
	double *jacobian;

	vector<int> freeDOF_vec;

	// utility variables
	int mappedInd, mappedInd1, mappedInd2;
	int row, col, offset;
	int NUMROWS;
	int jacobianLen;
	int nrhs;
  int *ipiv;
  int info;
  int ldb;

	void pardisoSolver();
};

#endif
