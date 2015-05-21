/*
 * main.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#include <aicMats.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <triads.hpp>
using namespace std;
using namespace Eigen;

double fRand(double fMin, double fMax) {
	// generate random numbers
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main() {
	// Test UVLMLib functions during development.

	// grid vars
	const int M=1;
	const int N=1;
	const int K=M*N;
	const int K_zeta = (M+1)*(N+1);

	Eigen::VectorXd zeta(3*K_zeta);
	Eigen::MatrixXd W(K,3*K_zeta);
	Eigen::MatrixXd Astar(3*K_zeta,K);
	Eigen::MatrixXd A(K,K);

//	zeta << 0.0, 0.0, 0.0,
//			0.0, 1.0, 0.0,
//			0.0, 2.0, 0.0,
//			1.0, 0.0, 0.0,
//			1.0, 1.0, 0.0,
//			1.0, 2.0, 0.0;

	zeta << 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 0.0, 0.0,
			1.0, 1.0, 0.0;


	genW(zeta,M,N,W);

	genAstar(zeta,M,N,zeta,M,N,Astar);

	// test dAgam_dzeta subroutines

	// dxHat_dx
	Matrix3d dX;
	Vector3d x(1.0,1.0,1.0);
	srand(time(NULL));
	Vector3d dx = Vector3d::Random()/(10.0*x.norm());
	Vector3d f = (x + dx).normalized();
	dX = dxHat_dx(x);
	Vector3d fApprox = x.normalized() + dX * dx;
	Vector3d diff = f-fApprox;
	cout << "dxHat_dx: ---------------" << endl << "dx:" << endl << dx
	     << endl << "norm (l-2) rel error:" << endl << diff.norm()/dx.norm()
	     << endl;

	// Init matrices for df_dGeom tests
	Vector3d f_r0, f_r1, f_r2, f_n;
	//df_dr0
	double r0[3] = {1.0,0.0,0.0};
	double r1[3] = {1.0,1.0,0.0};
	double r2[3] = {0.0,1.0,0.0};
	double n[3] = {0.0,0.0,1.0};
	// calc variations
	df_dgeom(r0,r1,r2,n,f_r0,f_r1,f_r2,f_n);
	double dr0[3] = {fRand(-1,1)/(10*NormTriad(r0)),
					  fRand(-1,1)/(10*NormTriad(r0)),
					  fRand(-1,1)/(10*NormTriad(r0))};
	cout << "df_dr0: ---------------" << endl << "dr0:"
		 << endl << Vector3d(dr0[0],dr0[1],dr0[2])
		 << endl;
	double r0Test[3];
	AddTriad(r0,dr0,r0Test);
	double fGeomExact = fGeom(r0Test,r1,r2,n);
	double fGeomAppDr0 = fGeom(r0,r1,r2,n) + f_r0.transpose()*Vector3d(dr0[0],dr0[1],dr0[2]);
	cout << "error:" << endl << fGeomExact - fGeomAppDr0
	     << endl;
}
