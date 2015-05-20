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
using namespace std;
using namespace Eigen;

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
	dxHat_dx(x,dX);
	Vector3d fApprox = x.normalized() + dX * dx;
	Vector3d diff = f-fApprox;
	cout << "dxHat_dx: ---------------" << endl << "dx:" << endl << dx
	     << endl << "norm (l-2) rel error:" << endl << diff.norm()/dx.norm()
	     << endl;
}
