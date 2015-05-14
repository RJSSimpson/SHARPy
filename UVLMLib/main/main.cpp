/*
 * main.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#include <aicMats.hpp>
#include <Eigen/Dense>
#include <iostream>

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

	std::cout << Astar << std::endl;

	A=W*Astar;

	std::cout << std::endl << A << std::endl;
}
