/*
 * wrapper.cpp
 *
 *  Created on: 29 Jan 2013
 *      Author: rjs10
 */

#include "vorticity.hpp"
#include <Eigen/Dense>
#include <iostream>


// forward declare some useful functions

Eigen::Map<Eigen::Vector3d> C2Eig_triad_map(double* x_triad) {
	return Eigen::Map<Eigen::Vector3d>(x_triad);
}

Eigen::Vector3d C2Eig_triad(const double* x_triad) {
	return Eigen::Vector3d(x_triad[0], x_triad[1], x_triad[2]);
}

void Eig2C_triad(const Eigen::Vector3d& Eig, double* x_triad) {
	x_triad[0] = Eig[0];
	x_triad[1] = Eig[1];
	x_triad[2] = Eig[2];
}

//wrapper interface must by in C-syntax
extern "C" {

void cpp_wrap_vorticity_biotsegment(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for BiotSegment.
	 *  @details c-style arrays are converted into Eigen::Vector3d.
	 */

	//call BiotSegment
	Eigen::Vector3d Uind = BiotSegment(C2Eig_triad(xp_triad),\
									   C2Eig_triad(x1_triad),\
									   C2Eig_triad(x2_triad),\
									   Gamma);

	//overwrite Uind with result
	Eig2C_triad(Uind,Uind_triad);

	return;
}

void cpp_wrap_vorticity_biotsegment_map(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for BiotSegment.
	 *  @details c-style arrays are converted into Eigen::Vector3d.
	 */

	//call BiotSegment
	BiotSegment_map(C2Eig_triad_map(xp_triad),\
					C2Eig_triad_map(x1_triad),\
					C2Eig_triad_map(x2_triad),\
					Gamma,\
					C2Eig_triad_map(Uind_triad) );

	//overwrite Uind with result
	//Eig2C_triad(Uind,Uind_triad);

	return;
}

void c_wrap_vorticity_biotsegment(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for C_BiotSegment.
	  * @details c-style pointers only.
	  */
	C_BiotSegment(xp_triad, x1_triad, x2_triad,\
				  Gamma, Uind_triad);
}

void cpp_wrap_test_biotsegment(const int& NumTests) {
	Eigen::Vector3d xp(0.0,0.0,-1.0);
	Eigen::Vector3d x1(-0.5,0.0,0.0);
	Eigen::Vector3d x2(0.5,0.0,0.0);
	double gam  = 1.0;
	Eigen::Vector3d Result(0.0,0.0,0.0);

	//call function
	for (int i = 0; i < NumTests; i++) {
		Result = BiotSegment(xp,x1,x2,gam);
	}
	return;
}

void c_wrap_test_biotsegment(const int& NumTests) {
	double xp[] = {0.0,0.0,-1.0};
	double x1[] = {-0.5,0.0,0.0};
	double x2[] = {0.5,0.0,0.0};
	double gam  = 1.0;
	double Result[] = {0.0,0.0,0.0};

	//call function
	for (int i = 0; i < NumTests; i++) {
		C_BiotSegment(xp,x1,x2,gam,Result);
	}
	return;
}

} // END extern C


