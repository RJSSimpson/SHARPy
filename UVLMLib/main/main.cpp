/*
 * main.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#include <vorticity.hpp>
#include <Eigen/Dense>
#include <iostream>
#include "time.h"

int main() {
	//test vortex line code segment
	// create points and gamma
	int TotalCalls = 160000;

	Eigen::Vector3d xp(0.0,0.0,-1.0);
	Eigen::Vector3d x1(-0.5,0.0,0.0);
	Eigen::Vector3d x2(0.5,0.0,0.0);
	double gam  = 1.0;
	Eigen::Vector3d Result(0.0,0.0,0.0);
	clock_t start, end, start2, end2;

	start = clock();
	//call function
	for (int i = 0; i < TotalCalls; i++) {
		Result = BiotSegment(xp,x1,x2,gam);
		//std::cout << Result << std::endl;
	}
	end = clock();
	std::cout << "done! (time elapsed: " << 1000*(end-start)/CLOCKS_PER_SEC << "ms)" << std::endl;

	double xp_ptr[] = {0.0,0.0,-1.0};
	double x1_ptr[] = {-0.5,0.0,0.0};
	double x2_ptr[] = {0.5,0.0,0.0};
	double Result_ptr[] = {0.0,0.0,0.0};

	//call function
	start2 = clock();
	for (int i = 0; i < TotalCalls; i++) {
		C_BiotSegment(xp_ptr,x1_ptr,x2_ptr,gam,Result_ptr);
		//std::cout << Result_ptr[0] << std::endl << Result_ptr[1] << std::endl << Result_ptr[2] << std::endl;
	}
	end2 = clock();
	std::cout << "done! (time elapsed: " << 1000*(end2-start2)/CLOCKS_PER_SEC << "ms)" << std::endl;

	//Print result
	std::cout << Result << std::endl;

	return 0;
}
