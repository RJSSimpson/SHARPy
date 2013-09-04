/*
 * PanelTools.cpp
 *
 *  Created on: 13 Feb 2013
 *      Author: rjs10
 */

#include <PanelTools.hpp>
#include <triads.hpp>

void PanelNormal(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* normal) {
	/** @brief Calculate panel normal based on corner points
	 *  @details Method from Katz and Plotkin.
	 */

	// Initialise temp triads
	double A[3];
	double B[3];
	double N[3];

	// define A and B triads
	SubTriad(p3,p4,A);
	SubTriad(p2,p4,B);

	// calculate cross product
	CrossTriad(A,B,N);

	//normalise and copy to normal
	NormaliseTriad(N,normal);
}

void PanelTau_c(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* Tau_c) {
	/** @brief Calculate panel chordwise unit vector
	 */

	double side1[3];
	double side2[3];

	// calculate panel side vectors
	SubTriad(p3,p2,side1);
	SubTriad(p4,p1,side2);

	//combine
	AddTriad(side1,side2,side1);

	//normalise and write to Tau_c
	NormaliseTriad(side1,Tau_c);
}

void PanelTau_s(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* Tau_s) {
	/** @brief Calculate panel sppanwise unit vector
	 */

	double side1[3];
	double side2[3];

	// calculate panel side vectors
	SubTriad(p2,p1,side1);
	SubTriad(p3,p4,side2);

	//combine
	AddTriad(side1,side2,side1);

	//normalise and write to Tau_c
	NormaliseTriad(side1,Tau_s);
}


double PanelDeltaC(const double* p1, const double* p2, \
				 const double* p3, const double* p4) {
	/** @brief Calculate panel chordwise extent.
	 */
	double side1[3];
	double side2[3];

	// calculate panel side vectors
	SubTriad(p4,p1,side1);
	SubTriad(p3,p2,side2);

	// combine to find mean
	AddTriad(side1,side2,side1);
	DivTriad(side1,2.0,side1);

	return NormTriad(side1);
}

double PanelDeltaS(const double* p1, const double* p2, \
				 const double* p3, const double* p4) {
	/** @brief Calculate panel spanwise extent.
	 */
	double side1[3];
	double side2[3];

	// calculate panel side vectors
	SubTriad(p2,p1,side1);
	SubTriad(p3,p4,side2);

	// combine to find mean
	AddTriad(side1,side2,side1);
	DivTriad(side1,2.0,side1);

	return NormTriad(side1);
}


double PanelArea(const double* p1, const double* p2, \
				 const double* p3, const double* p4) {
	/** @brief Calculate panel area.
	 */
	double diag1[3];
	double diag2[3];
	double cross[3];

	// calculate panel side vectors
	SubTriad(p4,p2,diag1);
	SubTriad(p3,p1,diag2);

	// find cross product
	CrossTriad(diag1,diag2,cross);

	return 0.5*NormTriad(cross);
}


