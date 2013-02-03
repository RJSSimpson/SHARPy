/**
 * vorticity.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */
#include <vorticity.hpp>
#include <triads.hpp>
#include <cmath>

Eigen::Vector3d BiotSegment(const Eigen::Vector3d& xp, \
							const Eigen::Vector3d& x1, \
							const Eigen::Vector3d& x2, \
							const double& gam) {
    /** @brief Calculate the velocity induced at point p by a vortex line.
    * with start and end points 1 and 2.
    * Note: if r is less than the vortex core radius then u = v = w = 0
    * @param xp the coordinates of point p
    * @param x1 the coordinates of point 1
    * @param x2 the coordinates of point 2
    * @param gam the circulation strength of the vortex line segment
    * @return u/v/w the induced velocity in x/y/z-direction
    */


    if (gam == 0) {return Eigen::Vector3d(0,0,0);}


    //define position vectors
    Eigen::Vector3d r0 = (x2-x1);
    Eigen::Vector3d r1 = (xp-x1);
    Eigen::Vector3d r2 = (xp-x2);
    Eigen::Vector3d rx = r1.cross(r2);


    // if the point P is within the core radius of p1 p2 or p1Xp2 return zeros
    if((r1.norm() <= 1.0e-5) || (r2.norm() <= 1.0e-5) || (rx.norm() <= 1.0e-5)){
        return Eigen::Vector3d(0,0,0);
    }


    //vector dot products
    double r0d1 = r0.dot(r1);
    double r0d2 = r0.dot(r2);

    //K defined as in Katz & Plotkin
    double K = (gam/(4*M_PI*pow(rx.norm(),2))) *
                    ((r0d1/r1.norm())-(r0d2/r2.norm()));

    return Eigen::Vector3d( K*rx(0), K*rx(1), K*rx(2) );
}

void BiotSegment_map(const Eigen::Map<Eigen::Vector3d>& xp, \
								const Eigen::Map<Eigen::Vector3d>& x1, \
								const Eigen::Map<Eigen::Vector3d>& x2, \
								const double& gam,\
								const Eigen::Map<Eigen::Vector3d>& Uind) {
    /** @brief Calculate the velocity induced at point p by a vortex line.
    * with start and end points 1 and 2.
    * Note: if r is less than the vortex core radius then u = v = w = 0
    * @param xp the coordinates of point p
    * @param x1 the coordinates of point 1
    * @param x2 the coordinates of point 2
    * @param gam the circulation strength of the vortex line segment
    * @return u/v/w the induced velocity in x/y/z-direction
    */

    if (gam == 0.0) {
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = 0.0;
    	return;
    }


    //define position vectors
    Eigen::Vector3d r0 = (x2-x1);
    Eigen::Vector3d r1 = (xp-x1);
    Eigen::Vector3d r2 = (xp-x2);
    Eigen::Vector3d rx = r1.cross(r2);


    // if the point P is within the core radius of p1 p2 or p1Xp2 return zeros
    if((r1.norm() <= 1.0e-5) || (r2.norm() <= 1.0e-5) || (rx.norm() <= 1.0e-5)){
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = 0.0;
    	return;
    }


    //vector dot products
    double r0d1 = r0.dot(r1);
    double r0d2 = r0.dot(r2);

    //K defined as in Katz & Plotkin
    double K = (gam/(4*M_PI*pow(rx.norm(),2))) *
                    ((r0d1/r1.norm())-(r0d2/r2.norm()));


	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = K*rx[0];
	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = K*rx[1];
	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = K*rx[2];
    return;
}

void C_BiotSegment(const double* xP,\
					  const double* x1,\
					  const double* x2,\
					  double& Gamma,\
					  double* Uind) {
	/** @brief Calculate the velocity induced at point p by a vortex line.
	  * with start and end points 1 and 2.
	  * Note: if r is less than the vortex core radius then u = v = w = 0
	  * @param xp the coordinates of point p
	  * @param x1 the coordinates of point 1
	  * @param x2 the coordinates of point 2
	  * @param gam the circulation strength of the vortex line segment
	  * @param Uind u/v/w the induced velocity in x/y/z-direction
	  * @details Solved using pointers and dereferencing only.
	  */

	if (Gamma == 0.0) {
		Uind[0] = 0.0;
		Uind[1] = 0.0;
		Uind[2] = 0.0;
		return;
	}

	double r0[] = {0.0,0.0,0.0};
	double r1[] = {0.0,0.0,0.0};
	double r2[] = {0.0,0.0,0.0};
	double rx[] = {0.0,0.0,0.0};

	SubTriad(x2,x1,r0);
	SubTriad(xP,x1,r1);
	SubTriad(xP,x2,r2);
	CrossTriad(r1,r2,rx);

	if((NormTriad(r1) <= 1.0e-5) || (NormTriad(r2) <= 1.0e-5) || \
			(NormTriad(rx) <= 1.0e-5)) {
		Uind[0] = 0.0;
		Uind[1] = 0.0;
		Uind[2] = 0.0;
	    return;
	}


	double r0d1 = DotTriad(r0,r1);
	double r0d2 = DotTriad(r0,r2);

	double K = (Gamma/(4*M_PI*pow(NormTriad(rx),2))) *
	                    ((r0d1/NormTriad(r1)-(r0d2/NormTriad(r2))));

	//overwrite elements of Uind
	Uind[0] = K*rx[0];
	Uind[1] = K*rx[1];
	Uind[2] = K*rx[2];
}
