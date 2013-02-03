/*
 * vorticity.hpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#ifndef VORTICITY_HPP_
#define VORTICITY_HPP_

#include <Eigen/Dense>

Eigen::Vector3d BiotSegment(const Eigen::Vector3d& xp, \
							const Eigen::Vector3d& x1, \
							const Eigen::Vector3d& x2, \
							const double& gam);

void BiotSegment_map(const Eigen::Map<Eigen::Vector3d>& xp, \
								const Eigen::Map<Eigen::Vector3d>& x1, \
								const Eigen::Map<Eigen::Vector3d>& x2, \
								const double& gam,
								const Eigen::Map<Eigen::Vector3d>& Uind_triad);

void C_BiotSegment(const double* xP,\
					  const double* x1,\
					  const double* x2,\
					  double& Gamma,\
					  double* Uind);



#endif /* VORTICITY_HPP_ */
