/*
 * aicMats.hpp
 *
 *  Created on: 13 May 2015
 *      Author: rjs10
 */

#ifndef AICMATS_HPP_
#define AICMATS_HPP_

#include<Eigen/Dense>
using namespace Eigen;

void genW(const VectorXd& zeta,
		   const int M,
		   const int N,
		   const MatrixXd& W_);

void getNormals(const VectorXd& zeta,
		          const int M,
		          const int N,
		          const VectorXd& normals_);

void genAstar(const VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const MatrixXd& Astar_);

void dxHat_dx(const Vector3d& x, const Matrix3d& dX_);

#endif /* AICMATS_HPP_ */
