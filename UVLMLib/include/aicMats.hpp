/*
 * aicMats.hpp
 *
 *  Created on: 13 May 2015
 *      Author: rjs10
 */

#ifndef AICMATS_HPP_
#define AICMATS_HPP_

#include<Eigen/Dense>

void genW(const Eigen::VectorXd& zeta,
		   const int M,
		   const int N,
		   const Eigen::MatrixXd& W_);

void getNormals(const Eigen::VectorXd& zeta,
		          const int M,
		          const int N,
		          const Eigen::VectorXd& normals_);

void genAstar(const Eigen::VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const Eigen::VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const Eigen::MatrixXd& Astar_);

#endif /* AICMATS_HPP_ */
