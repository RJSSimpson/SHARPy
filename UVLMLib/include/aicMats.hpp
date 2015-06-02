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

double fGeom(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  const double* n);

void df_dgeom(const double* r0,
		        const double* r1,
			 	const double* r2,
			 	const double* n,
			 	const Vector3d& f_r0_,
			 	const Vector3d& f_r1_,
			 	const Vector3d& f_r2_,
			 	const Vector3d& f_n_);

Matrix3d dxHat_dx(const Vector3d& x);

Matrix3d duHat2_dr1(const Vector3d& r1,
				     const Vector3d& r2);

Matrix3d duHat2_dr2(const Vector3d& r1,
				     const Vector3d& r2);

Matrix3d skew(const Vector3d& x);

Matrix3d dn_dd(const Vector3d& d, const Vector3d& e);

Matrix3d dn_de(const Vector3d& d, const Vector3d& e);

void dAgamma0_dZeta(const VectorXd& zetaSrc,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const VectorXd& gamma0,
					 const VectorXd& zetaTgt,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
					 const MatrixXd& dX_);

void AIC(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  double* dX_);

#endif /* AICMATS_HPP_ */
