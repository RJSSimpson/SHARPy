/**@brief      AIC and downwash matrices from new UVLM description.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       13/05/2015
 * @pre        None
 * @warning    None
 */

#include <aicMats.hpp>
#include <assert.h>
#include <indices.hpp>
#include <PanelTools.hpp>
#include <iostream>
#include <triads.hpp>
#include <stdexcept>
#include <vorticity.hpp>

void genW(const Eigen::VectorXd& zeta,
		  const int M,
		  const int N,
		  const Eigen::MatrixXd& W_) {
	/**@brief Populate W, downwash interpolation, matrix.
	 * @param zeta Vector of lattice vertices.
	 * @param W downwash interpolation matrix.
	 */

	// const cast Eigen outputs
	Eigen::MatrixXd& W = const_cast<Eigen::MatrixXd&> (W_);

	// check relative size of zeta and W
	assert(zeta.size()==W.cols());

	// temporary vars
	Eigen::VectorXd normals(3*W.rows());
	Eigen::MatrixXd norMat(W.rows(),3*W.rows());
	Eigen::MatrixXd Xi(3*W.rows(),W.cols());
	Eigen::Matrix3d XiKern;

	// calculate normal vectors
	getNormals(zeta,M,N,normals);

	// get matrix of normal vectors
	unsigned int K=M*N;
	for (unsigned int k = 0; k < K; k++){
		norMat.block<1,3>(k,3*k)=normals.block<3,1>(3*k,0).transpose();
	}

	// get interpolation matrix, \Xi
	for (unsigned int k = 0; k < K; k++){
		for (unsigned int q = 0; q < W.cols()/3; q++){
			XiKernel(k,q,N,0.5,0.5,XiKern);
			Xi.block<3,3>(3*k,3*q)=XiKern;
		}
	}

	// calculate product of normal and interpolation matrices
	W = norMat*Xi;
	return;
}

void getNormals(const Eigen::VectorXd& zeta,
		        const int M,
		        const int N,
		        const Eigen::VectorXd& normals_) {
	/**@brief Get vector containing panel normal vectors.
	 * @param zeta Vector containing lattice vertices.
	 * @param M chordwise panels.
	 * @param N spanwise panels.
	 * @param normals Vector containing normal vectors.
	 */

	// const cast Eigen outputs
	Eigen::VectorXd& normals = const_cast<Eigen::VectorXd&> (normals_);

	// check inputs
	assert(zeta.size()==3*(M+1)*(N+1));
	assert(normals.size()==3*M*N);

	// declare corner points
	double x1[3] = {0.0,0.0,0.0};
	double x2[3] = {0.0,0.0,0.0};
	double x3[3] = {0.0,0.0,0.0};
	double x4[3] = {0.0,0.0,0.0};
	// declare normal vectors
	double n[3] = {0.0,0.0,0.0};

	unsigned int K = M*N;
	for (unsigned int k = 0; k < K; k++){

		// get corner points
		AssignTriad(x1,zeta(3*q_k(k,N,1)),zeta(3*q_k(k,N,1)+1),zeta(3*q_k(k,N,1)+2));
		AssignTriad(x2,zeta(3*q_k(k,N,2)),zeta(3*q_k(k,N,2)+1),zeta(3*q_k(k,N,2)+2));
		AssignTriad(x3,zeta(3*q_k(k,N,3)),zeta(3*q_k(k,N,3)+1),zeta(3*q_k(k,N,3)+2));
		AssignTriad(x4,zeta(3*q_k(k,N,4)),zeta(3*q_k(k,N,4)+1),zeta(3*q_k(k,N,4)+2));

		// calculate normal
		PanelNormal(x1,x2,x3,x4,n);

		//save to normals
		normals.block<3,1>(3*k,0)=Eigen::Vector3d(n[0],n[1],n[2]);
	}
	return;
}

void genAstar(const Eigen::VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const Eigen::VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const Eigen::MatrixXd& Astar_) {
	/**@brief Get vector containing panel normal vectors.
	 * @param zetaSrc Vector containing lattice vertices of source.
	 * @param mSrc Spanwise panels in the source lattice.
	 * @param nSrc Spanwise panels in the source lattice.
	 * @param zetaTgt Vector containing lattice vertices of target.
	 * @param mTgt Spanwise panels in the source lattice.
	 * @param nTgt Spanwise panels in the target lattice.
	 * @param N spanwise panels.
	 * @param Astar AIC matrix from source to target.
	 */

	// const cast Eigen outputs
	Eigen::MatrixXd& Astar = const_cast<Eigen::MatrixXd&> (Astar_);

	// reused panel numbers
	const unsigned int kSrc = mSrc*nSrc;
	const unsigned int kTgt_zeta = (mTgt+1)*(nTgt+1);

	// check sizes
	assert(zetaSrc.size()==3*(mSrc+1)*(nSrc+1));
	assert(zetaTgt.size()==3*kTgt_zeta);
	assert(Astar.rows()==3*kTgt_zeta);
	assert(Astar.cols()==kSrc);

	// init temps
	const double gamma = 1.0;
	unsigned int segStart = 0;
	unsigned int segEnd = 0;
	double xP[3] = {0.0, 0.0, 0.0};
	double x1[3] = {0.0, 0.0, 0.0};
	double x2[3] = {0.0, 0.0, 0.0};
	double vel[3] = {0.0, 0.0, 0.0};
	double sumVel[3] = {0.0, 0.0, 0.0};

	// loop through grid points, panels, segments
	for (unsigned int q = 0; q < kTgt_zeta; q++) {
		for (unsigned int k = 0; k < kSrc; k++) {
			// reset sum
			AssignTriad(sumVel,0.0,0.0,0.0);
			for (unsigned int l = 1; l < 5; l++) {
				// get segment corner points
				switch(l){
				case 1 :
					segStart=1;
					segEnd=2;
					break;
				case 2 :
					segStart=2;
					segEnd=3;
					break;
				case 3 :
					segStart=3;
					segEnd=4;
					break;
				case 4 :
					segStart=4;
					segEnd=1;
					break;
				default:
					throw std::invalid_argument("Segment index invalid.");
				}
				// get target point
				AssignTriad(xP,zetaTgt(3*q),zetaTgt(3*q+1),zetaTgt(3*q+2));
				// start of segment
				AssignTriad(x1,
						    zetaSrc(3*q_k(k,nSrc,segStart)),
						    zetaSrc(3*q_k(k,nSrc,segStart)+1),
						    zetaSrc(3*q_k(k,nSrc,segStart)+2));
				// end of segment
				AssignTriad(x2,
							zetaSrc(3*q_k(k,nSrc,segEnd)),
							zetaSrc(3*q_k(k,nSrc,segEnd)+1),
							zetaSrc(3*q_k(k,nSrc,segEnd)+2));
				// get velocity
				C_BiotSegment(xP,x1,x2,gamma,vel);
				// add to total
				AddTriad(sumVel,vel,sumVel);
			}
			// save to aStar
			Astar.block<3,1>(3*q,k)=Eigen::Vector3d(sumVel[1],
													sumVel[2],
													sumVel[3]);
		}
	}
	return;
}
