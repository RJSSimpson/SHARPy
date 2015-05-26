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
using namespace Eigen;

void genW(const ::VectorXd& zeta,
		  const int M,
		  const int N,
		  const ::MatrixXd& W_) {
	/**@brief Populate W, downwash interpolation, matrix.
	 * @param zeta Vector of lattice vertices.
	 * @param W downwash interpolation matrix.
	 */

	// const cast Eigen outputs
	MatrixXd& W = const_cast<MatrixXd&> (W_);

	// check relative size of zeta and W
	assert(zeta.size()==W.cols());

	// temporary vars
	VectorXd normals(3*W.rows());
	MatrixXd norMat(W.rows(),3*W.rows());
	MatrixXd Xi(3*W.rows(),W.cols());
	Matrix3d XiKern;

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

void getNormals(const VectorXd& zeta,
		        const int M,
		        const int N,
		        const VectorXd& normals_) {
	/**@brief Get vector containing panel normal vectors.
	 * @param zeta Vector containing lattice vertices.
	 * @param M chordwise panels.
	 * @param N spanwise panels.
	 * @param normals Vector containing normal vectors.
	 */

	// const cast Eigen outputs
	VectorXd& normals = const_cast<VectorXd&> (normals_);

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
		normals.block<3,1>(3*k,0)=Vector3d(n[0],n[1],n[2]);
	}
	return;
}

void genAstar(const VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const MatrixXd& Astar_) {
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
	MatrixXd& Astar = const_cast<MatrixXd&> (Astar_);

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
			Astar.block<3,1>(3*q,k)=Vector3d(sumVel[1],
													sumVel[2],
													sumVel[3]);
		}
	}
	return;
}

double fGeom(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  const double* n) {
	/**@brief Calculate geometry influence function f.
	 * @param r0 Vector r0.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @param n Normal vector of target panel.
	 */

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSum[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);

	// get x
	MulTriad(u,1.0/pow(NormTriad(u),2.0),x);

	// get normalised triads
	NormaliseTriad(r1,r1hat);
	NormaliseTriad(r2,r2hat);

	// sum them
	AddTriad(r1hat,r2hat,normSum);

	return DotTriad(r0,normSum)*DotTriad(x,n);
}

void df_dgeom(const double* r0,
		        const double* r1,
			 	const double* r2,
			 	const double* n,
			 	const Vector3d& f_r0_,
			 	const Vector3d& f_r1_,
			 	const Vector3d& f_r2_,
			 	const Vector3d& f_n_) {
	/**@brief Calculate df/dr and df/dn, row vectors.
	 * @param r0 Vector r0.
     * @param r1 Vector r1.
     * @param r2 Vector r2.
     * @param n Normal vector of target panel.
     * @param f_r0 Row vector, gradient of f w.r.t r0.
     * @param f_r1 Row vector, gradient of f w.r.t r1.
     * @param f_r2 Row vector, gradient of f w.r.t r2.
     * @param f_n Row vector, gradient of f w.r.t n.
     */

	// const cast Eigen outputs
	Vector3d& f_r0 = const_cast<Vector3d&> (f_r0_);
	Vector3d& f_r1 = const_cast<Vector3d&> (f_r1_);
	Vector3d& f_r2 = const_cast<Vector3d&> (f_r2_);
	Vector3d& f_n = const_cast<Vector3d&> (f_n_);

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSum[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);

	// get x
	MulTriad(u,1.0/pow(NormTriad(u),2.0),x);
	Vector3d xV(x[0],x[1],x[2]);

	// get normalised triads
	NormaliseTriad(r1,r1hat);
	NormaliseTriad(r2,r2hat);

	// sum them
	AddTriad(r1hat,r2hat,normSum);

	// if within cut-off radius, r-vector derivatives are set to zero
	if (NormTriad(u) <= 1.0e-5) {

		f_r0.setZero();
		f_r1.setZero();
		f_r2.setZero();

	} else {

		// temps
		double f_r0_triad[3] = {0.0,0.0,0.0};

		// get f_r0
		MulTriad(normSum,DotTriad(x,n),f_r0_triad);
		f_r0(0) = f_r0_triad[0];
		f_r0(1) = f_r0_triad[1];
		f_r0(2) = f_r0_triad[2];

		// get f_r1
		Vector3d r0v(r0[0],r0[1],r0[2]);
		Vector3d r1v(r1[0],r1[1],r1[2]);
		Vector3d r2v(r2[0],r2[1],r2[2]);
		Vector3d nV(n[0],n[1],n[2]);
		// calc dr1Hat_dr1 terms
		Vector3d preMul;
		preMul = DotTriad(x,n)*r0v;
		// calc d(u/|u|^2)/dr1 terms
		Vector3d preMul2;
		preMul2 = DotTriad(r0,normSum)*nV;
		f_r1 = preMul.transpose()*dxHat_dx(r1v)
			 + preMul2.transpose()*duHat2_dr1(r1v,r2v);

		// get f_r2
		f_r2 = preMul.transpose()*dxHat_dx(r2v)
			 + preMul2.transpose()*duHat2_dr2(r1v,r2v);
	}

	// get f_n
	f_n = DotTriad(r0,normSum)*xV; // f_n needs to be transposed

	return;
}

Matrix3d dxHat_dx(const Vector3d& x) {
	/**@brief Calculate derivative of x/|x| w.r.t x.
	 * @param x Vector x.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	if (x.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX = (1.0/x.norm())*(Eye - (1.0/pow(x.norm(),2.0)) * x*x.transpose());
	}
	return dX;
}

Matrix3d duHat2_dr1(const Vector3d& r1,
				     const Vector3d& r2) {
	/**@brief Calculate derivative of (r1xr2)/|r1xr2|^2 w.r.t r1.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	// temps
	Vector3d u = r1.cross(r2);

	if (u.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX = - (1.0/pow(u.norm(),2.0))
			 * (Eye - (2.0/pow(u.norm(),2.0)) * u*u.transpose())
			 * skew(r2);
	}
	return dX;
}

Matrix3d duHat2_dr2(const Vector3d& r1,
				     const Vector3d& r2) {
	/**@brief Calculate derivative of (r1xr2)/|r1xr2|^2 w.r.t r2.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	// temps
	Vector3d u = r1.cross(r2);

	if (u.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX =   (1.0/pow(u.norm(),2.0))
			 * (Eye - (2.0/pow(u.norm(),2.0)) * u*u.transpose())
			 * skew(r1);
	}
	return dX;
}

Matrix3d skew(const Vector3d& x) {
	/**@brief Calculate skew-symmetric matrix of vector x.
	 */
	Matrix3d X;
	X  << 0.0,  -x(2), x(1),
		  x(2),  0.0, -x(0),
		 -x(1),  x(0), 0.0;
	return X;
}

Matrix3d dn_dd(const Vector3d& d, const Vector3d& e) {
	/**@brief Calculate derivative of normal vector w.r.t d (diagonal 1).
	 * @param d Diagonal vector 1 (corner 1 to corner 3).
	 * @param e Diagonal vector 2 (corner 4 to corner 2).
	 * @return dX Matrix output.
	 */
	return -dxHat_dx(d.cross(e))*skew(e);
}

Matrix3d dn_de(const Vector3d& d, const Vector3d& e) {
	/**@brief Calculate derivative of normal vector w.r.t e (diagonal 2).
	 * @param d Diagonal vector 1 (corner 1 to corner 3).
	 * @param e Diagonal vector 2 (corner 4 to corner 2).
	 * @return dX Matrix output.
	 */
	return dxHat_dx(d.cross(e))*skew(d);
}

void dAgamma0_dZeta(const VectorXd& zetaSrc,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const VectorXd& gamma0,
					 const VectorXd& zetaTgt,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
					 const MatrixXd& dX_) {
	/**@brief Calculate tensor-free derivative of (A gamma_0) w.r.t zeta.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spnwise panels on source lattice.
	 * @param gamma0 Reference circulation distribution on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spnwise panels on target lattice.
	 * @param dX K x 3K_{\zeta_{src}} matrix output.
	 */

	// const cast Eigen outputs
	MatrixXd& dX = const_cast<MatrixXd&> (dX_);

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int qTgt = 3*(mTgt+1)*(nTgt+1);
	unsigned int ll = 0;
	unsigned int llp1 = 0; //segment counters

	// loop through DoFs to make (1x3) submatrices
	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			for (unsigned int q = 0; q < qTgt; q++) {
				//
				for (unsigned int l = 1; l < 5; l++) {
					// roll around segment index
					if (l < 4) {
						ll = l;
						llp1 = l+1;
					} else if (l == 4) {
						ll = l;
						llp1= 1;
					}
					// contributions from targets

					// contributions from sources
				}
			}
		}
	}
	return;
}
