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
using namespace std;

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
	double normSub[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);
	if (NormTriad(u) <= 1.0e-5) {
		return 0.0;
	} else {
		// get x
		MulTriad(u,1.0/pow(NormTriad(u),2.0),x);

		// get normalised triads
		NormaliseTriad(r1,r1hat);
		NormaliseTriad(r2,r2hat);

		// sum them
		SubTriad(r1hat,r2hat,normSub);

		return DotTriad(r0,normSub)*DotTriad(x,n);
	}
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
	double normSub[3] = {0.0,0.0,0.0};
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
	SubTriad(r1hat,r2hat,normSub);

	// if within cut-off radius, r-vector derivatives are set to zero
	if (NormTriad(u) <= 1.0e-5) {

		f_r0.setZero();
		f_r1.setZero();
		f_r2.setZero();

	} else {

		// temps
		double f_r0_triad[3] = {0.0,0.0,0.0};

		// get f_r0
		MulTriad(normSub,DotTriad(x,n),f_r0_triad);
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
		preMul2 = DotTriad(r0,normSub)*nV;
		f_r1 = preMul.transpose()*dxHat_dx(r1v)
			 + preMul2.transpose()*duHat2_dr1(r1v,r2v);

		// get f_r2
		f_r2 = -preMul.transpose()*dxHat_dx(r2v)
			 + preMul2.transpose()*duHat2_dr2(r1v,r2v);
	}

	// get f_n
	f_n = DotTriad(r0,normSub)*xV; // f_n needs to be transposed

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

void dAgamma0_dZeta(const double* zetaSrc_,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const double* gamma0_,
					 const double* zetaTgt_,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
					 double* dX_) {
	/**@brief Calculate tensor-free derivative of (A gamma_0) w.r.t zeta.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param gamma0 Reference circulation distribution on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @param dX K x 3K_{\zeta_{src}} matrix output.
	 * @warning If the zeta arguments are the same they must be the same object,
	 * therefore in the function call the arguments must be previously
	 * instantiated objects, e.g (zeta+delZeta, ..., zeta+delZeta, ...) is
	 * invalid because a two temps are instantiated.
	 */

	// eigen map output matrix
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd gamma0(gamma0_,mSrc*nSrc);
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,mTgt*nTgt,3*(mSrc+1)*(nSrc+1));

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int qTgt = 3*(mTgt+1)*(nTgt+1);
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1, c2, c3, c4, cp; // panel corner points, collocation point
	Vector3d d, e, n; // panel diagonal and normal vectors
	Matrix3d n_d, n_e;
	Vector3d r0, r1, r2; // Biot-Savart kernel vectors
	Vector3d f_r0, f_r1, f_r2, f_n; //dvtvs required for target/source variation
	Matrix3d Xi; // interpolating matrix
	double a; // prefactor (\gamma_0(k2)/(4 \pi)).

	// loop through DoFs to make (1x3) submatrices
	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, dn_dd, dn_de, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// diagonals
		d = c3 - c1;
		e = c2 - c4;
		// calc n
		n = (d.cross(e)).normalized();
		// calc dn_dd, dn_de
		n_d = dn_dd(d,e);
		n_e = dn_de(d,e);
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			a = gamma0(k2)/(4*M_PI);
			for (unsigned int q = 0; q < qTgt; q++) {
				for (unsigned int l = 1; l < 5; l++) {
					// roll around segment index
					if (l < 4) {
						ll = l;
						llp1 = l+1;
					} else if (l == 4) {
						ll = l;
						llp1= 1;
					}
					// check if either k1 or k2 are active based on q
					if (q_k(k1,nTgt,1) == q || q_k(k1,nTgt,2) == q ||
						q_k(k1,nTgt,3) == q || q_k(k1,nTgt,4) == q ||
						q_k(k2,nSrc,ll) == q || q_k(k2,nSrc,llp1) == q) {

						// contributions from targets
						// calc r0, r1, r2
						r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
							-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
						// r1
						r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
						// r2
						r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);

						// calc f_r0, f_r1, f_r2, f_n
						df_dgeom(r0.data(),
								 r1.data(),
								 r2.data(),
								 n.data(),
								 f_r0,
								 f_r1,
								 f_r2,
								 f_n);

						// add effect to output matrix
						XiKernel(k1,q,nTgt,0.5,0.5,Xi);

						if (q_k(k1,nTgt,1) == q) {
							// add Xi, and Kronecker delta terms, c1
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   - f_n.transpose()*n_d);
						} else if (q_k(k1,nTgt,2) == q) {
							// add Xi, and Kronecker delta terms, c2
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   + f_n.transpose()*n_e);
						} else if (q_k(k1,nTgt,3) == q) {
							// add Xi, and Kronecker delta terms, c3
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   + f_n.transpose()*n_d);
						} else if (q_k(k1,nTgt,4) == q) {
							// add Xi, and Kronecker delta terms, c4
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   - f_n.transpose()*n_e);
						} // end if k1,q (targets)

						// contributions from sources
						if (zetaSrc.data() == zetaTgt.data()) {
							if (q_k(k2,nSrc,ll) == q) {
								// add Kronecker delta term, segment start
								dX.block<1,3>(k2,3*q) += -a*(f_r0+f_r1).transpose();
							} else if (q_k(k2,nSrc,llp1) == q) {
								dX.block<1,3>(k2,3*q) += a*(f_r0-f_r2).transpose();
							} // end if k2,q (sources)
						} // end if Src == Tgt

					} else {
						continue;
					} // end if k,q
				} // for l
			} // for q
		} // for k2
	} // for k1
	return;
}

void AIC(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  double* dX_) {
	/**@brief Calculate AIC matrix.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @param dX K_tgt x K_src matrix output.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,mTgt*nTgt,mSrc*nSrc);
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1, c2, c3, c4, cp; // panel corner points, collocation point
	Vector3d d, e, n; // panel diagonal and normal vectors
	Vector3d r0, r1, r2; // Biot-Savart kernel vectors

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// diagonals
		d = c3 - c1;
		e = c2 - c4;
		// calc n
		n = (d.cross(e)).normalized();
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			for (unsigned int l = 1; l < 5; l++) {
				// roll around segment index
				if (l < 4) {
					ll = l;
					llp1 = l+1;
				} else if (l == 4) {
					ll = l;
					llp1= 1;
				}

				// calc r0
				r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
					-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r1
				r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r2
				r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);
				// AIC entry
				dX(k1,k2)+=1.0/(4.0*M_PI)*fGeom(r0.data(),
										    r1.data(),
										    r2.data(),
										    n.data());
			}
		}
	}
	return;
}
