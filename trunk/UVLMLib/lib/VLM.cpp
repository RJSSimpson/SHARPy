/**@brief      quasi-steady VLM solution with fixed-wake.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       30/01/2013
 * @pre        None
 * @warning    None
 */

#include <stdio.h>
#include <datatypesx.hpp>
#include <triads.hpp>
#include <Eigen/Dense>
#include <vorticity.hpp>

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

void cpp_solver_vlm(const double* Zeta_Vec, const double* ZetaDot_Vec, \
				    const double* Uext_Vec, VMopts VMOPTS, double* Forces_Vec) {
	/* Calculate and save panel normal vectors, location of collocation point,
	 * velocity at collocation point, and the External fluid velocity at
	 * collocation points (using bilinear mapping on grid) and relative
	 *  normal-wash at collocation points */

	/** Here the contiguous data provided by Python - 'const double* Zeta_Vec' -
	 * needs to be accessed as if it were a 3D C-array -'double Zeta[M][N][3]'.
	 * To do this WITHOUT copying any data we create an appropriate pointer
	 * to the 3D array
	 * - 'double (*Zeta)[VMOPTS_N+1][3]' - that is initialised by 'casting'
	 * the contiguous data pointer 'Zeta_Vec' to the type required for our
	 * 3D array - '(*Zeta)[VMOPTS.N+1][3]'.
	 */
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;


	//Declare memory for local variables
	double ZetaCol[VMOPTS.M+1][VMOPTS.N+1][3];
	double ZetaDotCol[VMOPTS.M+1][VMOPTS.N+1][3];
	double UextCol[VMOPTS.M+1][VMOPTS.N+1][3];
	double Normal[VMOPTS.M+1][VMOPTS.N+1][3];
	double NormalWashCol[VMOPTS.M+1][VMOPTS.N+1][3];

	// Declare Eigen types for linear algebra
	Eigen::VectorXd RHS(VMOPTS.M*VMOPTS.N);
	Eigen::VectorXd Gamma(VMOPTS.M*VMOPTS.N);
	Eigen::MatrixXd AIC(VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	Eigen::MatrixXd BIC(VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	Eigen::MatrixXd AIC_wake(VMOPTS.M*VMOPTS.N,VMOPTS.N);

	//temporary variables
	double GammaOne = 1.0;
	double Temp1[3] = {0.0,0.0,0.0};
	double Temp2[3] = {0.0,0.0,0.0};
	double Temp3[3] = {0.0,0.0,0.0};
	Eigen::Matrix3d Proj;


	//loop through collocation points i,j
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//collocation points
			BilinearMapTriad(Zeta[i][j], Zeta[i][j+1], \
							 Zeta[i+1][j+1], Zeta[i+1][j], \
							 ZetaCol[i][j]);

			//relative inertial velocity at collocation points
			BilinearMapTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
							 ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
							 ZetaDotCol[i][j]);

			//External fluid velocities at collocation
			BilinearMapTriad(Uext[i][j], Uext[i][j+1], \
							 Uext[i+1][j+1], Uext[i+1][j], \
							 UextCol[i][j]);

			//panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
					    Zeta[i+1][j+1], Zeta[i+1][j],
					    Normal[i][j]);

			//Normal wash at Collocations points
			SubTriad(UextCol[i][j],ZetaCol[i][j],NormalWashCol[i][j]);

			// Fill RHS Vector
			RHS(i*VMOPTS.N + j) = - DotTriad(NormalWashCol[i][j],Normal[i][j]);

			// calc local lift vector for BIC calculation (as in Simpson, 2013)
			// calc orthogonal project operator


			// Loop through each vortex ring (panel) ii, jj
			for (unsigned int ii = 0; ii < VMOPTS.M; ii++) {
				for (unsigned int jj = 0; jj < VMOPTS.N; jj++) {

					/* Induced vel of ring ii,jj to panel i,j */

					// chordwise orientated vorticity only (for BIC)
					// segment 2 (point 2 -> point 3)
					C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj+1],Zeta[ii+1][jj+1],
								  GammaOne, Temp1);

					// segment 4 (point 4 -> point 1)
					C_BiotSegment(ZetaCol[i][j],Zeta[ii+1][jj],Zeta[ii][jj],
								  GammaOne, Temp2);

					// sum of segments 2 and 4
					AddTriad(Temp1,Temp2,Temp3);

					// element of BIC matrix
					//BIC(i*VMOPTS.N + j,ii*VMOPTS.N + jj) = \
							DotTriad(,)


				}
			}
		}
	}



}


