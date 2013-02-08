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
#include <boost/multi_array.hpp>

typedef boost::multi_array<double, 2> BoostArray2D;
typedef boost::multi_array<double, 3> BoostArray3D;

typedef Eigen::Map<Eigen::VectorXd> EigenMapVectXd;

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


void C_BiotSegment_ImageYZ(const double* xP,\
		  const double* x1,\
		  const double* x2,\
		  double& Gamma,\
		  double* Uind) {
	/** @brief calculate effect due to mirror of segment in YZ-plane.
	 *  @details Uind is overwritten.
	 */
	// create target image location (reverse x ordinate)
	double xImage[3] = {0.0,0.0,0.0};
	xImage[0] = -xP[0];
	xImage[1] = xP[1];
	xImage[2] = xP[2];

	//create result
	double VelImage[3] = {0.0,0.0,0.0};

	// segment 2 (point 2 -> point 3)
	C_BiotSegment(xImage, \
				  x1,x2, \
				  Gamma, VelImage);

	//reverse x velocity
	VelImage[0] = -VelImage[0];

	//copy to output
	CopyTriad(Uind,VelImage);
}


void BiotSavartSurf(const double* Zeta_Vec, const double* Gamma_Vec, \
					const double* TargetTriad,\
					const unsigned int Mstart, const unsigned int Nstart,
					const unsigned int Mend, const unsigned int Nend,
					const unsigned int Mfull, const unsigned int Nfull, \
					const bool ImageMethod, \
					double* Uout) {
	/**@brief Velocity induced by grid.
	 * @param Zeta_Vec Aero grid - size (M+1)*(N+1)*3.
	 * @param Gamma_Vec Vortex ring circulation strengths - size M*N.
	 * @param TargetTriad Point at which velocity is required - size 3.
	 */

	//create some useful pointers
	double (*Zeta)[Nfull+1][3] = (double (*)[Nfull+1][3]) Zeta_Vec;
	double (*Gamma)[Nfull] = (double (*)[Nfull]) Gamma_Vec;

	//slow method (every segment of every ring)
	//TODO: differencing for individual segment strengths, required for KJ methods.

	//set Uout to zero
	Uout[0] = 0.0;
	Uout[1] = 0.0;
	Uout[2] = 0.0;

	//define temp variables
	double Temp1[3] = {0.0,0.0,0.0};
	double Temp2[3] = {0.0,0.0,0.0};

	//TODO: to make this parallel use locally defined variables NOT pointers
	//to global data
	for (unsigned int i = Mstart; i <= Mend; i++) {
		for (unsigned int j = Nstart; j <= Nend; j++) {
			// loop through each panel
			// reset running total
			Temp2[0] = 0.0;
			Temp2[1] = 0.0;
			Temp2[2] = 0.0;

			// get effect of all segments
			// segment 1 (point 1 -> point 2)
			C_BiotSegment(TargetTriad,Zeta[i][j],Zeta[i][j+1],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 2 (point 2 -> point 3)
			C_BiotSegment(TargetTriad,Zeta[i][j+1],Zeta[i+1][j+1],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 3 (point 3 -> point 4)
			C_BiotSegment(TargetTriad,Zeta[i+1][j+1],Zeta[i+1][j],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 3 (point 4 -> point 1)
			C_BiotSegment(TargetTriad,Zeta[i+1][j],Zeta[i][j],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			if (ImageMethod == 1) {
				// get effect of all segments
				// segment 1 (point 1 -> point 2) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i][j],Zeta[i][j+1],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 2 (point 2 -> point 3) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i][j+1],Zeta[i+1][j+1],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 3 (point 3 -> point 4) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i+1][j+1],Zeta[i+1][j],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 3 (point 4 -> point 1) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i+1][j],Zeta[i][j],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);
			}

			//Add to Uout
			AddTriad(Uout,Temp2,Uout);
		} // END for j
	} //END for i
}


void KatzForces(const double* Zeta_Vec, const double* Gamma_Vec,\
				const double* ZetaStar_Vec, const double* GammaStar_Vec,\
				const double* ZetaDot_Vec, \
				const double* Uext_Vec, \
				VMopts VMOPTS,\
				const double* Gamma_tm1_Vec,\
				const double* Downwash,\
				const double* LiftVector_Vec,\
				double* Forces_Vec) {
	/** @brief Calculate panel forces at collocation points.
	 * @param Zeta_Vec Grid points.
	 * @param Gamma_Vec Bound vortex ring gammas.
	 * @param ZetaStar_Vec Wake grid points.
	 * @param GammaStar_Vec Wake vortex ring gammas.
	 * @param VMOPTS Simulation options.
	 * @param Gamma_tm1_Vec Bound vortex ring gammas from previous timestep.
	 * @param ForcesCol Unsteady Forces at collocation points.
	 */

	// cast vectors into useful pointers
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Gamma)[VMOPTS.N+1] = (double (*)[VMOPTS.N+1]) Gamma_Vec;
	double (*Gamma_tm1)[VMOPTS.N+1] = (double (*)[VMOPTS.N+1]) Gamma_tm1_Vec;
	double (*LiftVector)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) LiftVector_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	// declare local, automatically-managed dynamic memory using boost ...
	BoostArray2D alpha_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray3D Tau_c_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D Tau_s_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray2D Del_c_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D Del_s_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D Area_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D dGamma_dt_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray3D Uwake_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D ForcesCol_(boost::extents[VMOPTS.M][VMOPTS.N][3]);

	// ... and get useful pointers to that data
	double (*alpha)[VMOPTS.N] = (double (*)[VMOPTS.N]) alpha_.data();
	double (*Tau_c)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Tau_c_.data();
	double (*Tau_s)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Tau_s_.data();
	double (*Del_c)[VMOPTS.N] = (double (*)[VMOPTS.N]) Del_c_.data();
	double (*Del_s)[VMOPTS.N] = (double (*)[VMOPTS.N]) Del_s_.data();
	double (*Area)[VMOPTS.N] = (double (*)[VMOPTS.N]) Area_.data();
	double (*dGamma_dt)[VMOPTS.N] = (double (*)[VMOPTS.N]) dGamma_dt_.data();
	double (*Uwake)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Uwake_.data();
	double (*ForcesCol)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3])ForcesCol_.data();

	//Temporary variables
	double ZetaDotCol[3] = {0.0,0.0,0.0};
	double UextCol[3] = {0.0,0.0,0.0};
	double NormalWashCol[3] = {0.0,0.0,0.0};
	double Normal[3] = {0.0,0.0,0.0};
	double Collocation[3] = {0.0,0.0,0.0};
	double LiftVel[3] = {0.0,0.0,0.0};
	double TempC = 0.0;
	double TempS = 0.0;
	double TempGamma_i = 0.0;
	double TempGamma_j = 0.0;
	double DeltaP = 0.0;
	double LiftLocal = 0.0;
	double DragLocal1 = 0.0;
	double DragLocal2 = 0.0;
	double LiftTemp[3] = {0.0,0.0,0.0};
	double DragTemp[3] = {0.0,0.0,0.0};
	double NormNormalWashCol[3] = {0.0,0.0,0.0};

	//loop through collocation points i,j
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//calculate velocity of collocation
			BilinearMapTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
							 ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
							 ZetaDotCol);

			//External fluid velocities at collocation
			BilinearMapTriad(Uext[i][j], Uext[i][j+1], \
							 Uext[i+1][j+1], Uext[i+1][j], \
							 UextCol);

			// set incident velocity
			NormalWashCol[0] = -ZetaDotCol[0] + UextCol[0];
			NormalWashCol[1] = -ZetaDotCol[1] + UextCol[1];
			NormalWashCol[2] = -ZetaDotCol[2] + UextCol[2];


			//calculate panel collocation
			BilinearMapTriad(Zeta[i][j], Zeta[i][j+1], \
							 Zeta[i+1][j+1], Zeta[i+1][j], \
							 Collocation);


			//calculate panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
						Zeta[i+1][j+1], Zeta[i+1][j],
								    Normal);


			// calculate panel tangential vectors
			PanelTau_c(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j],
					   Tau_c[i][j]);

			PanelTau_s(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j],
					   Tau_s[i][j]);


			// calculate local AoA
			alpha[i][j] = atan2(DotTriad(NormalWashCol,Normal),\
							    DotTriad(NormalWashCol,Tau_c[i][j]));


			// calculate panel \Delta chord and \Delta span
			Del_c[i][j] = PanelDeltaC(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j]);

			Del_s[i][j] = PanelDeltaS(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel areas
			Area[i][j] = PanelArea(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel d(Gamma)/dt
			//dGamma_dt[i][j] = Gamma[i][j] - Gamma_tm1[i][j];
			dGamma_dt[i][j] = 0.0;


			// calculate wake induced velocity
			BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, Collocation, \
					0, 0, 0, VMOPTS.N, \
					1,VMOPTS.N, VMOPTS.ImageMethod,\
					Uwake[i][j]);


			// calculate local lift from Simpson (2013) plus external velocities
			// (-ZetaDot+Uwake+Uext)
			AddTriad(NormalWashCol,Uwake[i][j],LiftVel);
			AddTriad(LiftVel,Uext[i][j],LiftVel);

			// dot product with tangential vectors
			TempC = DotTriad(LiftVel,Tau_c[i][j]);
			TempS = DotTriad(LiftVel,Tau_s[i][j]);

			//calculate spatial Delta Gamma_i
			if ( i==0 ) {
				TempGamma_i = Gamma[i][j];
			} else {
				TempGamma_i = Gamma[i][j] - Gamma[i-1][j];
			}

			//calculate spatial Delta Gamma_j
			// firstly for if there is an image plane at wing root
			if (VMOPTS.ImageMethod == 1) {
				if (j == 0) {
					TempGamma_j = 0.0;
				} else if (j > 0) {
					TempGamma_j = Gamma[i][j] - Gamma[i][j-1];
					//if the whole wing is modelled then
				} else if (VMOPTS.ImageMethod == 0) {
				//TODO: this.
				}
			}

			// pressure jump
			DeltaP = (TempC * TempGamma_i / Del_c[i][j] + \
					TempS * TempGamma_j / Del_s[i][j] + \
					dGamma_dt[i][j]);


			// calculate local lift
			LiftLocal = DeltaP * Area[i][j] * cos(alpha[i][j]);

			// calculate local drag
			printf("\tDownwash = %f\n",Downwash[i*VMOPTS.N + j]);
			DragLocal1 = -Downwash[i*VMOPTS.N + j]*TempGamma_i*Del_s[i][j];
			DragLocal2 = dGamma_dt[i][j]*Area[i][j] * sin(alpha[i][j]);
			printf("\tdrag1:%f\tdrag2:%f\n",DragLocal1,DragLocal2);

			// total force
			//lift vector
			MulTriad(LiftVector[i][j],LiftLocal,LiftTemp);
			//PrintTriad(LiftVector[i][j]);
			//PrintTriad(LiftTemp);

			//drag vector
			NormaliseTriad(NormalWashCol,NormNormalWashCol);
			if (NormTriad(NormalWashCol) == 0.0) {
				std::cerr << "VLM: Warning incident velocity Zero!" \
										  << std::endl;
			}
			MulTriad(NormNormalWashCol,DragLocal1+DragLocal2,DragTemp);
			//PrintTriad(NormNormalWashCol);
			//PrintTriad(DragTemp);

			//combine and save
			AddTriad(LiftTemp,DragTemp,LiftTemp);

//			printf("\tlift:%f\tdrag:%f\tVector:",LiftLocal,DragLocal);
//			PrintTriad(LiftTemp);
//			printf("\n");

			CopyTriad(ForcesCol[i][j],LiftTemp);

			CopyTriad(Forces[i][j],ForcesCol[i][j]);
		}
	}

}


void cpp_solver_vlm(const double* Zeta_Vec, const double* ZetaDot_Vec, \
				    const double* Uext_Vec, double* ZetaStar_Vec, \
				    VMopts VMOPTS, \
				    double* Forces_Vec, \
				    double* Gamma_Vec, double* GammaStar_Vec) {
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

	//Declare dynamic memory for local variables using boost multi_array...
	/**@warning Once this function ends this data does not exist anymore */
	BoostArray3D ZetaCol_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);
	BoostArray3D ZetaDotCol_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);
	BoostArray3D UextCol_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);
	BoostArray3D Normal_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);
	BoostArray3D NormalWashCol_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);
	BoostArray3D LocalLift_(boost::extents[VMOPTS.M+1][VMOPTS.N+1][3]);


	// ... and useful array pointers to that memory.
	double (*ZetaCol)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaCol_.data();
	double (*ZetaDotCol)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDotCol_.data();
	double (*UextCol)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) UextCol_.data();
	double (*Normal)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Normal_.data();
	double (*NormalWashCol)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) NormalWashCol_.data();
	double (*LocalLift)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) LocalLift_.data();


	// Declare Eigen types for linear algebra
	Eigen::VectorXd RHS(VMOPTS.M*VMOPTS.N);
	EigenMapVectXd Gamma(Gamma_Vec,VMOPTS.M*VMOPTS.N);
	EigenMapVectXd GammaStar(GammaStar_Vec,VMOPTS.N);
	Eigen::MatrixXd AIC(VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	Eigen::MatrixXd BIC(VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	Eigen::VectorXd Downwash(VMOPTS.M*VMOPTS.N);

	//temporary variables
	double GammaOne = 1.0;
	double Temp1[3] = {0.0,0.0,0.0};
	double Temp2[3] = {0.0,0.0,0.0};
	double Temp3[3] = {0.0,0.0,0.0};
	Eigen::Matrix3d Proj;
	Eigen::Matrix3d Eye = Eigen::Matrix3d::Identity();
	Eigen::Vector3d Uincident(0.0,0.0,0.0);
	Eigen::Vector3d LocalLiftEig(0.0,0.0,0.0);
	Eigen::Vector3d NormalEig(0.0,0.0,0.0);

	//use Eigen to set GammaStar_Vec to 1.0 for all j
	GammaStar.setOnes();


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
							 &UextCol[i][j][0]);

			//panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
					    Zeta[i+1][j+1], Zeta[i+1][j],
					    Normal[i][j]);

			//Normal wash at Collocations points
			SubTriad(UextCol[i][j],ZetaCol[i][j],NormalWashCol[i][j]);

			// Fill RHS Vector
			RHS(i*VMOPTS.N + j) = - DotTriad(NormalWashCol[i][j],
											 Normal[i][j]);


			// calc local lift vector for BIC calculation (as in Simpson, 2013)
			// calc orthogonal project operator (Eigen types)
			Uincident(0) = -ZetaDotCol[i][j][0] + Uext[i][j][0];
			Uincident(1) = -ZetaDotCol[i][j][1] + Uext[i][j][1];
			Uincident(2) = -ZetaDotCol[i][j][2] + Uext[i][j][2];

			//check incident velocity and define projection operator
			if (Uincident.norm() == 0.0) {
				std::cerr << "VLM: Warning incident velocity Zero!" \
						  << std::endl;
				Proj = Eye;
			} else {
				Proj = Eye - ( Uincident.normalized() ) * \
							( Uincident.normalized().transpose() );
			} // END if, else

			// get local lift vector
			NormalEig(0) = Normal[i][j][0];
			NormalEig(1) = Normal[i][j][1];
			NormalEig(2) = Normal[i][j][2];
			LocalLiftEig = Proj*NormalEig;
			LocalLiftEig.normalize();
			// populate local lift LocalLift
			LocalLift[i][j][0] = LocalLiftEig[0];
			LocalLift[i][j][1] = LocalLiftEig[1];
			LocalLift[i][j][2] = LocalLiftEig[2];

			// Loop through each vortex ring (panel) ii, jj
			// TODO: if this is to be parallelised write to variables local to
			// the scope NOT the global pointers!
			for (unsigned int ii = 0; ii < VMOPTS.M; ii++) {
				for (unsigned int jj = 0; jj < VMOPTS.N; jj++) {
					//Induced vel of ring ii,jj to panel i,j

					// reset running total to zero
					Temp3[0] = 0.0;
					Temp3[1] = 0.0;
					Temp3[2] = 0.0;


					// chordwise orientated vorticity only (for BIC)

					// segment 2 (point 2 -> point 3)
					C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj+1],Zeta[ii+1][jj+1],
								  GammaOne, Temp1);

					// segment 4 (point 4 -> point 1)
					C_BiotSegment(ZetaCol[i][j],Zeta[ii+1][jj],Zeta[ii][jj],
								  GammaOne, Temp2);

					// sum of segments 2 and 4
					AddTriad(Temp1,Temp2,Temp3); //Temp3 overwritten here

					//image method
					if (VMOPTS.ImageMethod == 1) {
						// segment 2 (point 2 -> point 3) image
						C_BiotSegment_ImageYZ(ZetaCol[i][j],\
											  Zeta[ii][jj+1],Zeta[ii+1][jj+1], \
											  GammaOne, Temp1);

						AddTriad(Temp3,Temp1,Temp3);


						// segment 4 (point 4 -> point 1) image
						C_BiotSegment_ImageYZ(ZetaCol[i][j],\
											  Zeta[ii+1][jj],Zeta[ii][jj],\
											  GammaOne, Temp1);

						AddTriad(Temp3,Temp1,Temp3);
					}


					// if we're at the trailing edge then segment3 must be added
					// also the wake effect is added here
					if (ii == VMOPTS.M-1) {
						// add TE segment
						C_BiotSegment(ZetaCol[i][j],\
									  Zeta[ii+1][jj+1],\
									  Zeta[ii+1][jj],\
									  GammaOne,\
									  Temp1);

						AddTriad(Temp3,Temp1,Temp3);

						//image method
						if (VMOPTS.ImageMethod == 1) {
							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
												  Zeta[ii+1][jj+1],\
												  Zeta[ii+1][jj],\
												  GammaOne, Temp1);
							AddTriad(Temp3,Temp1,Temp3);
						} // END if Image Method

						// add wake effect
						BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, \
									   ZetaCol[i][j],\
									   0, jj, \
									   0, jj, \
									   1, VMOPTS.N, VMOPTS.ImageMethod,\
									   Temp1);
						AddTriad(Temp3,Temp1,Temp3);
					} // END if VMOPTS.M-1


					// element of BIC matrix
					BIC(i*VMOPTS.N + j,ii*VMOPTS.N + jj) = \
							DotTriad(Temp3,LocalLift[i][j]);


					//calculate induced velocity for AIC calculation
					// both spanwise-orientated segments must be added UNLESS
					// we are at the trailing edge. At the TE we add the only
					// remaining bound segment.

					if (ii < VMOPTS.M-1) {
						//add segments 1 and 3
						C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj],Zeta[ii][jj+1],
													GammaOne, Temp1);
						C_BiotSegment(ZetaCol[i][j],Zeta[ii+1][jj+1],Zeta[ii+1][jj],
													GammaOne, Temp2);
						AddTriad(Temp3,Temp1,Temp3);
						AddTriad(Temp3,Temp2,Temp3);

						if (VMOPTS.ImageMethod == 1) {
							//add segments 1 and 3 image
							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
												  Zeta[ii][jj],\
												  Zeta[ii][jj+1],\
												  GammaOne,\
												  Temp1);

							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
												  Zeta[ii+1][jj+1],\
												  Zeta[ii+1][jj],\
												  GammaOne,\
												  Temp2);
							AddTriad(Temp3,Temp1,Temp3);
							AddTriad(Temp3,Temp2,Temp3);

						} // END if Image Method

					} else if (ii == VMOPTS.M-1) {
						//add segments 1 only
						C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj],Zeta[ii][jj+1],
													GammaOne, Temp1);

						AddTriad(Temp3,Temp1,Temp3);
						if (VMOPTS.ImageMethod == 1) {
							//add segment 1 image
							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
											      Zeta[ii][jj],\
											      Zeta[ii][jj+1],
												  GammaOne, Temp1);

							AddTriad(Temp3,Temp1,Temp3);

						} // END if Image Method

					} //END if, else if VMOPTS.M-1

					// element of AIC matrix
					AIC(i*VMOPTS.N + j,ii*VMOPTS.N + jj) = \
												DotTriad(Temp3,Normal[i][j]);

				} //END for jj
			} //END for ii
		} //END for j
	} //END for i


	//solve for gamma
	Gamma = AIC.colPivHouseholderQr().solve(RHS);

	//set GammaStar to TE velocities
	GammaStar = Gamma.tail(VMOPTS.N);

	//calculate downwash
	Downwash = BIC * (Gamma);

	double* Downwash_ptr = Downwash.data();

	// Calculate forces
	BoostArray2D FooGamma_tm1_(boost::extents[VMOPTS.M][VMOPTS.N]);
	double* FooGamma_tm1_vec = FooGamma_tm1_.data();
	KatzForces(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,\
			   ZetaDot_Vec, Uext_Vec, VMOPTS, FooGamma_tm1_vec, Downwash_ptr,\
			   LocalLift_.data(), Forces_Vec);
}






