%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 MAIN PROGRAM                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              MEX ROUTINES                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set mexing options
!./mexopts.sh
disp('Mexing Fortran Routines')
disp('=======================')

%Input, initialization, geometric properties of elements and nodes
mex -g bgeom_elemcurv2_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g bgeom_elemframe_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g bgeom_elemlength_wrap.F90    -lxbeam -llapack -lblas -L"."
mex -g bgeom_nodeframe_wrap.F90     -lxbeam -llapack -lblas -L"."
disp('Finished mexing initialisation routines...')

%Matrices for the FE assembly
mex -g cbeam3_cgyr_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_cvel_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_fext_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_fgyr_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_fstif_wrap.F90    -lxbeam -llapack -lblas -L"."
mex -g cbeam3_kgeom_wrap.F90    -lxbeam -llapack -lblas -L"."
mex -g cbeam3_kgyr_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_kmat_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_mass_wrap.F90     -lxbeam -llapack -lblas -L"."
mex -g cbeam3_mvel_wrap.F90     -lxbeam -llapack -lblas -L"."
disp('Finished mexing FE routines...')

%Master/slave transformations
mex -g cbeam3_projm2s_wrap.F90  -lxbeam -llapack -lblas -L"."
mex -g cbeam3_projs2m_wrap.F90  -lxbeam -llapack -lblas -L"."
disp('Finished mexing master/slave routines...')

%Rigid-body dynamics
mex -g xbeam_cgyr_wrap.F90  -lxbeam -llapack -lblas -L"."
mex -g xbeam_crr_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g xbeam_fext_wrap.F90  -lxbeam -llapack -lblas -L"."
mex -g xbeam_fgyr_wrap.F90  -lxbeam -llapack -lblas -L"."
mex -g xbeam_kgyr_wrap.F90  -lxbeam -llapack -lblas -L"."
mex -g xbeam_kmass_wrap.F90 -lxbeam -llapack -lblas -L"."
mex -g xbeam_mrr_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g xbeam_mrs_wrap.F90   -lxbeam -llapack -lblas -L"."
disp('Finished mexing rigid-body routines...')

%Matrices for perturbed simulation (IFASD2011)
mex -g perturb_css_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g perturb_csr_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g perturb_crs_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g perturb_crr_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g perturb_kss_wrap.F90   -lxbeam -llapack -lblas -L"."
mex -g perturb_krs_wrap.F90   -lxbeam -llapack -lblas -L"."
disp('Finished mexing perturbed routines.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
