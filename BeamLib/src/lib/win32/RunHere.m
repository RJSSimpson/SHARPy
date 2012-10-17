%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 MAIN PROGRAM                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              MEX ROUTINES                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input, initialization, geometric properties of elements and nodes
mex -g bgeom_elemcurv2_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g bgeom_elemframe_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g bgeom_elemlength_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g bgeom_nodeframe_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."

%Matrices for the FE assembly
mex -g cbeam3_cgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_cvel_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_fext_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_dqext_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_fgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_fstif_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_kgeom_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_kgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_kmat_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_mass_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_mvel_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."

%Master/slave transformations
mex -g cbeam3_projm2s_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_projs2m_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."

%Rigid-body dynamics
mex -g xbeam_cgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_crr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_fext_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_fgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_kgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_kmass_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_mrr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_mrs_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."

%Lump masses
mex -g cbeam3_rbcgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_rbcvel_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_rbfgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_rbkgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_rbmass_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g cbeam3_rbmvel_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbcgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbcrr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbfgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbkgyr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbkmass_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbmrr_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."
mex -g xbeam_rbmrs_wrap.f90 xbeam.lib lapack.lib blas.lib -L"."

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%