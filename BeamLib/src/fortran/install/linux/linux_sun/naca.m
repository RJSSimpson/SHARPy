function foil=naca(ref,c,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% naca.m  Henrik Hesse 07/01/2011
% 
% interpolated NACA 2412 profile
% with position of reference line in percent, chord c and half thickness t
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foil= [         0         0         0
         0   -0.0125    0.2728
         0   -0.0250    0.3794
         0   -0.0500    0.5241
         0   -0.0750    0.6294
         0   -0.1000    0.7145
         0   -0.1500    0.8388
         0   -0.2000    0.9213
         0   -0.2500    0.9734
         0   -0.3000    1.0000
         0   -0.4000    0.9898
         0   -0.5000    0.9188
         0   -0.6000    0.8071
         0   -0.7000    0.6574
         0   -0.8000    0.4759
         0   -0.9000    0.2640
         0   -0.9500    0.1447
         0   -1.0000         0
         0   -0.9500   -0.0609
         0   -0.9000   -0.1041
         0   -0.8000   -0.1904
         0   -0.7000   -0.2716
         0   -0.6000   -0.3503
         0   -0.5000   -0.4239
         0   -0.4000   -0.4822
         0   -0.3000   -0.5228
         0   -0.2500   -0.5355
         0   -0.2000   -0.5368
         0   -0.1500   -0.5203
         0   -0.1000   -0.4759
         0   -0.0750   -0.4391
         0   -0.0500   -0.3820
         0   -0.0250   -0.2881
         0   -0.0125   -0.2094
         0         0         0];

foil(:,2)=c.*(foil(:,2)+ref/100); 
foil(:,3)=t.*foil(:,3);