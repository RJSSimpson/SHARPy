function res=check_CantTipdef()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_tipdef.m  Henrik Hesse 07/01/2011
% 
% Plot tip position wrt global a frame of coupled rigid-body and structural 
% dynamics simulation obtained from XBeam. 
%
% shape containes the elastic deformation wrt the global frame a and 
% rigid contains rigid-body velocities of the reference frame. The rotation
% matrix can be calculated from angular velocities or might be stored in
% rotation.txt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

NumN  =       51;
t0    =      0.0;
te    =      2.0;
h     =     0.05;

shape=load('GOLD_SOL212_shape.txt');
rigid=0*load('HALE_SOL910_rigid.txt');

S0(1:6) = 0;
S0(7:10)= [ 1; 0; 0; 0];

[Time,S] = runge_solver(@eom,S0,[t0 te],h,rigid(:,2:7));

for i=1:length(Time)-2,
    % Generate rotation matrix at each time step i
    q = S(i,7:10);
    R = Rotation(q);

    % plot aircraft
    right= shape((i-1)*NumN+1:i*NumN,2:4);
   
    for j=1:size(right,1),
        right(j,:)=(R*right(j,:)'+S(i,1:3)')';
    end
    
    tip(i,1)=Time(i);
    tip(i,2)=right(j,2);
    tip(i,3)=right(j,2);
end
    % plot velocities
    figure(1)      
    hold on
    plot(tip(:,1),tip(:,3),'r')
    plot(tip(:,1),tip(:,3),'r--')      
 
    
function dSdt = eom(t,S,V)
% Equations of motion in the form of:
%       dx/dt = v
%       dq/dt = -0.5*OMEGA*q
% where OMEGA is an extended skew-symmetric matrix consisting of the
% angular velocity components (see Shearer and Cesnik (2007))
dSdt(1:3)  = V(1:3);
dSdt(4:6)  = 0;
dSdt(7:10) = (0.5).*QuadSkew(V(4:6))*S(7:10);
dSdt = dSdt';


function [t,y] = runge_solver(F,y0,time,h,rigid)
% Solving ode's using the Runge-Kutta scheme
t = time(1):h:time(2);
y = [];
for i=1:length(y0), y(1,i) = y0(i); end
for n=1:length(t)-1,
    for i=1:length(y0),
        k1 = h*F(t(n),y(n,:)',rigid(n,:));
        k2 = h*F((t(n)+0.5*h),(y(n,:)'+0.5.*k1),rigid(n,:));
        k3 = h*F((t(n)+0.5*h),(y(n,:)'+0.5.*k2),rigid(n,:));
        k4 = h*F(t(n+1),(y(n,:)'+k3),rigid(n,:));
        y(n+1,:) = y(n,:)+(1/6).*(k1+2.*k2+2.*k3+k4)';
    end
end

function Omega = QuadSkew(a)
% Used to solve ode of quaternions
% See Shearer and Cesnik (2007)
Omega = [  0    a(1)  a(2)  a(3)
         -a(1)   0   -a(3)  a(2)
         -a(2)  a(3)   0   -a(1)
         -a(3) -a(2)  a(1)   0  ];
     
     
function R = Rotation(q)
% Definition of rotation operator in Euler quaternions
% See Aircraft Control and Simulation by Stevens, Lewis
R =[];

R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2; 
R(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
R(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

R(1,2) = 2*(q(2)*q(3) + q(1)*q(4));
R(2,1) = 2*(q(2)*q(3) - q(1)*q(4));
   
R(1,3) = 2*(q(2)*q(4) - q(1)*q(3));    
R(3,1) = 2*(q(2)*q(4) + q(1)*q(3));

R(2,3) = 2*(q(3)*q(4) + q(1)*q(2));        
R(3,2) = 2*(q(3)*q(4) - q(1)*q(2));       