function res=check_tipreldif()
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

NumN  =      13;
t0    =     0.0;
te    =      10;
h     =    0.05;

shape912=load('DANC_SOL912_shape.txt');
shape952=load('DANC_SOL952_shape.txt');

rigid912=0*load('DANC_SOL912_rigid.txt');
rigid952=0*load('DANC_SOL952_rigid.txt');

shaperel=shape912-shape952;
rigidrel=rigid912-rigid952;

S0(1:6) = 0;
S0(7:10)= [ 1; 0; 0; 0];

[Time,S] = runge_solver(@eom,S0,[t0 te],h,rigidrel(:,2:7));

for i=1:1:length(Time)-1,
    % Generate rotation matrix at each time step i
    q = S(i,7:10);
    R = Rotation(q);

    % plot aircraft
    lr = (NumN+1)/2;
    rightrel         = shaperel((i-1)*NumN+1:(i-1)*NumN+lr,2:4);
    leftrel(1:lr-1,:)= shaperel((i-1)*NumN+lr+1:i*NumN,2:4);
    leftrel(lr,:)    = shaperel((i-1)*NumN+1,2:4);

    right912         = shape912((i-1)*NumN+1:(i-1)*NumN+lr,2:4);
    left912(1:lr-1,:)= shape912((i-1)*NumN+lr+1:i*NumN,2:4);
    left912(lr,:)    = shape912((i-1)*NumN+1,2:4);
    
    for j=1:size(rightrel,1),
        rightrel(j,:)=(R*rightrel(j,:)'+S(i,1:3)')';
        leftrel(j,:) =(R*leftrel(j,:)' +S(i,1:3)')';

        right912(j,:)=(R*right912(j,:)'+S(i,1:3)')';
        left912(j,:) =(R*left912(j,:)' +S(i,1:3)')';
    end
    
    tip(i,1)=Time(i);
    tip(i,2)=rightrel(j,2)/(right912(j,2)+1e4)*1e4;
    tip(i,3)=leftrel(1,2) /(left912(1,2)+1e4)*1e4;
    tip(i,4)=rightrel(j,3)/(right912(j,3)+1e4)*1e4;
    tip(i,5)=leftrel(1,3) /(left912(1,3)+1e4)*1e4;
end

max(abs(tip(:,2:5)))

% plot velocities
figure(12)      
    hold on
    plot(tip(:,1),tip(:,2),'-')
    plot(tip(:,1),tip(:,3),'--')    
    plot(tip(:,1),tip(:,4),'r-')
    plot(tip(:,1),tip(:,5),'r--')  

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