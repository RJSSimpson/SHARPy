%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_RigStruc2.m  Henrik Hesse 07/01/2011
% 
% Animate coupled rigid-body and structural dynamics simulation obtained
% from XBeam. 
%
% shape containes the elastic deformation wrt the global frame a and 
% rigid contains rigid-body velocities of the reference frame. The rotation
% matrix can be calculated from angular velocities or might be stored in
% rotation.txt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=check_RigStruc2()
clear;
clc;

NumN  =       7;
t0    =      0.0;
te    =     10.0;
h     =     0.2;

shape=load('WINF_SOL912_shape.txt');
rigid=load('WINF_SOL912_rigid.txt');
%Coa=load('rotation.txt');

S0(1:6) = 0;
S0(7:10)= [ 1; 0; 0; 0];

[Time,S] = runge_solver(@eom,S0,[t0 te],h,rigid(:,2:7))

for i=1:1:length(Time)-1,
    % Generate rotation matrix at each time step i
    q = S(i,7:10)
    R = Rotation(q);
%    R = Coa((3*(i-1)+1:3*(i-1)+3),2:4);    
    
    % plot aircraft
    lr1 = (NumN-1)/3+1; 
    ll1 = (NumN-1)/3; 
    lr2 = (NumN-1)/6; 

    right1          = shape((i-1)*NumN+1:(i-1)*NumN+lr1,2:4);
    left1(1,:)      = shape((i-1)*NumN+1,2:4);
    left1(2:lr1,:)  = shape((i-1)*NumN+1+lr1:(i-1)*NumN+lr1+ll1,2:4);

    right2(1,:)     = shape((i-1)*NumN+lr1,2:4);
    right2(2:lr2+1,:) = shape((i-1)*NumN+1+lr1+ll1:(i-1)*NumN+lr1+ll1+lr2,2:4);

    left2(1,:)     = shape((i-1)*NumN+lr1+ll1,2:4);
    left2(2:lr2+1,:) = shape((i-1)*NumN+1+lr1+ll1+lr2:(i-1)*NumN+lr1+ll1+lr2+lr2,2:4);

    for j=1:size(right1,1),
        right1(j,:)=(R*right1(j,:)'+S(i,1:3)')';
        left1(j,:) =(R*left1(j,:)' +S(i,1:3)')';
    end
    for j=1:size(right2,1),
        right2(j,:)=(R*right2(j,:)'+S(i,1:3)')';
        left2(j,:) =(R*left2(j,:)' +S(i,1:3)')';
    end

    figure(4)
    plot3(right1(:,1),right1(:,2),right1(:,3),'-b','LineWidth',2)
    hold on
    plot3(right2(:,1),right2(:,2),right2(:,3),'-b','LineWidth',2)
    plot3(left1(:,1),left1(:,2),left1(:,3),'-r','LineWidth',2)
    plot3(left2(:,1),left2(:,2),left2(:,3),'-r','LineWidth',2)
    axis equal
    axis([-18 18 -18 18 -18 18])
    grid
    hold off
    pause(0.01)
%    F(i/1) = getframe(gcf);
end

% plot velocities
figure(1)      
hold on
plot(rigid(:,1),rigid(:,5),'k')
figure(2)
hold on
plot(rigid(:,1),rigid(:,6),'k')
figure(3)
hold on
plot(rigid(:,1),rigid(:,7),'k')
%movie2avi(F,'torque','compression','None')


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