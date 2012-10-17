%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_RigStruc.m  Henrik Hesse 07/01/2011
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
function res=check_CantStruc()
% clear;
% clc;

NumN = 11;
Time = 0:0.1:6;

shape=load('SIMO_SOL912_shape.txt');
rigid=load('SIMO_SOL912_rigid.txt');
%Coa=load('rotation.txt');

S0(1:6) = 0;
S0(7:10)= [ 1; 0; 0; 0];

[S] = runge(@eom,S0,Time,rigid(:,2:7));

for i=1:1:length(Time)-2,
    % Generate rotation matrix at each time step i
    q = S(i,7:10);
    R = Rotation(q); 
    
    % plot aircraft
    right= shape((i-1)*NumN+1:i*NumN,2:4);
   
    for j=1:size(right,1),
        right(j,:)=(R*right(j,:)'+S(i,1:3)')';
    end

    figure(33)
    hold on
%    plot3(right(:,1),right(:,2),right(:,3),'-b')
    plot(right(:,1),right(:,3),'-r','LineWidth',2)
%    plot(right(1,2),right(1,3),'.r')
    set(gcf,'Color',[1,1,1])
%     axis equal
%     axis([-10 10 -2 10])
    grid
    pause(0.05)
    hold on
%    F(i/1) = getframe(gcf);
end

% % plot velocities
% figure(1)      
% hold on
% plot(rigid(:,1),rigid(:,5),'k')
% figure(2)
% hold on
% plot(rigid(:,1),rigid(:,6),'k')
% figure(3)
% hold on
% plot(rigid(:,1),rigid(:,7),'k')
% %movie2avi(F,'torque','compression','None')


function dSdt = eom(S,V)
% Equations of motion in the form of:
%       dx/dt = v
%       dq/dt = -0.5*OMEGA*q
% where OMEGA is an extended skew-symmetric matrix consisting of the
% angular velocity components (see Shearer and Cesnik (2007))
dSdt(1:3)  = V(1:3);
dSdt(4:6)  = V(4:6);
dSdt(7:10) = (0.5).*QuadSkew(V(4:6))*S(7:10)';
dSdt = dSdt';

function [y] = runge_solver(F,y0,t,rigid)
% Solving ode's using the Runge-Kutta scheme
y = [];
h = t(end)-t(end-1);

% BC
y(1,:)=y0'; 

for n=1:length(t)-1,
    k1 = h*F(y(n,:),rigid(n,:));
    k2 = h*F((y(n,:)+0.5.*k1'),rigid(n,:));
    k3 = h*F((y(n,:)+0.5.*k2'),rigid(n,:));
    k4 = h*F((y(n,:)+k3'),rigid(n,:));
    y(n+1,:) = y(n,:)+(1/6).*(k1+2.*k2+2.*k3+k4)';
end

function [Y] = runge(F,Y0,t,rigid)
% Initialize 
Y(1,:) = Y0;
step = round((t(end)-t(end-1))*1000)/1000;

for n=1:length(t)-1,
    % Interpolate
    Ytmp(n,:)  =F(Y(n,:),rigid(n,:));
    Ytmp(n+1,:)=F(Y(n,:),rigid(n+1,:));
    
    tinterp=t(n):step/2:t(n+1);
    Yinterp=interp1(t(n:n+1),Ytmp(n:n+1,:),tinterp,'spline');

    a = Yinterp(1,:);
    b = Yinterp(2,:);
    d = Yinterp(3,:);
    
    Y(n+1,:) = Y(n,:) + 1/6*(a+4*b+d)*step;

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