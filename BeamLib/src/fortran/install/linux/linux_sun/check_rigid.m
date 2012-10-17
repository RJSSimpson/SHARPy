%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK_RIGID.M  Henrik Hesse 07/01/2011
% 
% Check rigid-body velocities (velocities of a frame) coming from 
% XBeam simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
sharp=load('CANT_SOL912_rigid.txt');

figure(10)
hold on
grid on
plot(sharp(:,1),sharp(:,2),'r-')
figure(10)                          
hold on                           
grid on
plot(sharp(:,1),sharp(:,3),'r--')
figure(10)                      
hold on
grid on
plot(sharp(:,1),sharp(:,4),'r-.')
box on

figure(11)
hold on
grid on
plot(sharp(:,1),sharp(:,5),'r-')
figure(11)                          
hold on                           
grid on
plot(sharp(:,1),sharp(:,6),'r--')
figure(11)                      
hold on
grid on
plot(sharp(:,1),sharp(:,7),'r-.')
box on
