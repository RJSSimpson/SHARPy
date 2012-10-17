%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_shape.m  Henrik Hesse 07/01/2011
% 
% Animate structural dynamics simulation obtained from XBeam. *shape.txt
% containes the elastic deformation wrt the global frame a. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

NumEl =    2;
Nodes =    2;
dt    =  0.1;
NumN=NumEl*(Nodes-1)+1;

shape=load('WING_SOL212_shape.txt');

for i=1:size(shape,1)/NumN,
    lr = (NumN+1)/2; 
    right        = shape((i-1)*NumN+1:(i-1)*NumN+lr,2:4);
    left(1,:)    = shape((i-1)*NumN+1,2:4);
    left(2:lr,:) = shape((i-1)*NumN+lr+1:i*NumN,2:4);

    figure(4)
    plot3(right(:,1),right(:,2),right(:,3),'-o')
%    plot(right(:,1),right(:,3),'-xr')
    hold on
    axis([-10 10 -10 10 -10 10])
    plot3(left(:,1),left(:,2),left(:,3),'-o')
%    plot(left(:,1),left(:,3),'-xr')
    axis equal
    grid
    hold off
    pause(0.01)
%}    
end

right=0; left=0;
t=dt:dt:size(shape,1)*dt/NumN;
for i=1:size(shape,1)/NumN,
    lr = (NumN+1)/2; 
    right(i) = shape(NumN*(i-1)+(NumN+1)/2,2);
    left (i) = shape(NumN*i,2);
end

length(t)
length(right)

figure(6)
hold on
plot(t,right)
plot(t,left,'r-')