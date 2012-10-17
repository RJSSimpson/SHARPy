%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK_RIGID.M  Henrik Hesse 07/01/2011
% 
% Check rigid-body velocities (velocities of a frame) coming from 
% XBeam simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rigid912=load('DANC_SOL912_rigid.txt');
rigid952=load('DANC_SOL952_rigid.txt');

for ii=1:size(rigid912,1)
    for jj=2:7
        rigid_rel(ii,jj)=(rigid912(ii,jj)-rigid952(ii,jj));
        rigid_rel(ii,jj)=rigid_rel(ii,jj)/(rigid912(ii,jj)+1e4)*1e4;
    end
end
rigid_rel(:,1)=rigid912(:,1);

figure(11)
hold on
grid on
plot(rigid_rel(:,1),rigid_rel(:,5),'-')
figure(11)                          
hold on                           
grid on
plot(rigid_rel(:,1),rigid_rel(:,6),'r--')
figure(11)                      
hold on
grid on
plot(rigid_rel(:,1),rigid_rel(:,7),'k--')
box on

max(abs(rigid_rel(:,2:7)))