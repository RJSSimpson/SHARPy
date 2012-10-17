function zeph_shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeph_shape.m  Henrik Hesse 07/01/2011
% 
% 3D rendering of Zephyr-type of aircraft with BeamLength1 for wing span,
% BeamLength2 for fuselage length and BeamLength3 for tail wing. fw and fl
% are mesh density parameters for wings and fuselage respectively. 
%
% For now only implemented for static deformation. Look more into dynamic
% animations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Geometry
BeamLength1 = 20;
BeamLength2 = 16;
BeamLength3 = BeamLength1/4;
fw=2;
fl=1;

NumEl =   fl*(BeamLength2+BeamLength3)+fw*2*(BeamLength1+BeamLength3);
Nodes =   2;
shape=load('fine/ZEPH_SOL112_def10.txt');

foilRef=25;
foilChord=4;
foilT=0.3;

fuseD=0.15;

% Init
fusek=1;
wingk=1;
tailk=1;
boomk=1;
lift1=[];

for i=1:2:NumEl*Nodes
    % Reference lines
%    plot3(shape(i:i+1,3),shape(i:i+1,4),shape(i:i+1,5),'-')
%    hold on
        
    % Rotation matrix from CRV
    Psi=shape(i,6:8)';
    psi=norm(Psi);
    if psi>1E-6
        R=eye(3) + sin(psi)./psi.*Skew(Psi) + (1-cos(psi))./(psi^2).*(Skew(Psi))^2;
    else
        R=eye(3);
    end
    
    if (i <= fl*BeamLength2*Nodes)
    % Fuselage
        circle=[];
        circle(:,2)=(fuseD.*sin(0:pi/8:2*pi))';
        circle(:,3)=(fuseD.*cos(0:pi/8:2*pi))';
        for j=1:size(circle,1)
            circle(j,:)=(R*circle(j,:)')';
            circle(j,:)=circle(j,:)+shape(i,3:5);
        end
        plot3(circle(:,1),circle(:,2),circle(:,3),'k')
        hold on
        
        fuseX(fusek,:)=circle(:,1)';
        fuseY(fusek,:)=circle(:,2)';
        fuseZ(fusek,:)=circle(:,3)';
        fusek=fusek+1;
                
        % Include last node
        if (i == (fl*BeamLength2*Nodes-1))
            circle=[];
            circle(:,2)=(fuseD.*sin(0:pi/8:2*pi))';
            circle(:,3)=(fuseD.*cos(0:pi/8:2*pi))';
            for j=1:size(circle,1)
                circle(j,:)=(R*circle(j,:)')';
                circle(j,:)=circle(j,:)+shape(i+1,3:5);
            end
            plot3(circle(:,1),circle(:,2),circle(:,3),'k')
        
            fuseX(fusek,:)=circle(:,1)';
            fuseY(fusek,:)=circle(:,2)';
            fuseZ(fusek,:)=circle(:,3)';
            fusek=fusek+1;      
        end
        
    elseif (i <= (fl*BeamLength2+fw*2*BeamLength1)*Nodes)
    % Main wing    
        foil=naca(foilRef,foilChord,foilT);
        for j=1:size(foil,1)
            foil(j,:)=(R*foil(j,:)')';
            foil(j,:)=foil(j,:)+shape(i,3:5);
        end
        plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
        wingX(wingk,:)=foil(:,1)';
        wingY(wingk,:)=foil(:,2)';
        wingZ(wingk,:)=foil(:,3)';
        wingk=wingk+1;
        
        % Include last node
        if (i == (fl*BeamLength2+fw*2*BeamLength1)*Nodes-1)
            foil=naca(foilRef,foilChord,foilT);
            for j=1:size(foil,1)
                foil(j,:)=(R*foil(j,:)')';
                foil(j,:)=foil(j,:)+shape(i+1,3:5);
            end
            plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
            wingX(wingk,:)=foil(:,1)';
            wingY(wingk,:)=foil(:,2)';
            wingZ(wingk,:)=foil(:,3)';
            wingk=wingk+1;            
        end
        
        % Lift distribution
        X(1,:)=0.5.*(shape(i,3:5)+shape(i+1,3:5));
        if (i==fl*BeamLength2*Nodes+1)
            lift1(1,:)=shape(fl*BeamLength2*Nodes+1,3:5);
            Z=(R*[0;0;1.5])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift1(wingk,:)=X(2,:);
        elseif (i == (fl*BeamLength2+fw*2*BeamLength1)*Nodes-1)
            lift1(wingk,:)=shape((fl*BeamLength2+fw*2*BeamLength1)*Nodes,3:5);
            Z=(R*[0;0;1.5])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift1(wingk-1,:)=X(2,:);
        else
            Z=(R*[0;0;2.5])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift1(wingk,:)=X(2,:);    
        end


        
    elseif (i <= (fl*BeamLength2+fw*2*(BeamLength1+BeamLength3))*Nodes)
    % Tail wing    
        foil=naca(foilRef,foilChord/2,foilT);
        for j=1:size(foil,1)
            foil(j,:)=(R*foil(j,:)')';
            foil(j,:)=foil(j,:)+shape(i,3:5);
        end
        plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
        tailX(tailk,:)=foil(:,1)';
        tailY(tailk,:)=foil(:,2)';
        tailZ(tailk,:)=foil(:,3)';
        tailk=tailk+1;   
        
        % Include last node
        if (i == (fl*BeamLength2+fw*2*(BeamLength1+BeamLength3))*Nodes-1)
            foil=naca(foilRef,foilChord/2,foilT);
            for j=1:size(foil,1)
                foil(j,:)=(R*foil(j,:)')';
                foil(j,:)=foil(j,:)+shape(i+1,3:5);
            end
            plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
            tailX(tailk,:)=foil(:,1)';
            tailY(tailk,:)=foil(:,2)';
            tailZ(tailk,:)=foil(:,3)';
            tailk=tailk+1;  
        end
                
        % Lift distribution
        X(1,:)=0.5.*(shape(i,3:5)+shape(i+1,3:5));
        if (i==(fl*BeamLength2+fw*2*BeamLength1)*Nodes+1)
            lift2(1,:)=shape((fl*BeamLength2+fw*2*BeamLength1)*Nodes+1,3:5);
            Z=(R*[0;0;1])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift2(tailk,:)=X(2,:);
        elseif (i == (fl*BeamLength2+fw*2*(BeamLength1+BeamLength3))*Nodes-1)
            lift2(tailk,:)=shape((fl*BeamLength2+fw*2*(BeamLength1+BeamLength3))*Nodes,3:5);
            Z=(R*[0;0;1])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift2(tailk-1,:)=X(2,:);
        else
            Z=(R*[0;0;2])';
            X(2,:)=X(1,:)+Z;
            plot3(X(:,1),X(:,2),X(:,3),'r')
            lift2(tailk,:)=X(2,:);    
        end

    elseif (i <= (fl*(BeamLength2+BeamLength3)+fw*2*(BeamLength1+BeamLength3)-1)*Nodes)
    % Tail boom    
        foil=naca(foilRef,foilChord/2,foilT);
        for j=1:size(foil,1)
            foil(j,:)=(R*foil(j,:)')';
            foil(j,:)=foil(j,:)+shape(i,3:5);
        end
        plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
        boomX(boomk,:)=foil(:,1)';
        boomY(boomk,:)=foil(:,2)';
        boomZ(boomk,:)=foil(:,3)';
        boomk=boomk+1;       
        
        % Include last node
        if (i == (fl*(BeamLength2+BeamLength3)+fw*2*(BeamLength1+BeamLength3))*Nodes-1)
            foil=naca(foilRef,foilChord/2,foilT);
            for j=1:size(foil,1)
                foil(j,:)=(R*foil(j,:)')';
                foil(j,:)=foil(j,:)+shape(i+1,3:5);
            end
            plot3(foil(:,1),foil(:,2),foil(:,3),'k')
        
            boomX(boomk,:)=foil(:,1)';
            boomY(boomk,:)=foil(:,2)';
            boomZ(boomk,:)=foil(:,3)';
            boomk=boomk+1;     
        end
        
    end
end

surfl(fuseX,fuseY,fuseZ)
surfl(wingX,wingY,wingZ)
surfl(tailX,tailY,tailZ)
surfl(boomX,boomY,boomZ)

shading interp
colormap(gray);
lighting phong
axis equal
grid on
axis([-22 22 -20 10 -5 15])

% Interpolate lift and plot
lift1int(:,1)=lift1(1,1):.1:lift1(end,1);
lift1int(:,3)=spline(lift1(:,1),lift1(:,3),lift1int(:,1));
plot3(lift1int(:,1),lift1int(:,2),lift1int(:,3),'r')

lift2int(:,1)=lift2(1,1):.1:lift2(end,1);
lift2int(:,2)=-3*BeamLength2/4;
lift2int(:,3)=spline(lift2(:,1),lift2(:,3),lift2int(:,1));
plot3(lift2int(:,1),lift2int(:,2),lift2int(:,3),'r')


function Omega = Skew(a)
Omega = [  0    -a(3)  a(2)
          a(3)    0   -a(1) 
         -a(2)   a(1)   0  ];