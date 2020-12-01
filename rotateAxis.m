function [x_t] = rotateAxis(x,Theta,Phi)
%Rotates x from Matrad coordinate system according to couch angle Theta and
%gantry angle Phi in the Topas coordinate system

%Swap axis (Matrad/Topas conversion)
x = x([2 1 3],:);

%Rotation matrices
Theta=Theta*pi/180;
Phi=Phi*pi/180;
rotTheta=[cos(Theta) 0 sin(Theta); 0 1 0; -sin(Theta) 0 cos(Theta)];
rotPhi=[cos(Phi) -sin(Phi) 0; sin(Phi) cos(Phi) 0; 0 0 1];

%Perform rotation
x_t=rotTheta*rotPhi*x;

%Swap axis back
x_t = x_t([2 1 3],:);
end

