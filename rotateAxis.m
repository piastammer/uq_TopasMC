function [x_t] = rotateAxis(x,Theta,Phi)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Theta=Theta*pi/180;
Phi=Phi*pi/180;
rotTheta=[cos(Theta) 0 sin(Theta); 0 1 0; -sin(Theta) 0 cos(Theta)];
rotPhi=[cos(Phi) -sin(Phi) 0; sin(Phi) cos(Phi) 0; 0 0 1];
x_t=rotTheta*rotPhi*x;
end

