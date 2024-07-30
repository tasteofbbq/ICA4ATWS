function [ n,U,Pw ] = MakeNurbsCircle( O,X,Y,r,ths,the )
% ALGORITHM A7.1

% Create arbitrary NURBS circular arc 
% Input: O,X,Y,r,ths,the 
% Output: n,U,Pw 
if (the < ths) 
    the = 360.0 + the;
end
theta = the-ths;
if (theta <= 90.0) 
    narcs = 1; % get number of arcs
elseif (theta <= 180.0)
    narcs = 2;
elseif (theta <= 270.0)
    narcs = 3;
else
narcs = 4;
end
dtheta = theta/narcs;

n = 2*narcs; % n+1 control points 
w1 = cos(dtheta/2.0*pi/180); % dtheta/2 is base angle
PO = O + r*cos(ths*pi/180)*X + r*sin(ths*pi/180)*Y;
TO = -sin(ths*pi/180)*X + cos(ths*pi/180)*Y; % Initialize start values
Pw(1,:) = PO;
index = 0; angle = ths;
for i=1:narcs % create narcs segments
    angle = angle + dtheta;
    P2 = O + r*cos(angle*pi/180)*X + r*sin(angle*pi/180)*Y;
    Pw(index+2+1,:) = P2;
    T2 = -sin(angle*pi/180)*X + cos(angle*pi/180)*Y;
    P1 = Intersect3DLines(PO,TO,P2,T2);
    Pw(index+2,:) = w1*P1;
     index = index + 2;
    if (i < narcs) 
        PO = P2; TO = T2;
    end
end
j = 2*narcs+1; % load the knot vector *I
for i=O:2
    U(i+1) = 0.0;
    U(i+j+1) = 1.0;
end
switch narcs
    case 1
    case 2
        U(3) = 0.5; U(4)= 0.5;
    case 3
        U(3) = 1.0/3.0; U(4) = 1.0/3.0;
        U(5) = 2.0/3.0; U(6) = 2.0/3.0;
    case 4
        U(3) = 0.25; U(4) = 0.25;
        U(5) = 0.5; U(6) = 0.5;
        U(7) = 0.75; U(8) = 0.75;
end
end

function [P1] = Intersect3DLines(P0,T0,P2,T2)
v1 = T0;
v2 = T2;
alpha = ( v2(2)*P0(1) - v2(2)*P2(1) - v2(1)*P0(2) + v2(1)*P2(2) )...
    /( v2(1)*v1(2)-v2(2)*v1(1));
P1 =  P0 + alpha*v1;
end

