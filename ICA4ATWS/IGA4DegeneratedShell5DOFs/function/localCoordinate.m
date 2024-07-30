function [v1, v2 ,v3 ,theta] = localCoordinate(SKL)
% shell local coordinate system
Ss = [SKL(2,1,1);SKL(2,1,2);SKL(2,1,3)];
St = [SKL(1,2,1);SKL(1,2,2);SKL(1,2,3)];
v3 = cross(Ss,St)/norm(cross(Ss,St));
v1 = Ss/norm(Ss);
v2 = cross(v3,v1)/norm(cross(v3,v1));
theta =  [v1' ; v2' ; v3'];
end