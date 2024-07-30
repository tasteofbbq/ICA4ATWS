function [v1s,v1t,v2s,v2t,v3s,v3t] = dVdst(SKL)
%     SKL = [ 0  v   v2
%             u  uv  uv2
%             u2 u2v u2v2]
Ss = [SKL(2,1,1);SKL(2,1,2);SKL(2,1,3)];
St = [SKL(1,2,1);SKL(1,2,2);SKL(1,2,3)];
Sss = [SKL(3,1,1);SKL(3,1,2);SKL(3,1,3)];
Sst = [SKL(2,2,1);SKL(2,2,2);SKL(2,2,3)];
Stt = [SKL(1,3,1);SKL(1,3,2);SKL(1,3,3)];

v3 = cross(Ss,St)/norm(cross(Ss,St));
v1 = Ss/norm(Ss);
v2 = cross(v3,v1)/norm(cross(v3,v1));

%% v3s
A = cross(Ss,St);
derA = cross(Sss,St) + cross(Ss,Sst);
% A = 3*1   derA = 3*1
v3s = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
%% v3t
A = cross(Ss,St);
derA = cross(Sst,St) + cross(Ss,Stt);
% A = 3*1   derA = 3*1
v3t = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
%% v1s
A = Ss;
derA = Sss;
% A = 3*1   derA = 3*1
v1s = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
%% v1t
A = Ss;
derA = Sst;
% A = 3*1   derA = 3*1
v1t = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
%% v2s
A = cross(v3,v1);
derA = cross(v3s,v1) + cross(v3,v1s);
% A = 3*1   derA = 3*1
v2s = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
%% v2t
A = cross(v3,v1);
derA = cross(v3t,v1) + cross(v3,v1t);
% A = 3*1   derA = 3*1
v2t = (derA*sqrt(A'*A)-(0.5*A*(A'*A)^(-0.5)*(derA'*A+A'*derA)))/(A'*A);
end