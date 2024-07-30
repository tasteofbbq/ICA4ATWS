function [ C ] = CurvePoint_w( n,p,U,Pw,u )
%-----------------------------------------------------------------
%   Detailed explanation goes here
%function [ C ] = CurvePoint( n,p,U,Pw,u )
% NURBS-Book (algorithm A4.1) (modified)
% calculate point on a curve
%INPUT:
% n         : number ob basis functions -1 !  
% p          : degree of the basis functions
% U          : knotvector 
% P          : control points with weight
% u          : x-coordinate
%OUTPUT:
% C          : coordinates of the point on the Curve
%--------------------------------------------------------------
span =  FindSpan(n,p,u,U);
N = BasisFuns(span,u,p,U);
dim=size(Pw,2);
C=zeros(1,dim);
for i=1:p+1
    C = C+N(i)*Pw(span-p+i,:);
end
C=C(1:dim-1)/C(dim);
end

