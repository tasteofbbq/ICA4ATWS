%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ C ] = CurvePoint( n,p,U,P,u )
%-----------------------------------------------------------------
% NURBS-Book (algorithm A3.1) (modified)
% calculate point on a curve
%INPUT:
% n         : number ob basis functions -1 !  - x-direction
% p          : degree of the basis functions - x-direction
% U          : knotvector - x-direction
% P          : control points
% u          : x-coordinate
%OUTPUT:
% C          : coordinates of the point on the Curve
%--------------------------------------------------------------
span =  FindSpan(n,p,u,U);
N = BasisFuns(span,u,p,U);
C=zeros(1,2);
for i=1:p+1
    C = C+N(i)*P(span-p+i,:);
end

end

