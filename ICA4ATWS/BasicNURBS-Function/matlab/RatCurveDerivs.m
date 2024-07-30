function [ CK ] = RatCurveDerivs( n,p,U,Pw,u,d )
% NURBS-Book (algorithm A4.2) (modified)
% calculate derivatives of a curve
%INPUT:
% n          : number ob basis functions -1 !  - x-direction
% p          : degree of the basis functions - x-direction
% U          : knotvector - x-direction
% P          : control points
% u          : x-coordinate
% d          : 
%OUTPUT:
% CK         : derivatives of non-rational B-spline curve
[ CKw ] = CurvedDerivsAlg1( n,p,U,Pw,u,d );
wders = CKw(:,end);
Aders = CKw(:,1:end-1);
CK = zeros(d+1,size(Pw,2)-1);
[ Bin ] = Binomial(d);
for k=0:d
    v = Aders(k+1,:);
    for i=1:k
        v = v - Bin(k+1,i+1)*wders(i+1)*CK(k-i+1,:);
    end
    CK(k+1,:) = v/wders(1);
end
end

function [ Bin ] = Binomial(n)
%  Purpose: caculate binominal
Bin=zeros(n+1,n+1);
for i=0:n
    for j=0:i
        if (j==0 || j==i) 
        Bin(i+1,j+1)=1;
        else
            Bin(i+1,j+1) = Bin(i,j) + Bin(i,j+1);
        end 
    end    
end
end

