%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ S ] = SurfacePoint( n,p,U,m,q,V,P,u,v )
% NURBS-Book (algorithm A4.3) (modified)
% calculate surface point
%-----------------------------------------------------
dim = size(P,3);
uspan = FindSpan( n,p,u,U );
vspan = FindSpan( m,q,v,V );
Nu = BasisFuns( uspan,u,p,U );
Nv = BasisFuns( vspan,v,q,V );
temp=zeros(1,q+1,dim);
for l=0:q
    for k=0:p
        temp(1,l+1,:)=temp(1,l+1,:)+Nu(k+1)*P(uspan-p+k+1,vspan-q+l+1,:);
    end
end
S = zeros(1,1,dim);
for l=0:q
    S = S+Nv(l+1)*temp(1,l+1,:);
end

end

