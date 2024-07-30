%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ S ] = SurfacePoint_w( n,p,U,m,q,V,Pw,u,v )
% NURBS-Book (algorithm A4.3) (modified)
% calculate surface point
%-----------------------------------------------------
dim = size(Pw,3);
uspan = FindSpan( n,p,u,U );
vspan = FindSpan( m,q,v,V );
Nu = BasisFuns( uspan,u,p,U );
Nv = BasisFuns( vspan,v,q,V );
temp=zeros(1,q+1,dim);
for l=0:q
    for k=0:p
        temp(1,l+1,:)=temp(1,l+1,:)+Nu(k+1)*Pw(uspan-p+k+1,vspan-q+l+1,:);
    end
end
Sw = zeros(1,1,dim);
for l=0:q
    Sw = Sw+Nv(l+1)*temp(1,l+1,:);
end
S=Sw(1,1,1:dim-1)/Sw(1,1,dim);

end

