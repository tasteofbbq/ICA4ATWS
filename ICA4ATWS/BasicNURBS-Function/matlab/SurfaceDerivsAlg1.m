function [ SKL ] = SurfaceDerivsAlg1( n,p,U,m,q,V,P,u,v,d )
%  ALGORITHM A3.6
%  Compute surface derivatives
%  Input: n,p,U,m,q,V,P,u,v,d 
%  Output: SKL 
DIM = size(P,3);
SKL =zeros(d+1,d+1,DIM);
temp = zeros(1,p+1,DIM);
du = min(d,p);
for k=p+1: d
    for l=0: d-k
        SKL(k+1,l+1,1:DIM) = 0.0;
    end 
end
dv = min(d,q);
for l=q+1: d
    for k=0: d-l
        SKL(k+1,l+1,1:DIM) = 0.0;
    end
end
uspan = FindSpan(n,p,u,U);
Nu = DersBasisFuns( uspan,u,p,du,U );
vspan = FindSpan( m,q,v,V );
Nv = DersBasisFuns( vspan,v,q,dv,V );
for k=0: du
    for s=0: q
        temp(1,s+1,1:DIM) = 0.0;
        for r=0: p
            temp(1,s+1,1:DIM) = temp(1,s+1,1:DIM) + Nu(k+1,r+1)*P(uspan-p+r+1,vspan-q+s+1,1:DIM);
        end 
    end 
    dd = min(d-k,dv);
    for l=0: dd
        SKL(k+1,l+1,1:DIM) = 0.0;
        for s=0: q
            SKL(k+1,l+1,1:DIM) = SKL(k+1,l+1,1:DIM) + Nv(l+1,s+1)*temp(1,s+1,1:DIM);
        end 
    end 
end 

end

