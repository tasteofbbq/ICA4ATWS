function [ SKL ] = RatSurfaceDerivs( nurbsSurf,u,v,d )
%     SKL = [ 0  v   v2
%             u  uv  uv2
%             u2 u2v u2v2]
n = nurbsSurf.number(1,1) - 1 ; m = nurbsSurf.number(1,2) - 1;
p = nurbsSurf.order(1,1) - 1  ; q = nurbsSurf.order(1,2) - 1; 
U = nurbsSurf.knots{1,1}      ; V = nurbsSurf.knots{1, 2};
Pw = permute(nurbsSurf.coefs,[2 3 1]);

[ SKLw ] = SurfaceDerivsAlg1( n,p,U,m,q,V,Pw,u,v,d );
Aders = SKLw(:,:,1:end-1);
wders = SKLw(:,:,end);
SKL = zeros(d+1,d+1,size(Pw,3)-1);
[ Bin ] = Binomial(d);
for k=0:d
    for l=0:d-k
        v_ = Aders(k+1,l+1,:);
        for j=1:l
            v_ = v_ - Bin(l+1,j+1)*wders(1,j+1)*SKL(k+1,l-j+1,:);
        end
        for i=1:k
            v_ = v_ - Bin(k+1,i+1)*wders(i+1,1)*SKL(k-i+1,l+1,:);
            v2 = 0.0;
            for j=1:l
               v2 = v2 + Bin(l+1,j+1)*wders(i+1,j+1)*SKL(k-i+1,l-j+1,:); 
            end
            v_ = v_ - Bin(k+1,i+1)*v2;
        end
        SKL(k+1,l+1,:) = v_/wders(1,1);
    end
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

