function [ CK ] = CurvedDerivsAlg1( curve,u,d )
%   CK = [ 0  0  0
%          u  u  u
%          u2 u2 u2]
n = curve.number - 1;
p = curve.order - 1;
U = curve.knots;
P = permute(curve.coefs,[2,1]);
% Algotithm A3.2
CK=zeros(1,size(P,2));
du=min(d,p);
for k=p+1:d
    CK(k,:)=0;
end
span = FindSpan(n,p,u,U);
nders = DersBasisFuns(span,u,p,du,U);
for i=0:du
    CK(i+1,:) = 0.0;
    for j=0:p
        CK(i+1,:) = CK(i+1,:) +nders(i+1,j+1)*P(span-p+1+j,:);
    end
end

if (p==1&& d>=2)
    CK(3:d+1,:) = 0;
end

end