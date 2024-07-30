function [ element_CP ] = elementCP_Surface( n,p,U,m,q,V,center,noElement )

% n = nurbsSurf.number(1,1); m = nurbsSurf.number(1,2);
% p = nurbsSurf.order(1,1); q = nurbsSurf.order(1,2);
% U = nurbsSurf.knots{1,1}; V = nurbsSurf.knots{1,2};

element_CP = zeros(noElement,(p+1)*(q+1));

for i=1:noElement
    U_span = FindSpan(n,p,center{1,i}(1,1),U);
    V_span = FindSpan(m,q,center{1,i}(1,2),V);
    ID = relavant_CP( U_span, V_span, n, p, m, q );
    element_CP(i,:) = ID;
end
end