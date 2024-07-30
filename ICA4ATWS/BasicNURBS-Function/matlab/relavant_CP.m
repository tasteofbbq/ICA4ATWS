function [ ID ] = relavant_CP( U_span, V_span ,n, p, m, q )
ID = zeros(1,(p+1)*(q+1));

for i = 0:q
    for j=0:p
        i_p=U_span-p; i_q=V_span-q;
        ID(1,j+1+i*(p+1)) = (i+i_q)*(n+1)+j+1+i_p;
    end
end

end

