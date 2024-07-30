function [ Nip ] = OneBasisFun( p,m,U,i,u )
%  Data dictionary: declare local variable
NN = zeros(1,p+1);

%  MAIN
if ( (i == 1 && u == U(1)) || (i == m - p && u == U(m+1)) )
    Nip = 1.0;
    return
end

if ( u < U(i) || u >= U(i + p + 1) ) 
    Nip = 0.0;
    return
end 

for j = 0:p
    if ( u >= U(i + j ) && u < U(i + j + 1) ) 
        NN(1,j+1) = 1.0;
    else
        NN(1,j+1) = 0.0;
    end 
end 

for k = 1:p
    if ( NN(1,1) == 0.0 ) 
        saved = 0.0;
    else
        saved = ((u - U(i))*NN(1,1)) / (U(i + k) - U(i));
    end
    
    for  j = 0:p - k
        Uleft = U(i + j + 1);
        Uright = U(i + j + k + 1);
        if ( NN(1,j + 2) == 0.0 )
            NN(1,j+1) = saved;
            saved = 0;
        else
            temp = NN(1,j + 2) / (Uright - Uleft);
            NN(1,j+1) = saved + (Uright - u)*temp;
            saved = (u - Uleft)*temp;
        end
    end
		
end
Nip = NN(1,1);
end

