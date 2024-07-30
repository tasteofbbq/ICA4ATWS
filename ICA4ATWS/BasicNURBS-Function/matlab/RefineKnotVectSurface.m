function [ Ubar,Vbar,Qw ] = RefineKnotVectSurface( n,p,U,m,q,V,Pw,X,r,dir )
% NURBS-Book (algorithm A5.5) (modified)
%  Refine surface knot vector
%  Input: n,p,U,m,q,V,Pw,X,r,dir
%  Output: Ubar,Vbar,Qw
%  Data dictionary: declare variable type & definition

dim=size(Pw,3);

if( strcmp(dir,'U_DIRECTION') )
    Qw = zeros(n+r+2,m+1,dim);
    Ubar = zeros(1,n+p+2+r+1);
    %
    a=FindSpan(n,p,X(1),U);
    b=FindSpan(n,p,X(r+1),U);
    b=b+1;
    Vbar=V;  %copy V to Vbar
    for j=0:a-p
        Qw(j+1,:,:) = Pw(j+1,:,:);
    end
    for j=b-1:n
        Qw(j+r+2,:,:) = Pw(j+1,:,:);
    end
    for j=0:a
        Ubar(j+1)= U(j+1);
    end
    for j=b+p :n+p+1
        Ubar(j+r+2) = U(j+1);
    end
    i=b+p-1;
    k=b+p+r;
    for j=r:-1:0
        while (X(j+1) <= U(i+1) && i>a)
            Ubar(k+1) = U(i+1);
            Qw(k-p,:,:) = Pw(i-p,:,:);
            k=k-1;
            i=i-1;
        end
        Qw(k-p,:,:) = Qw(k-p+1,:,:);
        for l=1:p
            ind = k-p+l;
            alfa = Ubar(k+l+1) - X(j+1);
            if (abs(alfa) == 0)
                Qw(ind,:,:) = Qw(ind+1,:,:);
            else
                alfa = alfa/(Ubar(k+l+1) - U(i-p+l+1));
                Qw(ind,:,:) = alfa* Qw(ind,:,:) + (1-alfa)* Qw(ind+1,:,:);
            end
        end
        Ubar(k+1) = X(j+1);
        k=k-1;
    end
end

if( strcmp(dir,'V_DIRECTION') )
    Qw = zeros(n+1,m+r+2,dim);
    Vbar = zeros(1,m+q+2+r+1);
    %
    a=FindSpan(m,q,X(1),V);
    b=FindSpan(m,q,X(r+1),V);
    b=b+1;
    Ubar=U;  %copy U to Ubar
    for j=0:a-q
        Qw(:,j+1,:) = Pw(:,j+1,:);
    end
    for j=b-1:m
        Qw(:,j+r+2,:) = Pw(:,j+1,:);
    end
    for j=0:a
        Vbar(j+1)= V(j+1);
    end
    for j=b+q :m+q+1
        Vbar(j+r+2) = V(j+1);
    end
    i=b+q-1;
    k=b+q+r;
    for j=r:-1:0
        while (X(j+1) <= V(i+1) && i>a)
            Vbar(k+1) = V(i+1);
            Qw(:,k-q,:) = Pw(:,i-q,:);
            k=k-1;
            i=i-1;
        end
        Qw(:,k-q,:) = Qw(:,k-q+1,:);
        for l=1:q
            ind = k-q+l;
            alfa = Vbar(k+l+1) - X(j+1);
            if (abs(alfa) == 0)
                Qw(:,ind,:) = Qw(:,ind+1,:);
            else
                alfa = alfa/(Vbar(k+l+1) - V(i-q+l+1));
                Qw(:,ind,:) = alfa* Qw(:,ind,:) + (1-alfa)* Qw(:,ind+1,:);
            end
        end
        Vbar(k+1) = X(j+1);
        k=k-1;
    end
end


end

