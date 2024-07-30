function [ nh,Uh,Qw ] = DegreeElevateCurve( n,p,U,Pw,t )


%  Record of revision:
%      Date      Programmer        Description of change
%      ====      ==========        =====================
%     15/5/24                    U to UU
%
%  Degree elevate a curve t times
%  Input: n,p,U,Pw,t 
%  Output: nh,Uh,Qw 
%  Data dictionary: declare variable type & definition
DIMENSION = size(Pw);
DIM = DIMENSION(end)-1;

Uh = zeros(1,n+p+1+(n-p+2)*t+1);
Qw = zeros(n+(n-p+1)*t+1,DIM+1);

bezalfs = zeros(p+t+1,p+1);
alfs = zeros(p-1,1);
Nextbpts = zeros(p-1,DIM+1);
ebpts = zeros(p+t+1,DIM+1);
bpts = zeros(p+1,DIM+1);

[ Bin ] = binomial(p+t);


m = n+p+1; ph = p+t; ph2 = fix(ph/2);
%  Compute Bezier degree elevation coefficients
bezalfs(1,1) = 1.0;
bezalfs(ph+1,p+1) = 1.0;

for i=1:ph2
    inv = 1.0/Bin(ph+1,i+1);
    mpi = min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = inv*Bin(p+1,j+1)*Bin(t+1,i-j+1);
    end
end

for i=ph2+1:ph-1
    mpi = min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);
    end 
end

mh = ph; kind_ = ph+1; r=-1; a=p;
b = p+1; cind = 1; ua = U(1);
Qw (1,1:DIM+1) = Pw (1,1:DIM+1);
for i=0:ph
    Uh(i+1) = ua;
end
% Initialize first Bezier seg
for i=0:p
    bpts(i+1,1:DIM+1) = Pw(i+1,1:DIM+1);
end

while (b < m)  %  Big loop thru knot vector
    i = b;
    while (b < m && U(b+1) == U(b+2)) 
        b = b+1;
    end 
    mul = b-i+1;
    mh = mh+mul+t;
    ub = U(b+1);
    oldr = r;
    r = p-mul;
    
    %  Insert knot u(b) r times 
    if (oldr > 0) 
        lbz = (oldr+2)/2;
    else
        lbz = 1;
    end
    
    if (r > 0) 
        rbz = ph-(r+1)/2;
    else
        rbz = ph;
    end 
    
    if (r > 0) 
        %  Insert knot to get Bezier segment
        numer = ub-ua;
        for k=p:mul+1
            alfs(k-mul) = numer/(U(a+k+1)-ua);
        end 
        for j=1:r
            saved= r-j;
            s = mul+j;
            for k=p:s
                bpts(k+1,1:DIM+1) = alfs(k-s+1)*bpts(k+1,1:DIM+1) + (1.0-alfs(k-s+1))*bpts(k,1:DIM+1);
            end
            Nextbpts(saved+1,1:DIM+1) = bpts(p+1,1:DIM+1);
        end
    end   %  End of "insert knot"
    for i=lbz:ph  %  Degree elevate Bezier *I
        %Only points lbz, ... ,ph are used below *I
        ebpts(i+1,1:DIM+1) = 0.0;
        mpi = min(p,i);
        for j=max(0,i-t):mpi
            ebpts(i+1,1:DIM+1) = ebpts(i+1,1:DIM+1) + bezalfs(i+1,j+1)*bpts(j+1,1:DIM+1);
        end 
    end   % End of degree elevating Bezier *I
    if (oldr > 1) 
        %  Must remove knot u=U(a) oldr times *I
        first = kind_-2;
        last = kind_;
        den = ub-ua;
        bet= (ub-Uh(kind_))/den;
        for tr=1:oldr-1
            %  Knot removal loop *I
            i = first;
            j = last;
            kj = j-kind_+1;
            while (j-i > tr)  %  Loop and compute the new *I
                %  control points for one removal step *I
                if (i < cind) 
                    alf = (ub-Uh(i+1))/(ua-Uh(i+1));
                    Qw(i+1,1:DIM+1) = alf*Qw(i+1,1:DIM+1) + (1.0-alf)*Qw(i,1:DIM+1);
                end 
                if ( j >= lbz) 
                    if ( j-tr <= kind_-ph+oldr ) 
                        gam= (ub-Uh(j-tr+1))/den;
                        ebpts(kj+1,1:DIM+1) = gam*ebpts(kj+1,1:DIM+1)+(1.0-gam)*ebpts(kj+2,1:DIM+1);
                    else
                        ebpts(kj+1,1:DIM+1) = bet*ebpts(kj+1,1:DIM+1)+(1.0-bet)*ebpts(kj+2,1:DIM+1);
                    end 
                end 
                i = i+1;
                j = j-1;
                kj = kj-1;
            end 
            first = first-1;
            last = last+1;
        end 
    end   %  End of removing knot, u=U(a)
    if (a ~= p)   %  Load the knot ua 
        for i=0:ph-oldr-1
            Uh(kind_+1) = ua;
            kind_ = kind_+1;
        end 
    end 
    for j=lbz: rbz  %  Load ctrl pts into Qw
            Qw(cind+1,1:DIM+1) = ebpts(j+1,1:DIM+1);
            cind = cind+1;
    end 
    if (b < m)   % Set up for next
        for j=0:r-1
            bpts(j+1,1:DIM+1) = Nextbpts(j+1,1:DIM+1);
        end
        for j=r: p
            bpts(j+1,1:DIM+1) = Pw(b-p+j+1,1:DIM+1);
        end 
        a=b;
        b=b+1;
        ua=ub;
    else
        %  End knot
        for i=0:ph
            Uh(kind_+i+1) = ub;
        end 
    end 
end   %   End of while-loop (b < m)
nh = mh - ph -1;

end

function [ Bin ] = binomial(n)
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
