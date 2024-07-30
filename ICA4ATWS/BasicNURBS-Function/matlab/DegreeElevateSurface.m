function [ nh,Uh,mh,Vh,Qw ] = DegreeElevateSurface( n,p,U,m,q,V,Pw,t1,t2 )
%
%  Record of revision:
%      Date      Programmer        Description of change
%      ====      ==========        =====================
%     15/11/30                    
%  Degree elevate a curve t times
%  Input: n,p,U,Pw,t 
%  Output: nh,Uh,Qw 
%  Data dictionary: declare variable type & definition
DIM = size(Pw,3)-1;
L_U = length(unique(U));
L_V = length(unique(V));
Uh = zeros(1,n+p+2+L_U*t1);
Vh = zeros(1,m+q+2+L_V*t2);
tempQw = zeros(n+(L_U-1)*t1+1,m+1,DIM+1);
Qw = zeros(n+(L_U-1)*t1+1,m+(L_V-1)*t2+1,DIM+1);

for ii=1:m+1
bezalfs = zeros(p+t1+1,p+1);
alfs = zeros(p-1,1);
Nextbpts = zeros(p-1,DIM+1);
ebpts = zeros(p+t1+1,DIM+1);
bpts = zeros(p+1,DIM+1);

[ Bin ] = binomial(p+t1);


noCPsX = n+p+1; ph = p+t1; ph2 = fix(ph/2);
%  Compute Bezier degree elevation coefficients
bezalfs(1,1) = 1.0;
bezalfs(ph+1,p+1) = 1.0;

for i=1:ph2
    inv = 1.0/Bin(ph+1,i+1);
    mpi = min(p,i);
    for j=max(0,i-t1):mpi
        bezalfs(i+1,j+1) = inv*Bin(p+1,j+1)*Bin(t1+1,i-j+1);
    end
end

for i=ph2+1:ph-1
    mpi = min(p,i);
    for j=max(0,i-t1):mpi
        bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);
    end 
end

nph = ph; kind_ = ph+1; r=-1; a=p;
b = p+1; cind = 1; ua = U(1);
tempQw (1,ii,1:DIM+1) = Pw (1,ii,1:DIM+1);
for i=0:ph
    Uh(i+1) = ua;
end
% Initialize first Bezier seg
for i=0:p
    bpts(i+1,1:DIM+1) = Pw(i+1,ii,1:DIM+1);
end

while (b < noCPsX)  %  Big loop thru knot vector
    i = b;
    while (b < noCPsX && U(b+1) == U(b+2)) 
        b = b+1;
    end 
    mul = b-i+1;
    nph = nph+mul+t1;
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
        for j=max(0,i-t1):mpi
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
                    tempQw(i+1,ii,1:DIM+1) = alf*tempQw(i+1,ii,1:DIM+1) + (1.0-alf)*tempQw(i,ii,1:DIM+1);
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
            tempQw(cind+1,ii,1:DIM+1) = ebpts(j+1,1:DIM+1);
            cind = cind+1;
    end 
    if (b < noCPsX)   % Set up for next
        for j=0:r-1
            bpts(j+1,1:DIM+1) = Nextbpts(j+1,1:DIM+1);
        end
        for j=r: p
            bpts(j+1,1:DIM+1) = Pw(b-p+j+1,ii,1:DIM+1);
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
nh = nph - ph -1;

end


for ii=1:nh+1
bezalfs = zeros(q+t2+1,q+1);
alfs = zeros(q-1,1);
Nextbpts = zeros(q-1,DIM+1);
ebpts = zeros(q+t2+1,DIM+1);
bpts = zeros(q+1,DIM+1);

[ Bin ] = binomial(q+t2);


noCPsY = m+q+1; qh = q+t2; qh2 = fix(qh/2);
%  Compute Bezier degree elevation coefficients
bezalfs(1,1) = 1.0;
bezalfs(qh+1,q+1) = 1.0;

for i=1:qh2
    inv = 1.0/Bin(qh+1,i+1);
    mpi = min(q,i);
    for j=max(0,i-t2):mpi
        bezalfs(i+1,j+1) = inv*Bin(q+1,j+1)*Bin(t2+1,i-j+1);
    end
end

for i=qh2+1:qh-1
    mpi = min(q,i);
    for j=max(0,i-t2):mpi
        bezalfs(i+1,j+1) = bezalfs(qh-i+1,q-j+1);
    end 
end

mqh = qh; kind_ = qh+1; r=-1; a=q;
b = q+1; cind = 1; va = V(1);
Qw (ii,1,1:DIM+1) = tempQw (ii,1,1:DIM+1);
for i=0:qh
    Vh(i+1) = va;
end
% Initialize first Bezier seg
for i=0:q
    bpts(i+1,1:DIM+1) = tempQw(ii,i+1,1:DIM+1);
end

while (b < noCPsY)  %  Big loop thru knot vector
    i = b;
    while (b < noCPsY && V(b+1) == V(b+2)) 
        b = b+1;
    end 
    mul = b-i+1;
    mqh = mqh+mul+t2;
    vb = V(b+1);
    oldr = r;
    r = q-mul;
    
    %  Insert knot u(b) r times 
    if (oldr > 0) 
        lbz = (oldr+2)/2;
    else
        lbz = 1;
    end
    
    if (r > 0) 
        rbz = qh-(r+1)/2;
    else
        rbz = qh;
    end 
    
    if (r > 0) 
        %  Insert knot to get Bezier segment
        numer = vb-va;
        for k=q:mul+1
            alfs(k-mul) = numer/(V(a+k+1)-va);
        end 
        for j=1:r
            saved= r-j;
            s = mul+j;
            for k=q:s
                bpts(k+1,1:DIM+1) = alfs(k-s+1)*bpts(k+1,1:DIM+1) + (1.0-alfs(k-s+1))*bpts(k,1:DIM+1);
            end
            Nextbpts(saved+1,1:DIM+1) = bpts(q+1,1:DIM+1);
        end
    end   %  End of "insert knot"
    for i=lbz:qh  %  Degree elevate Bezier *I
        %Only points lbz, ... ,ph are used below *I
        ebpts(i+1,1:DIM+1) = 0.0;
        mpi = min(q,i);
        for j=max(0,i-t2):mpi
            ebpts(i+1,1:DIM+1) = ebpts(i+1,1:DIM+1) + bezalfs(i+1,j+1)*bpts(j+1,1:DIM+1);
        end 
    end   % End of degree elevating Bezier *I
    if (oldr > 1) 
        %  Must remove knot u=U(a) oldr times *I
        first = kind_-2;
        last = kind_;
        den = vb-va;
        bet= (vb-Vh(kind_))/den;
        for tr=1:oldr-1
            %  Knot removal loop *I
            i = first;
            j = last;
            kj = j-kind_+1;
            while (j-i > tr)  %  Loop and compute the new *I
                %  control points for one removal step *I
                if (i < cind) 
                    alf = (vb-Vh(i+1))/(va-Vh(i+1));
                    Qw(ii,i+1,1:DIM+1) = alf*Qw(ii,i+1,1:DIM+1) + (1.0-alf)*Qw(ii,i,1:DIM+1);
                end 
                if ( j >= lbz) 
                    if ( j-tr <= kind_-qh+oldr ) 
                        gam= (vb-Vh(j-tr+1))/den;
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
    if (a ~= q)   %  Load the knot ua 
        for i=0:qh-oldr-1
            Vh(kind_+1) = va;
            kind_ = kind_+1;
        end 
    end 
    for j=lbz: rbz  %  Load ctrl pts into Qw
            Qw(ii,cind+1,1:DIM+1) = ebpts(j+1,1:DIM+1);
            cind = cind+1;
    end 
    if (b < noCPsY)   % Set up for next
        for j=0:r-1
            bpts(j+1,1:DIM+1) = Nextbpts(j+1,1:DIM+1);
        end
        for j=r: q
            bpts(j+1,1:DIM+1) = tempQw(ii,b-q+j+1,1:DIM+1);
        end 
        a=b;
        b=b+1;
        va=vb;
    else
        %  End knot
        for i=0:qh
            Vh(kind_+i+1) = vb;
        end 
    end 
end   %   End of while-loop (b < m)
mh = mqh - qh -1;

end

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

