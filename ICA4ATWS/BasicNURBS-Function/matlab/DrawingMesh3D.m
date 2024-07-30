%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawingMesh3D( n, p, U, m, q, V, Pw )
DIM = size(Pw,3);

Unique_U = unique(U(p+1:end-p));
Unique_V = unique(V(q+1:end-q));

u=linspace(U(p+1),U(end-p),101);
v=linspace(V(q+1),V(end-q),101);
SIZEU = size(u,2);
SIZEV = size(v,2);
SIZEP = DIM-1;

for i = 1:size(Unique_V,2)
    S = zeros(SIZEU,1,SIZEP);
    for j = 1:size(u,2)
        S(j,1,1:DIM-1)  = SurfacePoint_w( n,p,U,m,q,V,Pw,u(j),Unique_V(i) );
    end
    x = S(:,:,1);
    y = S(:,:,2);
    z = S(:,:,3);
    plot3(x,y,z,'-','LineWidth',1.5);
    hold on
end

for i = 1:size(Unique_U,2)
    S = zeros(1,SIZEV,SIZEP);
    for j = 1:size(v,2)
        S(1,j,1:DIM-1)  = SurfacePoint_w( n,p,U,m,q,V,Pw,Unique_U(i),v(j) );
    end
    x = S(:,:,1);
    y = S(:,:,2);
    z = S(:,:,3);
    plot3(x,y,z,'-','LineWidth',1.5);
    hold on
end
axis equal
hold off
end

