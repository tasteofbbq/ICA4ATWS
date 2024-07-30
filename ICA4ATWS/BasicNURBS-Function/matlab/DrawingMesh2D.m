%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawingMesh2D( n, p, U, m, q, V, Pw )
DIM = size(Pw,3);

Unique_U = unique(U);
Unique_V = unique(V);

u=(0:0.01:1);
v=(0:0.01:1);
uu = (0:0.001:1);
vv = (0:0.001:1);


for i = 1:size(Unique_V,2)
    S = zeros(size(u,2),1,size(Pw,3)-1);
    for j = 1:size(u,2)
        S(j,1,1:DIM-1)  = SurfacePoint_w( n,p,U,m,q,V,Pw,u(j),Unique_V(i) );
    end
    x = spline(u,S(:,:,1),uu);
    y = spline(u,S(:,:,2),uu);
    plot(x,y,'-','LineWidth',1.5);
    hold on
end

for i = 1:size(Unique_U,2)
    S = zeros(1,size(v,2),size(Pw,3)-1);
    for j = 1:size(v,2)
        S(1,j,1:DIM-1)  = SurfacePoint_w( n,p,U,m,q,V,Pw,Unique_U(i),v(j) );
    end
    x = spline(u,S(:,:,1),uu);
    y = spline(u,S(:,:,2),uu);
    plot(x,y,'-','LineWidth',1.5);
    hold on
end

axis equal
hold off



end

