function drawEigenmodes2D( WWW, NN, N, D, L,n,p,U,m,q,V )
u=0:0.01:1;v=u;
[U_u,V_v] = meshgrid(u,v);
WWWW = zeros(numel(u),numel(v));
for mode = 1:size(WWW,2)
for i=1:numel(u)
    for j=1:numel(v)
WWWW(i,j) = SurfacePoint(n,p,U,m,q,V,WWW(:,mode),u(i),v(j));
    end
end
figure
colormap(jet)
surf(U_u*L,V_v*L,WWWW);
shading interp,axis ( [0,L,0,L,-1,1] ) 
end
end

