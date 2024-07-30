function [ n, p, U, m, q, V, Pw, w ] = p_h_Refine_for_Surf( n, p, U, m, q,...
    V, Pw, t1, t2,newKnotsU,newKnotsV )

%% p-refine
[ nh,Uh,mh,Vh,Qw ] = DegreeElevateSurface( n,p,U,m,q,V,Pw,t1,t2 );
n = nh; m =mh; U = Uh; V = Vh; Pw = Qw;
clear nh;clear mh;
clear Qw; clear Uh; clear Vh;
p = p+t1; q = q + t2;

w = Pw(:,:,end);
CP_x = Pw(:,:,1)./w; CP_y = Pw(:,:,2)./w; CP_z = Pw(:,:,3)./w;

figure
plot3(CP_x,CP_y,CP_z,'.','color','r','MarkerSize',25);
hold on
plot3(CP_x,CP_y,CP_z,'--','color','k','LineWidth',1.5);
hold on
plot3(CP_x',CP_y',CP_z','--','color','k','LineWidth',1.5);
axis equal

% DrawingMesh3D( n, p, U, m, q, V, Pw );

%% h-Refine

% U direction
r = numel(newKnotsU)-1;
[ Ubar,Vbar,Qw ] = RefineKnotVectSurface( n,p,U,m,q,V,Pw,newKnotsU,r,'U_DIRECTION' );
n=n+r+1;
clear U;clear Pw;
U = Ubar;
Pw = Qw;
clear Qw; 

% V direction
r = numel(newKnotsV)-1;
[ Ubar,Vbar,Qw ] = RefineKnotVectSurface( n,p,U,m,q,V,Pw,newKnotsV,r,'V_DIRECTION' );
m=m+r+1;
clear V;clear Pw
V = Vbar;
Pw = Qw;
clear Qw;clear Ubar;clear Vbar;
clear CP_x;clear CP_y;clear CP_z;clear w; 

% figure
 w = Pw(:,:,4);
CP_x = Pw(:,:,1)./w;
CP_y = Pw(:,:,2)./w;
CP_z = Pw(:,:,3)./w;
plot3(CP_x,CP_y,CP_z,'.','color','r','MarkerSize',18);
hold on
plot3(CP_x,CP_y,CP_z,'--','color','k','LineWidth',1.5);
hold on
plot3(CP_x',CP_y',CP_z','--','color','k','LineWidth',1.5);
axis equal

figure
% DrawingMesh3D( n, p, U, m, q, V, Pw );


end

