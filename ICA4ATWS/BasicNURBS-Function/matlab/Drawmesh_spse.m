function Drawmesh_spse( ieplane,F_geometry)
if(F_geometry~=0)
    n = ieplane.number(1,1) - 1 ; m = ieplane.number(1,2) - 1;
    p = ieplane.order(1,1) - 1  ; q = ieplane.order(1,2) - 1;
    U = ieplane.knots{1,1}; V = ieplane.knots{1, 2};
    Pw = permute(ieplane.coefs,[2 3 1]);
    figure(F_geometry)
    w = Pw(:,:,4);
    CP_x = Pw(:,:,1)./w;
    CP_y = Pw(:,:,2)./w;
    CP_z = Pw(:,:,3)./w;
    plot3(CP_x,CP_y,CP_z,'.','color','r','MarkerSize',20);
    hold on
    plot3(CP_x,CP_y,CP_z,'--','color','k','LineWidth',1.5);
    hold on
    plot3(CP_x',CP_y',CP_z','--','color','k','LineWidth',1.5);
    axis equal
    xlabel('X');ylabel('Y');zlabel('Z')
    set(0,'defaultfigurecolor','w')
    set(gca,'DataAspectRatio',[1 1 1],'color',[1,1,1])
    DrawingMesh3D( n, p, U, m, q, V, Pw );
end
end

