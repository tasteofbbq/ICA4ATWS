function [ieplane] = H_Refine(ieplane,Nu,Nv)
%% U direction
L = 1;
nkntsU = linspace(ieplane.knots{1,1}(1),ieplane.knots{1,1}(end),Nu);
for i = 1:size(nkntsU,2)
    lib = abs(ieplane.knots{1,1}-nkntsU(i));
    if min(lib)>0.001 
        newKnotsU(L) = nkntsU(i);
        L = L+1;
    end
end

%% V direction
L = 1;
nkntsV = linspace(ieplane.knots{1,2}(1),ieplane.knots{1,2}(end),Nv);
for i = 1:size(nkntsV,2)
    lib = abs(ieplane.knots{1,2}-nkntsV(i));
    if min(lib)>0.001 
        newKnotsV(L) = nkntsV(i);
        L = L+1;
    end
end

%% 
iknots{1,1} = newKnotsU;
iknots{1,2} = newKnotsV;
ieplane = nrbkntins(ieplane,iknots);
end