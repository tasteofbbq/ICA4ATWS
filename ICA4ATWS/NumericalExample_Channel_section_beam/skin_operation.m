%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Zhengyang Zhang <zhangzhengyang4top@outlook.com>
% Date:   July 29, 2024
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    The GNU General Public License, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Skin ] = skin_operation(nurbsSurf,Nu,Nv,F_geometry,F_parameter)
%% P-refine
pre_p_order = 3;
pre_q_order = 3;
order_u = nurbsSurf.order(1);
order_v = nurbsSurf.order(2);
p = pre_p_order-order_u;
q = pre_q_order-order_v;
if(p<=0)
    p = 0;
end
if(q<=0)
    q = 0;
end
nurbsSurf = nrbdegelev(nurbsSurf, [p q]);

%% H-refine
nurbsSurf = H_Refine(nurbsSurf,Nu,Nv);
ieplane = nurbsSurf;

%% P-refine
n = ieplane.number(1,1) - 1 ; m = ieplane.number(1,2) - 1;
p = ieplane.order(1,1) - 1  ; q = ieplane.order(1,2) - 1; 
U = ieplane.knots{1,1}; V = ieplane.knots{1, 2};
Pw = permute(ieplane.coefs,[2 3 1]);

%% H-refine
Drawmesh_spse( ieplane,F_geometry )

%% generate skin
Skin.Pw1 = Pw;
Pw  = permute(Pw,[3,1,2]);
Skin.ieplane.coefs  = Pw;
Skin.ieplane.order  = [p+1,q+1];
Skin.ieplane.number = [n+1,m+1];
Skin.ieplane.knots  = [{U},{V}];
Skin.ieplane.dim  = 4;
Skin.ieplane.unique_knots{1} = unique(Skin.ieplane.knots{1});
Skin.ieplane.unique_knots{2} = unique(Skin.ieplane.knots{2});

%% pre-processing
Skin.Ncops = Skin.ieplane.number(1,1) * Skin.ieplane.number(1,2);
%  Number of DOFs in the whole model
Skin.Gdofs = Skin.Ncops*5;
disp(Skin.Gdofs)
% Homogeneous coordinates are transformed into Cartesian coordinates
Skin.conpsX = reshape(Skin.ieplane.coefs(1,:,:)./Skin.ieplane.coefs(4,:,:),Skin.Ncops,1);
Skin.conpsY = reshape(Skin.ieplane.coefs(2,:,:)./Skin.ieplane.coefs(4,:,:),Skin.Ncops,1);
Skin.conpsZ = reshape(Skin.ieplane.coefs(3,:,:)./Skin.ieplane.coefs(4,:,:),Skin.Ncops,1);

%% control point
[ Skin.CP ] = ControlPoint(Skin.conpsX,Skin.conpsY,Skin.conpsZ,Skin.Ncops);
[ Skin.CPw ] = Project_Surface_greville( Skin );
[ Skin.elementNode , Skin.noElement ] = GetParaMesh_Surface( U,V,p,q,F_parameter);
Skin.elementNode = Skin.elementNode; % Coordinates of the parameter domain for the center point of each element
[ Skin.element_CP ] = elementCP_Surface(  n,p,U,m,q,V,Skin.elementNode.center,Skin.noElement );
[ Skin ] = Local_Coordinate_Surface( Skin );

%% get weight of control point
coefs_temp = reshape(Skin.ieplane.coefs,Skin.ieplane.dim,[])';
Skin.CP_wghts = coefs_temp(:,end);

%% Cycle each element
knots = Skin.ieplane.knots;
order = Skin.ieplane.order;
uknots_uni = unique(knots{1}(order(1):end-order(1)+1));
vknots_uni = unique(knots{2}(order(2):end-order(2)+1));
EL_number = [length(uknots_uni),length(vknots_uni)]-1;
count = 0;
for index_v = 1:EL_number(2)
    for index_u = 1:EL_number(1)
        count = count+1;
        % Element boundary knot
        u1 = uknots_uni(index_u); u2 = uknots_uni(index_u+1);
        v1 = vknots_uni(index_v); v2 = vknots_uni(index_v+1);
        % Element center coordinates in the parameter space
        EL_center(:,count) = [0.5*(u1+u2);0.5*(v1+v2);0];
    end
end

%% Element control point indices
number = length(knots{1})-order(1);
uspans = FindSpans(knots{1},order(1)-1,EL_center(1,:));
vspans = FindSpans(knots{2},order(2)-1,EL_center(2,:));
initId = reshape((repmat((1:order(1)).',1,order(2))+number*(0:order(2)-1)),1,[]);
Skin.EL_nodeId = initId+number*(vspans-order(2))+uspans-order(1);