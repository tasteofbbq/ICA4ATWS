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

function [interfacial_stiffness] = interfacial_stiffness_matrices4Mortar_self_tie(Skin_total, temp_information)
patch_index1         = temp_information(1,1);
patch_index2         = temp_information(1,2);
local_surface_index1 = temp_information(2,1);
local_surface_index2 = temp_information(2,2);
NoCP1 = (Skin_total{patch_index1}.ieplane.number(1,1))*(Skin_total{patch_index1}.ieplane.number(1,2));
dof1 = 5*NoCP1;
NoCP2 = (Skin_total{patch_index2}.ieplane.number(1,1))*(Skin_total{patch_index2}.ieplane.number(1,2));
dof2 = 5*NoCP2;

ieplane1  = Skin_total{patch_index1}.ieplane;
down1  = 1:ieplane1.number(1,1);
up1    = (ieplane1.number(1,2)-1)*(ieplane1.number(1,1))+1:(ieplane1.number(1,1))*(ieplane1.number(1,2));
left1  = 1:(ieplane1.number(1,1)):(ieplane1.number(1,2)-1)*(ieplane1.number(1,1))+1;
right1 = (ieplane1.number(1,1)):(ieplane1.number(1,1)):(ieplane1.number(1,1))*(ieplane1.number(1,2));

ieplane2  = Skin_total{patch_index2}.ieplane;
down2  = 1:ieplane2.number(1,1);
up2    = (ieplane2.number(1,2)-1)*(ieplane2.number(1,1))+1:(ieplane2.number(1,1))*(ieplane2.number(1,2));
left2  = 1:(ieplane2.number(1,1)):(ieplane2.number(1,2)-1)*(ieplane2.number(1,1))+1;
right2 = (ieplane2.number(1,1)):(ieplane2.number(1,1)):(ieplane2.number(1,1))*(ieplane2.number(1,2));

switch local_surface_index1
    case 4
        coupling_edge1 = right1;
    case 6
        coupling_edge1 = left1;
    case 2
        coupling_edge1 = up1;
    case 1
        coupling_edge1 = down1;
end
switch local_surface_index2
    case 4
        coupling_edge2 = right2;
    case 6
        coupling_edge2 = left2;
    case 2
        coupling_edge2 = up2;
    case 1
        coupling_edge2 = down2;
end

Ks_lower_left  = sparse(length(coupling_edge1)*5,dof1);
Ks_lower_right = sparse(length(coupling_edge1)*5,dof2);
for index_CE = 1:length(coupling_edge1)
    index_CEP1 = coupling_edge1(index_CE);
    index_CEP2 = coupling_edge2(index_CE);
    Ks_lower_left(index_CE*5-4:index_CE*5,index_CEP1*5-4:index_CEP1*5) = eye(5);
    Ks_lower_right(index_CE*5-4:index_CE*5,index_CEP2*5-4:index_CEP2*5) = eye(5);
end

Ks_lower_left  = sparse(Ks_lower_left);
Ks_lower_right = sparse(Ks_lower_right);

if(patch_index1<=patch_index2)
    Ks = [-Ks_lower_left  Ks_lower_right];
else
    Ks = [Ks_lower_right -Ks_lower_left;];
end
interfacial_stiffness = Ks;

end