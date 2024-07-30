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

function [save_p, dele_p, multi_K_new, f_new, activeDof_new] = control_point_bonding(Total_Skin, coupling_information, multi_K, f, activeDof)
Skin_total = Total_Skin;
save_p = [];
dele_p = [];
for i = 1:length(Total_Skin)
    temp_Skin        = Total_Skin{i};
    Skin_index       = temp_Skin.index;
    Skin{Skin_index} = temp_Skin;
    NoCP = (Skin{Skin_index}.ieplane.number(1,1))*(Skin{Skin_index}.ieplane.number(1,2));
    dof_total(Skin_index) = 5*NoCP;
end
all_dof = 1:sum(dof_total);
for loop_CI = 1:length(coupling_information)
    temp_information     = coupling_information{loop_CI};
    patch_index1         = temp_information(1,1);
    patch_index2         = temp_information(1,2);
    local_surface_index1 = temp_information(2,1);
    local_surface_index2 = temp_information(2,2);
    tolerance            = temp_information(3,1);
    self_tie             = temp_information(4,1);

    if self_tie~=1
        continue
    end

    NoCP1 = (Skin_total{patch_index1}.ieplane.number(1,1))*(Skin_total{patch_index1}.ieplane.number(1,2));
    dof1 = 5*NoCP1;
    NoCP2 = (Skin_total{patch_index2}.ieplane.number(1,1))*(Skin_total{patch_index2}.ieplane.number(1,2));
    dof2 = 5*NoCP2;
    % 1 patch
    ieplane1  = Skin_total{patch_index1}.ieplane;
    % 2 patch
    ieplane2  = Skin_total{patch_index2}.ieplane;
    % 1 patch
    down1  = 1:ieplane1.number(1,1);
    up1    = (ieplane1.number(1,2)-1)*(ieplane1.number(1,1))+1:(ieplane1.number(1,1))*(ieplane1.number(1,2));
    left1  = 1:(ieplane1.number(1,1)):(ieplane1.number(1,2)-1)*(ieplane1.number(1,1))+1;
    right1 = (ieplane1.number(1,1)):(ieplane1.number(1,1)):(ieplane1.number(1,1))*(ieplane1.number(1,2));
    % 2 patch
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

    % There may be issue
    coupling_edge1 = sort([coupling_edge1*5-4 coupling_edge1*5-3 coupling_edge1*5-2 coupling_edge1*5-1 coupling_edge1*5]);
    coupling_edge2 = sort([coupling_edge2*5-4 coupling_edge2*5-3 coupling_edge2*5-2 coupling_edge2*5-1 coupling_edge2*5]);

    sum1 = 0;
    for loop1 = 1:patch_index1-1
        sum1 = sum1+dof_total(loop1);
    end

    sum2 = 0;
    for loop2 = 1:patch_index2-1
        sum2 = sum2+dof_total(loop2);
    end

    save_p = [save_p coupling_edge1+sum1];
    dele_p = [dele_p coupling_edge2+sum2];

end

multi_K(save_p,:) = multi_K(save_p,:) + multi_K(dele_p,:);
f(save_p,:) = f(save_p,:) + f(dele_p,:);
multi_K(:,save_p) = multi_K(:,save_p) + multi_K(:,dele_p);

activeDof = setdiff(activeDof,dele_p);
multi_K_new = multi_K(activeDof,activeDof);
f_new = f(activeDof,:);
activeDof_new = activeDof;

end