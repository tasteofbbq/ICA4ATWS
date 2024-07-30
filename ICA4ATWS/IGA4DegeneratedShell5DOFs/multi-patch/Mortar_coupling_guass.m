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

function [multi_K,multi_K_b] = Mortar_coupling_guass(Total_Skin, information, coupling_points)
multi_K = [];
multi_KC_assmble = [];
for i = 1:length(Total_Skin)
    temp_Skin        = Total_Skin{i};
    Skin_index       = temp_Skin.index;
    Skin{Skin_index} = temp_Skin;
    [t,~,~]          = Thickness(temp_Skin.anchors,0,temp_Skin);
    Skin{Skin_index}.ti = t;
    K{Skin_index}    = K_for_Skin(temp_Skin, temp_Skin.DL_Shell, t,0);
    NoCP = (Skin{Skin_index}.ieplane.number(1,1))*(Skin{Skin_index}.ieplane.number(1,2));
    dof_total(Skin_index) = 5*NoCP;
end
for Skin_index = 1:length(Total_Skin)
    multi_K = blkdiag(multi_K, K{Skin_index});
end
multi_K_b = multi_K;

if(~isempty(information))
    for i = 1:length(information)
        temp_information     = information{i};
        temp_coupling_points = coupling_points{i};
        self_tie = temp_information(4,1);
        switch self_tie
            case 0
                [interfacial_stiffness] = interfacial_stiffness_matrices4Mortar_guass(Skin, temp_information,temp_coupling_points);
            case 1
                continue;
            case 2
                [interfacial_stiffness] = interfacial_stiffness_matrices4Mortar_self_tie(Skin, temp_information);
        end
        patch_index1         = temp_information(1,1);
        patch_index2         = temp_information(1,2);
        dof1 = dof_total(patch_index1);
        dof2 = dof_total(patch_index2);
        if(patch_index1>patch_index2)
            t    = dof1;
            dof1 = dof2;
            dof2 = t;
            t            = patch_index1;
            patch_index1 = patch_index2;
            patch_index2 = t;
        end
        sum1 = 0;
        sum2 = 0;
        for loop1 = 1:patch_index1-1
            sum1 = sum1+dof_total(loop1);
        end
        for loop2 = 1:patch_index2-1
            sum2 = sum2+dof_total(loop2);
        end
        patchDof1 = sum1+1:sum1+dof1;
        patchDof2 = sum2+1:sum2+dof2;
        patchDof = [patchDof1 patchDof2];

        non_zero_rows = any(interfacial_stiffness, 2);
        interfacial_stiffness = interfacial_stiffness(non_zero_rows, :);

        multi_KC = sparse(size(interfacial_stiffness,1),size(multi_K_b,1));
        if(patch_index1~=patch_index2)
            multi_KC(:,patchDof) = interfacial_stiffness;
        end
        if(patch_index1==patch_index2)
            multi_KC(:,patchDof1) = interfacial_stiffness(:,1:dof1)+interfacial_stiffness(:,1+dof1:dof1+dof1);
        end
        multi_KC_assmble = [multi_KC_assmble;multi_KC];
        fprintf('No.%d coupling(Mortar) has been done!\n', i)
    end

    multi_K = [multi_K          multi_KC_assmble';
               multi_KC_assmble sparse(size(multi_KC_assmble,1),size(multi_KC_assmble,1))];
end

end