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

function [activeDof] = Boundary_condition_multi(Total_Skin,coupling_information)
if(isempty(coupling_information)&&length(Total_Skin)==1)
    Skin = Total_Skin{end};
    [activeDof,fix] = Boundary_condition_spse_single(Skin);
else
    activeDof = [];
    for i = 1:length(Total_Skin)
        temp_Skin        = Total_Skin{i};
        Skin_index       = temp_Skin.index;
        Skin{Skin_index} = temp_Skin;
        NoCP = (Skin{Skin_index}.ieplane.number(1,1))*(Skin{Skin_index}.ieplane.number(1,2));
        dof_total(Skin_index) = 5*NoCP;
    end

    for Skin_index = 1:length(Skin)
        [activeDof_temp,fix_temp] = Boundary_condition_spse_multi(Skin{Skin_index});
        sum = 0;
        for loop = 1:Skin_index-1
            sum = sum+dof_total(loop);
        end
        activeDof = [activeDof activeDof_temp'+sum];
    end
end

end