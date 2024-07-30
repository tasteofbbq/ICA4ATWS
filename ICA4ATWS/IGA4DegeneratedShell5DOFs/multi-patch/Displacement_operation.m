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

function Total_Skin = Displacement_operation(Total_Skin, Displacement)
for i = 1:length(Total_Skin)
    temp_Skin        = Total_Skin{i};
    Skin_index       = temp_Skin.index;
    NoCP = (temp_Skin.ieplane.number(1,1))*(temp_Skin.ieplane.number(1,2));
    dof_total(Skin_index) = 5*NoCP;
end

for i = 1:length(Total_Skin)
    Skin_index = Total_Skin{i}.index;
    sum = 0;
    for loop = 1:Skin_index-1
        sum = sum+dof_total(loop);
    end
    dofs = sum + 1:sum + dof_total(Skin_index);
    Total_Skin{i}.Displacement = Displacement(dofs);
end
end