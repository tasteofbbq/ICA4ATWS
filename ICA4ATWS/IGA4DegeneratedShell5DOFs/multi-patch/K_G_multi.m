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

function multi_K_G = K_G_multi(Total_Skin)

multi_K_G = [];
for i = 1:length(Total_Skin)
    temp_Skin        = Total_Skin{i};
    Skin_index       = temp_Skin.index;
    [t,~,~]          = Thickness(temp_Skin.anchors,0,temp_Skin);
    K_G{Skin_index}    = K_G_for_Skin(temp_Skin, temp_Skin.DL_Shell, t,0,temp_Skin.Displacement);
end
for Skin_index = 1:length(Total_Skin)
    multi_K_G = blkdiag(multi_K_G, K_G{Skin_index});
end
end