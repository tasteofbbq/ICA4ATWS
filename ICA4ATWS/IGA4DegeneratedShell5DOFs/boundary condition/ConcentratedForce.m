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

function forcevector = ConcentratedForce(Skin,para,f_global,forcevector)
ieplane=Skin.ieplane;
[R, Rindex] = nrbbasisfun ({para(1),para(2)},ieplane);
for i = 1:size(R,2)
    temp_f=[R(i) 0 0;0 R(i) 0;0 0 R(i)]*f_global;
    active=[Rindex(i)*5-4,Rindex(i)*5-3,Rindex(i)*5-2];
    forcevector(active) = forcevector(active) + temp_f;
end
end