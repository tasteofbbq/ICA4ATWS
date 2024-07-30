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

function Gausspoints=writeGauss(Skin)
elementNode = Skin.elementNode;
noElement = Skin.noElement;
elementType = zeros(1,noElement);

fullEle = find(elementType == 0);
nofullEle = size(fullEle,2);
for i = 1:nofullEle
    ind = fullEle(i);
    s1 = elementNode.vertex{1, ind}(1,1);
    s2 = elementNode.vertex{1, ind}(2,1);
    t1 = elementNode.vertex{1, ind}(2,2);
    t2 = elementNode.vertex{1, ind}(3,2);
    
    [ GaussPoint,Weight_st ] = GaussQ9( s1, s2, t1, t2 );
    Gausspoints{ind,1} = GaussPoint;
    Gausspoints{ind,2} = Weight_st;
end
end