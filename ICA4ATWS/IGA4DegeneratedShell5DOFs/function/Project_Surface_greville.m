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

function CPw =Project_Surface_greville(Parameter)
% anchors is formed as (number of control points)*2 matrix
ieplane = Parameter.ieplane;
% get the knot vector
U=ieplane.knots{1};V=ieplane.knots{2};
% U=unique(U);V=unique(V);
nopu=ieplane.number(1);nopv=ieplane.number(2);
p1=ieplane.order(1);p2=ieplane.order(2);
% p1=ieplane.order(1)-1;p2=ieplane.order(2)-1;
anchorsU=zeros(nopu);anchorsV=zeros(nopv);
for i=1:nopu
    for j=1:p1-1
        anchorsU(i)=anchorsU(i)+U(i+j);
    end
end
anchorsU=anchorsU./(p1-1);
for i=1:nopv
    for j=1:p2-1
        anchorsV(i)=anchorsV(i)+V(i+j);
    end
end
anchorsV=anchorsV./(p2-1);
for j=1:nopv
    for i=1:nopu
        anchors(i+(j-1)*nopu,1)=anchorsU(i);
        anchors(i+(j-1)*nopu,2)=anchorsV(j);
    end
end
for loop = 1:size(anchors,1)
    CPw{1,loop} = anchors(loop,:);
end
end