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
function [u,v,w] = Evaluate_Displacement(para, Skin)

ieplane = Skin.ieplane;

xi   = para(1);
eta  = para(2);
zeta = para(3);

eye_3 = eye(3);

ndof         = ieplane.order(1)*ieplane.order(2);

[ti,~,~]     = Thickness(Skin.anchors,0,Skin);
Displacement = Skin.Displacement;

[R, Rindex]             = nrbbasisfun ({ xi,eta },ieplane);

for iii = 1:ndof
    ID         = Rindex(iii);
    lxn        = Skin.LocalCoor{1,ID}(:,1);
    lyn        = Skin.LocalCoor{1,ID}(:,2);
    V1N(:,iii) = lxn;
    V2N(:,iii) = lyn;
end

elementDof = sort([5*Rindex,5*Rindex-4,5*Rindex-3,5*Rindex-2,5*Rindex-1]);

N = zeros(3,5*ndof);
for iii = 1:ndof
    N(1:3,(5*iii-4):5*iii) = R(iii)*[eye_3 zeta*ti/2*V1N(:,iii) -zeta*ti/2*V2N(:,iii)];
end

disp = N*Displacement(elementDof);

u = disp(1);
v = disp(2);
w = disp(3);

end