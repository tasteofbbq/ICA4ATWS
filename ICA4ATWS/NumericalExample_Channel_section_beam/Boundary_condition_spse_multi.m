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

function [activeDof,fix] = Boundary_condition_spse_multi(Skin)
conpsX = Skin.conpsX;
conpsY = Skin.conpsY;
conpsZ = Skin.conpsZ;
NoCP   = (Skin.ieplane.number(1,1))*(Skin.ieplane.number(1,2));
Gdofs  = 5*NoCP;
Skin_index = Skin.index;

down  = 1:Skin.ieplane.number(1,1);
up    = (Skin.ieplane.number(1,2)-1)*(Skin.ieplane.number(1,1))+1:(Skin.ieplane.number(1,1))*(Skin.ieplane.number(1,2));
left  = 1:(Skin.ieplane.number(1,1)):(Skin.ieplane.number(1,2)-1)*(Skin.ieplane.number(1,1))+1;
right = (Skin.ieplane.number(1,1)):(Skin.ieplane.number(1,1)):(Skin.ieplane.number(1,1))*(Skin.ieplane.number(1,2));


node1=intersect(left,down);
node2=intersect(right,down);
node3=intersect(left,up);
node4=intersect(right,up);
corner = [node1 node2 node3 node4 ]';

fixedNodeU  = [ ];
fixedNodeV  = [ ];
fixedNodeW  = [ ];
fixedNodeTX = [ ];
fixedNodeTY = [ ];

switch Skin_index
    case 1
        fixedNodeU  = [down ];
        fixedNodeV  = [down ];
        fixedNodeW  = [down ];
        fixedNodeTX = [down ];
        fixedNodeTY = [down ];

    case 2
        fixedNodeU  = [down ];
        fixedNodeV  = [down ];
        fixedNodeW  = [down ];
        fixedNodeTX = [down ];
        fixedNodeTY = [down ];

    case 3
        fixedNodeU  = [down ];
        fixedNodeV  = [down ];
        fixedNodeW  = [down ];
        fixedNodeTX = [down ];
        fixedNodeTY = [down ];

end
fix = sort([ 5*fixedNodeU-4 5*fixedNodeV-3 5*fixedNodeW-2 5*fixedNodeTX-1 5*fixedNodeTY]);
fix = unique(fix);
activeDof = setdiff([1:Gdofs]',[fix]);
end