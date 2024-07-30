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

function [Skin] = skin_operation_TTO_multipatch(nurbsSurf,NuandNv,E,miu,Skinindex)
F_geometry = 0;
F_parameter = 0;

Nu = NuandNv(1);
Nv = NuandNv(2);

Skin     = skin_operation(nurbsSurf,Nu + 1,Nv + 1,F_geometry,F_parameter);
Skin.e1 = E;
Skin.miu = miu;
Skin.DL_Shell = Elastic_Shell(Skin.e1,Skin.miu);
Skin.DL_Solid = Elastic_Solid(Skin.e1,Skin.miu);
anchors  = Skin.CPw;anchors=cell2mat(anchors);
Skin.anchors = reshape(anchors,Skin.Ncops,2);

Skin.Gausspoints          = writeGauss(Skin);
Skin.NURBSbasisofGeometry = writeNURBSbasisofGeometry(Skin);
Skin.index = Skinindex;
Skin.fi = [];

NoCP = (Skin.ieplane.number(1,1))*(Skin.ieplane.number(1,2));
dof_total = 5*NoCP;
Skin.Displacement = zeros(dof_total,1);

end