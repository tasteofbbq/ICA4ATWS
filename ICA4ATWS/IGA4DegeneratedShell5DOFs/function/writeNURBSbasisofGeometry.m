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

function NURBSbasisofGeometry=writeNURBSbasisofGeometry(Skin)
ieplane=Skin.ieplane;
Gausspoints=Skin.Gausspoints;
empty_cell = cell(size(Gausspoints,1),size(Gausspoints{1,1},1));
NURBSbasisofGeometryRindex = empty_cell;
NURBSbasisofGeometryR = empty_cell;
NURBSbasisofGeometryRs = empty_cell;
NURBSbasisofGeometryRt = empty_cell;
NURBSbasisofGeometrySKL = empty_cell;
NURBSbasisofGeometryv3 = empty_cell;
NURBSbasisofGeometrytheta = empty_cell;
NURBSbasisofGeometryv3s = empty_cell;
NURBSbasisofGeometryv3t = empty_cell;
NURBSbasisofGeometryT = empty_cell;
NURBSbasisofGeometryH = empty_cell;
NURBSbasisofGeometryT_G = empty_cell;
a = size(Gausspoints,1);
b = size(Gausspoints{1,1},1);
parfor i=1:a
    for numGP=1:b
        temp_para = Gausspoints{i,1};
        [R, Rindex]   = nrbbasisfun    ({ temp_para(numGP,1),temp_para(numGP,2) },ieplane);
        [Rs,Rt]       = nrbbasisfunder ({ temp_para(numGP,1),temp_para(numGP,2) },ieplane);
        [SKL] = RatSurfaceDerivs( ieplane,temp_para(numGP,1),temp_para(numGP,2),2);
        [~, ~ ,v3 ,theta] = localCoordinate(SKL);
        [~,~,~,~,v3s,v3t] = dVdst(SKL);
        [T,H] = local2global_voigt( theta );
        [T_G] = local2global_voigt_G( theta );
        NURBSbasisofGeometryRindex{i,numGP} = Rindex;
        NURBSbasisofGeometryR{i,numGP} = R;
        NURBSbasisofGeometryRs{i,numGP} = Rs;
        NURBSbasisofGeometryRt{i,numGP} = Rt;
        NURBSbasisofGeometrySKL{i,numGP} = SKL;
        NURBSbasisofGeometryv3{i,numGP} = v3;
        NURBSbasisofGeometrytheta{i,numGP} = theta;
        NURBSbasisofGeometryv3s{i,numGP} = v3s;
        NURBSbasisofGeometryv3t{i,numGP} = v3t;
        NURBSbasisofGeometryT{i,numGP} = T;
        NURBSbasisofGeometryH{i,numGP} = H;
        NURBSbasisofGeometryT_G{i,numGP} = T_G;
    end
end
NURBSbasisofGeometry.Rindex = NURBSbasisofGeometryRindex;
NURBSbasisofGeometry.R = NURBSbasisofGeometryR;
NURBSbasisofGeometry.Rs = NURBSbasisofGeometryRs;
NURBSbasisofGeometry.Rt = NURBSbasisofGeometryRt;
NURBSbasisofGeometry.SKL = NURBSbasisofGeometrySKL;
NURBSbasisofGeometry.v3 = NURBSbasisofGeometryv3;
NURBSbasisofGeometry.theta = NURBSbasisofGeometrytheta;
NURBSbasisofGeometry.v3s = NURBSbasisofGeometryv3s;
NURBSbasisofGeometry.v3t = NURBSbasisofGeometryv3t;
NURBSbasisofGeometry.T = NURBSbasisofGeometryT;
NURBSbasisofGeometry.H = NURBSbasisofGeometryH;
NURBSbasisofGeometry.T_G = NURBSbasisofGeometryT_G;

end