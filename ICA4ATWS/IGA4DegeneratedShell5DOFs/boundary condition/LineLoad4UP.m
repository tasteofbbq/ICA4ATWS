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

function forcevector = LineLoad4UP(Skin,f_global,forcevector)
ieplane=Skin.ieplane;
%% v = 1,w = 0,u->[0,1]
u = unique(ieplane.knots{1});
v = 1;
w = 0;
% Homogeneous coordinates are transformed into Cartesian coordinates
conpsX = Skin.conpsX;
conpsY = Skin.conpsY;
conpsZ = Skin.conpsZ;
% Number of control points in a element
ndof = ieplane.order(1)*ieplane.order(2);
for ii = 1:size(u,2)-1
    [s, w] = Gausspoint(u(ii),u(ii+1),2);
    for numGP = 1:size(s,1)
        [R, Rindex]   = nrbbasisfun ({ s(numGP),v },ieplane);
        [Rs,Rt] = nrbbasisfunder ({ s(numGP),v },ieplane);
        [SKL] = RatSurfaceDerivs( ieplane,s(numGP),v,2);
        [~, ~ ,v3 ,theta] = localCoordinate(SKL);
        [~,~,~,~,v3s,v3t] = dVdst(SKL);
        [T,H] = local2global_voigt( theta );
        for i = 1:ndof
            ID = Rindex(i);
            lxn = Skin.LocalCoor{1,ID}(:,1);
            lyn = Skin.LocalCoor{1,ID}(:,2);
            V1N(:,i) = lxn;
            V2N(:,i) = lyn;
        end
        % Control point coordinate in a element
        Cps = [conpsX(Rindex),conpsY(Rindex),conpsZ(Rindex)];
        JS(1,:) = Rs*Cps;
        N = zeros(3,5*ndof);
        for i = 1:ndof
            N(1:3,(5*i-4):5*i)=[ R(i)    0      0    0 0 
                                  0     R(i)    0    0 0
                                  0      0    R(i)   0 0];
        end
    temp_f=N'*f_global*norm(JS)*w(numGP);
    active=sort([Rindex*5-4,Rindex*5-3,Rindex*5-2,Rindex*5-1,Rindex*5]);
    forcevector(active) = forcevector(active) + temp_f;
    end
end
end