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

function [ stiffness ] = K_Element_Skin_STO( Skin,GP,Weight,DL,t,thetaofthickness,ii)
ieplane = Skin.ieplane;
NURBSbasisofGeometry=Skin.NURBSbasisofGeometry;

% Homogeneous coordinates are transformed into Cartesian coordinates
conpsX = Skin.conpsX;
conpsY = Skin.conpsY;
conpsZ = Skin.conpsZ;

% Number of control points in a element
ndof = ieplane.order(1)*ieplane.order(2);
stiffness = zeros(5*ndof,5*ndof);

Rindex = NURBSbasisofGeometry.Rindex{ii,1};
ti=t; % Thickness

Cps = [conpsX(Rindex),conpsY(Rindex),conpsZ(Rindex)];
%%
% Calculate the gauss points and its weight
[gs3,w3] = Gausspoint(-1,1,2);
for numGP = 1:size(GP,1)       % Eta circulation
    % Create shapefunctions and its derivative of NURBS and serial number of control points in a element
    Rindex = NURBSbasisofGeometry.Rindex{ii,numGP};
    R = NURBSbasisofGeometry.R{ii,numGP};
    Rs = NURBSbasisofGeometry.Rs{ii,numGP};
    Rt = NURBSbasisofGeometry.Rt{ii,numGP};
    v3 = NURBSbasisofGeometry.v3{ii,numGP};
    v3s = NURBSbasisofGeometry.v3s{ii,numGP};
    v3t = NURBSbasisofGeometry.v3t{ii,numGP};
    T = NURBSbasisofGeometry.T{ii,numGP};
    H = NURBSbasisofGeometry.H{ii,numGP};

    t  = ti;
    ts = 0;
    tt = 0;
    rou = 1;

    for i = 1:ndof
        ID = Rindex(i);
        lxn = Skin.LocalCoor{1,ID}(:,1);
        lyn = Skin.LocalCoor{1,ID}(:,2);
        V1N(:,i) = lxn;
        V2N(:,i) = lyn;
    end

    for qw = 1:size(gs3,1)
        % Calculate the overall coordinate of zeta via the local
        gs33 = gs3(qw);
        % Calculate the Jacob matrix
        JS(1,:) = Rs*Cps + gs33*v3s'*(t/2) + gs33*v3'*(ts/2);
        JS(2,:) = Rt*Cps + gs33*v3t'*(t/2) + gs33*v3'*(tt/2);
        JS(3,:) = v3'*(t/2);
        invJS = inv(JS);
        F = blkdiag(invJS,invJS,invJS);
        fac = det(JS)*Weight(numGP)*w3(qw,1);
        G = zeros(9,5*ndof);
        for i = 1:ndof
            Ri    = R(i);
            Rsi   = Rs(i);
            Rti   = Rt(i);
            V1N1i = V1N(1,i);
            V1N2i = V1N(2,i);
            V1N3i = V1N(3,i);
            V2N1i = V2N(1,i);
            V2N2i = V2N(2,i);
            V2N3i = V2N(3,i);
            tii   = ti;
            G(1:9,(5*i-4):5*i)=[ Rsi   0    0   gs33*(V1N1i*Rsi)*tii/2   -gs33*(V2N1i*Rsi)*tii/2 
                                 Rti   0    0   gs33*(V1N1i*Rti)*tii/2   -gs33*(V2N1i*Rti)*tii/2
                                  0      0    0               V1N1i*Ri*tii/2   -V2N1i*Ri*tii/2
                                  0    Rsi  0   gs33*(V1N2i*Rsi)*tii/2   -gs33*(V2N2i*Rsi)*tii/2 
                                  0    Rti  0   gs33*(V1N2i*Rti)*tii/2   -gs33*(V2N2i*Rti)*tii/2
                                  0      0    0               V1N2i*Ri*tii/2   -V2N2i*Ri*tii/2
                                  0      0  Rsi gs33*(V1N3i*Rsi)*tii/2   -gs33*(V2N3i*Rsi)*tii/2 
                                  0      0  Rti gs33*(V1N3i*Rti)*tii/2   -gs33*(V2N3i*Rti)*tii/2
                                  0      0    0               V1N3i*Ri*tii/2   -V2N3i*Ri*tii/2];
        end
        B = T*H*F*G;
        stiffness  =  stiffness +  B'*DL*B*fac*rou;
    end
end
end