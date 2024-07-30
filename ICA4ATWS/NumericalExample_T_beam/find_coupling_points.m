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

function coupling_points = find_coupling_points(Total_Skin, coupling_information)
for i = 1:length(Total_Skin)
    temp_Skin        = Total_Skin{i};
    Skin_index       = temp_Skin.index;
    Skin_total{Skin_index} = temp_Skin;
end

coupling_points = cell(1,length(coupling_information));
for loop_coup_info = 1:length(coupling_information)
    temp_information = coupling_information{loop_coup_info};
    patch_index1         = temp_information(1,1);
    patch_index2         = temp_information(1,2);
    local_surface_index1 = temp_information(2,1);
    local_surface_index2 = temp_information(2,2);

    ieplane1  = Skin_total{patch_index1}.ieplane;
    conpsX1   = Skin_total{patch_index1}.conpsX;
    conpsY1   = Skin_total{patch_index1}.conpsY;
    conpsZ1   = Skin_total{patch_index1}.conpsZ;

    num_GS_Thi = 2;
    num_GS = 2;
    switch local_surface_index1
        case 1
            u = unique(ieplane1.knots{1});
            v = 0;
            [wGS, wW] = Gausspoint(-1,1,num_GS_Thi);
            coupling_points_guass = cell(size(u,2)-1,num_GS,num_GS_Thi);
            for ii = 1:size(u,2)-1
                [s, w] = Gausspoint(u(ii),u(ii+1),num_GS);
                for numGP = 1:num_GS
                    [R1, Rindex1]   = nrbbasisfun ({ s(numGP),v },ieplane1);
                    [SKL1] = RatSurfaceDerivs( ieplane1,s(numGP),v,2);
                    [~, ~ ,v31 ,theta1] = localCoordinate(SKL1);
                    [t1,ts1,tt1] = Thickness([s(numGP),v],0,Skin_total{patch_index1});
                    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];
                    for qw = 1:num_GS_Thi
                        gs1 = wGS(qw);
                        X = R1*Cps1+gs1*v31'*(t1/2);
                        [found_para, distance] = find_Matlab(X,Skin_total{patch_index2},local_surface_index2);
                        w1 = w(numGP);
                        w2 = wW(qw,1);
                        para = [s(numGP) v gs1];
                        coupling_points_guass{ii,numGP,qw} = [para;found_para;w1 w2 distance];
                    end
                end
            end
            coupling_points_guass = reshape(coupling_points_guass,1,(size(u,2)-1)*num_GS*num_GS_Thi);

        case 2
            u = unique(ieplane1.knots{1});
            v = 1;
            [wGS, wW] = Gausspoint(-1,1,num_GS_Thi);
            coupling_points_guass = cell(size(u,2)-1,num_GS,num_GS_Thi);
            for ii = 1:size(u,2)-1
                [s, w] = Gausspoint(u(ii),u(ii+1),num_GS);
                for numGP = 1:num_GS
                    [R1, Rindex1]   = nrbbasisfun ({ s(numGP),v },ieplane1);
                    [SKL1] = RatSurfaceDerivs( ieplane1,s(numGP),v,2);
                    [~, ~ ,v31 ,theta1] = localCoordinate(SKL1);
                    [t1,ts1,tt1] = Thickness([s(numGP),v],Skin_total{patch_index1}.theta,Skin_total{patch_index1});
                    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];
                    for qw = 1:num_GS_Thi
                        gs1 = wGS(qw);
                        X = R1*Cps1+gs1*v31'*(t1/2);
                        [found_para, distance] = find_Matlab(X,Skin_total{patch_index2},local_surface_index2);
                        w1 = w(numGP);
                        w2 = wW(qw,1);
                        para = [s(numGP) v gs1];
                        coupling_points_guass{ii,numGP,qw} = [para;found_para;w1 w2 distance];
                    end
                end
            end
            coupling_points_guass = reshape(coupling_points_guass,1,(size(u,2)-1)*num_GS*num_GS_Thi);

        case 3

        case 4
            u = 1;
            v = unique(ieplane1.knots{2});
            [wGS, wW] = Gausspoint(-1,1,num_GS_Thi);
            coupling_points_guass = cell(size(v,2)-1,num_GS,num_GS_Thi);
            for ii = 1:size(v,2)-1
                [s, w] = Gausspoint(v(ii),v(ii+1),num_GS);
                for numGP = 1:num_GS
                    [R1, Rindex1]   = nrbbasisfun ({ u,s(numGP) },ieplane1);
                    [SKL1] = RatSurfaceDerivs( ieplane1,u,s(numGP),2);
                    [~, ~ ,v31 ,theta1] = localCoordinate(SKL1);
                    [t1,ts1,tt1] = Thickness([u,s(numGP)],Skin_total{patch_index1}.theta,Skin_total{patch_index1});
                    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];
                    for qw = 1:num_GS_Thi
                        gs1 = wGS(qw);
                        X = R1*Cps1+gs1*v31'*(t1/2);
                        [found_para, distance] = find_Matlab(X,Skin_total{patch_index2},local_surface_index2);
                        w1 = w(numGP);
                        w2 = wW(qw,1);
                        para = [u,s(numGP) gs1];
                        coupling_points_guass{ii,numGP,qw} = [para;found_para;w1 w2 distance];
                    end
                end
            end
            coupling_points_guass = reshape(coupling_points_guass,1,(size(v,2)-1)*num_GS*num_GS_Thi);

        case 5

        case 6
            u = 0;
            v = unique(ieplane1.knots{2});
            [wGS, wW] = Gausspoint(-1,1,num_GS_Thi);
            coupling_points_guass = cell(size(v,2)-1,num_GS,num_GS_Thi);
            for ii = 1:size(v,2)-1
                [s, w] = Gausspoint(v(ii),v(ii+1),num_GS);
                for numGP = 1:num_GS
                    [R1, Rindex1]   = nrbbasisfun ({ u,s(numGP) },ieplane1);
                    [SKL1] = RatSurfaceDerivs( ieplane1,u,s(numGP),2);
                    [~, ~ ,v31 ,theta1] = localCoordinate(SKL1);
                    [t1,ts1,tt1] = Thickness([u,s(numGP)],Skin_total{patch_index1}.theta,Skin_total{patch_index1});
                    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];
                    for qw = 1:num_GS_Thi
                        gs1 = wGS(qw);
                        X = R1*Cps1+gs1*v31'*(t1/2);
                        [found_para, distance] = find_Matlab(X,Skin_total{patch_index2},local_surface_index2);
                        w1 = w(numGP);
                        w2 = wW(qw,1);
                        para = [u,s(numGP) gs1];
                        coupling_points_guass{ii,numGP,qw} = [para;found_para;w1 w2 distance];
                    end
                end
            end
            coupling_points_guass = reshape(coupling_points_guass,1,(size(v,2)-1)*num_GS*num_GS_Thi);
    end
    coupling_points{loop_coup_info} = coupling_points_guass;
    fprintf('No.%d coupling information has been searched!\n', loop_coup_info)
end
end