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

function [interfacial_stiffness] = interfacial_stiffness_matrices4Penalty_guass(Skin_total, temp_information,temp_coupling_points,alpha)
patch_index1         = temp_information(1,1);
patch_index2         = temp_information(1,2);
local_surface_index1 = temp_information(2,1);
local_surface_index2 = temp_information(2,2);
tolerance            = temp_information(3,1);

NoCP1 = (Skin_total{patch_index1}.ieplane.number(1,1))*(Skin_total{patch_index1}.ieplane.number(1,2));
dof1 = 5*NoCP1;
NoCP2 = (Skin_total{patch_index2}.ieplane.number(1,1))*(Skin_total{patch_index2}.ieplane.number(1,2));
dof2 = 5*NoCP2;

Ks_upper_left  = sparse(dof1,dof1);
Ks_upper_right = sparse(dof1,dof2);
Ks_lower_left  = sparse(dof2,dof1);
Ks_lower_right = sparse(dof2,dof2);

ieplane1  = Skin_total{patch_index1}.ieplane;
conpsX1   = Skin_total{patch_index1}.conpsX;
conpsY1   = Skin_total{patch_index1}.conpsY;
conpsZ1   = Skin_total{patch_index1}.conpsZ;
ndof1     = ieplane1.order(1)*ieplane1.order(2);
ti1 = Skin_total{patch_index1}.ti;

ieplane2  = Skin_total{patch_index2}.ieplane;
conpsX2   = Skin_total{patch_index2}.conpsX;
conpsY2   = Skin_total{patch_index2}.conpsY;
conpsZ2   = Skin_total{patch_index2}.conpsZ;
ndof2     = ieplane2.order(1)*ieplane2.order(2);
ti2 = Skin_total{patch_index2}.ti;

V1N1_all = zeros(3,NoCP1);
V2N1_all = zeros(3,NoCP1);
V1N2_all = zeros(3,NoCP2);
V2N2_all = zeros(3,NoCP2);
Skin_total_patch_index1_LocalCoor = Skin_total{patch_index1}.LocalCoor;
Skin_total_patch_index2_LocalCoor = Skin_total{patch_index2}.LocalCoor;

for i = 1:NoCP1
    lxn = Skin_total_patch_index1_LocalCoor{1,i}(:,1);
    lyn = Skin_total_patch_index1_LocalCoor{1,i}(:,2);
    V1N1_all(:,i) = lxn;
    V2N1_all(:,i) = lyn;
end

for i = 1:NoCP2
    lxn = Skin_total_patch_index2_LocalCoor{1,i}(:,1);
    lyn = Skin_total_patch_index2_LocalCoor{1,i}(:,2);
    V1N2_all(:,i) = lxn;
    V2N2_all(:,i) = lyn;
end

elementDof1_son = cell(length(temp_coupling_points));
elementDof2_son = cell(length(temp_coupling_points));
Ks_upper_left_son  = cell(length(temp_coupling_points));
Ks_upper_right_son = cell(length(temp_coupling_points));
Ks_lower_left_son  = cell(length(temp_coupling_points));
Ks_lower_right_son = cell(length(temp_coupling_points));

for loop_guass = 1:length(temp_coupling_points)

    fac1 = 0;

    para = temp_coupling_points{loop_guass}(1,:);
    found_para = temp_coupling_points{loop_guass}(2,:);
    w1 = temp_coupling_points{loop_guass}(3,1);
    w2 = temp_coupling_points{loop_guass}(3,2);
    distance = temp_coupling_points{loop_guass}(3,3);
    if(distance>tolerance)
        elementDof1_son{loop_guass} = [];
        elementDof2_son{loop_guass} = [];
        Ks_upper_left_son{loop_guass}  = 0;
        Ks_upper_right_son{loop_guass} = 0;
        Ks_lower_left_son{loop_guass}  = 0;
        Ks_lower_right_son{loop_guass} = 0;
        continue;
    end

    [R1, Rindex1]   = nrbbasisfun ({ para(1),para(2) },ieplane1);
    [Rs1,Rt1] = nrbbasisfunder ({ para(1),para(2) },ieplane1);
    [SKL1] = RatSurfaceDerivs( ieplane1,para(1),para(2),2);
    [~, ~ ,v31 ,theta1] = localCoordinate(SKL1);
    [~,~,~,~,v3s1,v3t1] = dVdst(SKL1);

    [t1,ts1,tt1] = Thickness([para(1),para(2)],0,Skin_total{patch_index1});
    V1N1 = V1N1_all(:,Rindex1);
    V2N1 = V2N1_all(:,Rindex1);
    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];

    JS11 = Rs1*Cps1 + para(3)*v3s1'*(t1/2) + para(3)*v31'*(ts1/2);
    JS12 = Rt1*Cps1 + para(3)*v3t1'*(t1/2) + para(3)*v31'*(tt1/2);
    JS13 = v31'*(t1/2);
    JS1  = [JS11;JS12;JS13];
    switch local_surface_index1
        case 1
            fac1 = sqrt((JS1(1,2)*JS1(3,3)+JS1(3,2)*JS1(1,3))^2+(JS1(1,3)*JS1(3,1)+JS1(3,3)*JS1(1,1))^2+(JS1(1,1)*JS1(3,2)+JS1(3,1)*JS1(1,2))^2)*w1*w2;
        case 2
            fac1 = sqrt((JS1(1,2)*JS1(3,3)+JS1(3,2)*JS1(1,3))^2+(JS1(1,3)*JS1(3,1)+JS1(3,3)*JS1(1,1))^2+(JS1(1,1)*JS1(3,2)+JS1(3,1)*JS1(1,2))^2)*w1*w2;
        case 3
        case 4
            fac1 = sqrt((JS1(2,2)*JS1(3,3)+JS1(3,2)*JS1(2,3))^2+(JS1(2,3)*JS1(3,1)+JS1(3,3)*JS1(2,1))^2+(JS1(2,1)*JS1(3,2)+JS1(3,1)*JS1(2,2))^2)*w1*w2;
        case 5
        case 6
            fac1 = sqrt((JS1(2,2)*JS1(3,3)+JS1(3,2)*JS1(2,3))^2+(JS1(2,3)*JS1(3,1)+JS1(3,3)*JS1(2,1))^2+(JS1(2,1)*JS1(3,2)+JS1(3,1)*JS1(2,2))^2)*w1*w2;
    end

    N1 = zeros(3,5*ndof1);
    for i = 1:ndof1
        N1(1:3,(5*i-4):5*i)=R1(i)*[eye(3) para(3)*ti1/2*V1N1(:,i) -para(3)*ti1/2*V2N1(:,i)];
    end

    [R2, Rindex2]   = nrbbasisfun ({ found_para(1),found_para(2) },ieplane2);

    V1N2 = V1N2_all(:,Rindex2);
    V2N2 = V2N2_all(:,Rindex2);

    N2 = zeros(3,5*ndof2);
    gs2 = found_para(3);
    for i = 1:ndof2
        N2(1:3,(5*i-4):5*i)=R2(i)*[eye(3) gs2*ti2/2*V1N2(:,i) -gs2*ti2/2*V2N2(:,i)];
    end

    elementDof1_son{loop_guass} = sort([5*Rindex1,5*Rindex1-4,5*Rindex1-3,5*Rindex1-2,5*Rindex1-1]);
    elementDof2_son{loop_guass} = sort([5*Rindex2,5*Rindex2-4,5*Rindex2-3,5*Rindex2-2,5*Rindex2-1]);
    Ks_upper_left_son{loop_guass}  = transpose(N1)*N1*fac1;
    Ks_upper_right_son{loop_guass} = transpose(N1)*N2*fac1;
    Ks_lower_left_son{loop_guass}  = transpose(N2)*N1*fac1;
    Ks_lower_right_son{loop_guass} = transpose(N2)*N2*fac1;
end

for loop_guass = 1:length(temp_coupling_points)
    elementDof1 = elementDof1_son{loop_guass};
    elementDof2 = elementDof2_son{loop_guass};
    Ks_upper_left(elementDof1,elementDof1)  = Ks_upper_left(elementDof1,elementDof1)+Ks_upper_left_son{loop_guass};
    Ks_upper_right(elementDof1,elementDof2) = Ks_upper_right(elementDof1,elementDof2)+Ks_upper_right_son{loop_guass};
    Ks_lower_left(elementDof2,elementDof1)  = Ks_lower_left(elementDof2,elementDof1)+Ks_lower_left_son{loop_guass};
    Ks_lower_right(elementDof2,elementDof2) = Ks_lower_right(elementDof2,elementDof2)+Ks_lower_right_son{loop_guass};
end
Ks_upper_left  = sparse(Ks_upper_left);
Ks_upper_right = sparse(Ks_upper_right);
Ks_lower_left  = sparse(Ks_lower_left);
Ks_lower_right = sparse(Ks_lower_right);

if(patch_index1<=patch_index2)
    Ks = [alpha*Ks_upper_left -alpha*Ks_upper_right;
        -alpha*Ks_lower_left  alpha*Ks_lower_right];
else
    Ks = [alpha*Ks_lower_right -alpha*Ks_lower_left;
        -alpha*Ks_upper_right  alpha*Ks_upper_left];
end
interfacial_stiffness = Ks;
end