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

function [interfacial_stiffness] = interfacial_stiffness_matrices4Nitsche_guass(Skin_total, temp_information,temp_coupling_points,gamma,alpha)
patch_index1         = temp_information(1,1);
patch_index2         = temp_information(1,2);
local_surface_index1 = temp_information(2,1);
local_surface_index2 = temp_information(2,2);
tolerance            = temp_information(3,1);

NoCP1 = (Skin_total{patch_index1}.ieplane.number(1,1))*(Skin_total{patch_index1}.ieplane.number(1,2));
dof1 = 5*NoCP1;
NoCP2 = (Skin_total{patch_index2}.ieplane.number(1,1))*(Skin_total{patch_index2}.ieplane.number(1,2));
dof2 = 5*NoCP2;

Kn_upper_left  = sparse(dof1,dof1);
Kn_upper_right = sparse(dof1,dof2);
Kn_lower_left  = sparse(dof2,dof1);
Kn_lower_right = sparse(dof2,dof2);
Ks_upper_left  = sparse(dof1,dof1);
Ks_upper_right = sparse(dof1,dof2);
Ks_lower_left  = sparse(dof2,dof1);
Ks_lower_right = sparse(dof2,dof2);

ieplane1  = Skin_total{patch_index1}.ieplane;
conpsX1   = Skin_total{patch_index1}.conpsX;
conpsY1   = Skin_total{patch_index1}.conpsY;
conpsZ1   = Skin_total{patch_index1}.conpsZ;
ndof1     = ieplane1.order(1)*ieplane1.order(2);

DL_Solid1 = Skin_total{patch_index1}.DL_Solid;
ti1 = Skin_total{patch_index1}.ti;

ieplane2  = Skin_total{patch_index2}.ieplane;
conpsX2   = Skin_total{patch_index2}.conpsX;
conpsY2   = Skin_total{patch_index2}.conpsY;
conpsZ2   = Skin_total{patch_index2}.conpsZ;
ndof2     = ieplane2.order(1)*ieplane2.order(2);

DL_Solid2 = Skin_total{patch_index2}.DL_Solid;
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
Kn_upper_left_son  = cell(length(temp_coupling_points));
Kn_upper_right_son = cell(length(temp_coupling_points));
Kn_lower_left_son  = cell(length(temp_coupling_points));
Kn_lower_right_son = cell(length(temp_coupling_points));
Ks_upper_left_son  = cell(length(temp_coupling_points));
Ks_upper_right_son = cell(length(temp_coupling_points));
Ks_lower_left_son  = cell(length(temp_coupling_points));
Ks_lower_right_son = cell(length(temp_coupling_points));

for loop_guass = 1:length(temp_coupling_points)

    fac1 = 0;
    n1 = 0;

    para = temp_coupling_points{loop_guass}(1,:);
    found_para = temp_coupling_points{loop_guass}(2,:);
    w1 = temp_coupling_points{loop_guass}(3,1);
    w2 = temp_coupling_points{loop_guass}(3,2);
    distance = temp_coupling_points{loop_guass}(3,3);
    if(distance>tolerance)
        elementDof1_son{loop_guass} = [];
        elementDof2_son{loop_guass} = [];
        Kn_upper_left_son{loop_guass}  = 0;
        Kn_upper_right_son{loop_guass} = 0;
        Kn_lower_left_son{loop_guass}  = 0;
        Kn_lower_right_son{loop_guass} = 0;
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
    [T1,H1] = local2global_voigt( theta1 );
    [t1,ts1,tt1] = Thickness([para(1),para(2)],0,Skin_total{patch_index1});
    V1N1 = V1N1_all(:,Rindex1);
    V2N1 = V2N1_all(:,Rindex1);
    Cps1 = [conpsX1(Rindex1),conpsY1(Rindex1),conpsZ1(Rindex1)];

    rou1 = 1;

    C1 = rou1*DL_Solid1;
    JS11 = Rs1*Cps1 + para(3)*v3s1'*(t1/2) + para(3)*v31'*(ts1/2);
    JS12 = Rt1*Cps1 + para(3)*v3t1'*(t1/2) + para(3)*v31'*(tt1/2);
    JS13 = v31'*(t1/2);
    JS1  = [JS11;JS12;JS13];
    switch local_surface_index1
        case 1
            fac1 = sqrt((JS1(1,2)*JS1(3,3)+JS1(3,2)*JS1(1,3))^2+(JS1(1,3)*JS1(3,1)+JS1(3,3)*JS1(1,1))^2+(JS1(1,1)*JS1(3,2)+JS1(3,1)*JS1(1,2))^2)*w1*w2;
            dsd1 = [SKL1(2,1,1);SKL1(2,1,2);SKL1(2,1,3)];
            n1 = cross(dsd1,v31);
        case 2
            fac1 = sqrt((JS1(1,2)*JS1(3,3)+JS1(3,2)*JS1(1,3))^2+(JS1(1,3)*JS1(3,1)+JS1(3,3)*JS1(1,1))^2+(JS1(1,1)*JS1(3,2)+JS1(3,1)*JS1(1,2))^2)*w1*w2;
            dsd1 = [SKL1(2,1,1);SKL1(2,1,2);SKL1(2,1,3)];
            n1 = -cross(dsd1,v31);
        case 3
            n1 = v31;
        case 4
            fac1 = sqrt((JS1(2,2)*JS1(3,3)+JS1(3,2)*JS1(2,3))^2+(JS1(2,3)*JS1(3,1)+JS1(3,3)*JS1(2,1))^2+(JS1(2,1)*JS1(3,2)+JS1(3,1)*JS1(2,2))^2)*w1*w2;
            dsd2 = [SKL1(1,2,1);SKL1(1,2,2);SKL1(1,2,3)];
            n1 = cross(dsd2,v31);
        case 5
            n1 = -v31;
        case 6
            fac1 = sqrt((JS1(2,2)*JS1(3,3)+JS1(3,2)*JS1(2,3))^2+(JS1(2,3)*JS1(3,1)+JS1(3,3)*JS1(2,1))^2+(JS1(2,1)*JS1(3,2)+JS1(3,1)*JS1(2,2))^2)*w1*w2;
            dsd2 = [SKL1(1,2,1);SKL1(1,2,2);SKL1(1,2,3)];
            n1 = -cross(dsd2,v31);
    end

    n1 = n1/norm(n1);
    n = [n1(1) 0     0     n1(2) 0     n1(3);
         0     n1(2) 0     n1(1) n1(3) 0;
         0     0     n1(3) 0     n1(2) n1(1)];
    
    invJS1 = inv(JS1);
    F1 = blkdiag(invJS1,invJS1,invJS1);
    G1 = zeros(9,5*ndof1);
    N1 = zeros(3,5*ndof1);
    for i = 1:ndof1
        G1(1:9,(5*i-4):5*i)=[ Rs1(i)   0    0   para(3)*(V1N1(1,i)*Rs1(i))*ti1/2   -para(3)*(V2N1(1,i)*Rs1(i))*ti1/2
                              Rt1(i)   0    0   para(3)*(V1N1(1,i)*Rt1(i))*ti1/2   -para(3)*(V2N1(1,i)*Rt1(i))*ti1/2
                              0      0    0               V1N1(1,i)*R1(i)*ti1/2   -V2N1(1,i)*R1(i)*ti1/2
                              0    Rs1(i)  0   para(3)*(V1N1(2,i)*Rs1(i))*ti1/2   -para(3)*(V2N1(2,i)*Rs1(i))*ti1/2
                              0    Rt1(i)  0   para(3)*(V1N1(2,i)*Rt1(i))*ti1/2   -para(3)*(V2N1(2,i)*Rt1(i))*ti1/2
                              0      0    0               V1N1(2,i)*R1(i)*ti1/2   -V2N1(2,i)*R1(i)*ti1/2
                              0      0  Rs1(i) para(3)*(V1N1(3,i)*Rs1(i))*ti1/2   -para(3)*(V2N1(3,i)*Rs1(i))*ti1/2
                              0      0  Rt1(i) para(3)*(V1N1(3,i)*Rt1(i))*ti1/2   -para(3)*(V2N1(3,i)*Rt1(i))*ti1/2
                              0      0    0               V1N1(3,i)*R1(i)*ti1/2   -V2N1(3,i)*R1(i)*ti1/2];
        N1(1:3,(5*i-4):5*i)=R1(i)*[eye(3) para(3)*ti1/2*V1N1(:,i) -para(3)*ti1/2*V2N1(:,i)];
    end

    B1 = H1*F1*G1;

    [R2, Rindex2]   = nrbbasisfun ({ found_para(1),found_para(2) },ieplane2);
    [Rs2,Rt2] = nrbbasisfunder ({ found_para(1),found_para(2) },ieplane2);
    [SKL2] = RatSurfaceDerivs( ieplane2,found_para(1),found_para(2),2);
    [~, ~ ,v32 ,theta2] = localCoordinate(SKL2);
    [~,~,~,~,v3s2,v3t2] = dVdst(SKL2);
    [T2,H2] = local2global_voigt( theta2 );
    [t2,ts2,tt2] = Thickness([found_para(1),found_para(2)],0,Skin_total{patch_index2});
    V1N2 = V1N2_all(:,Rindex2);
    V2N2 = V2N2_all(:,Rindex2);
    Cps2 = [conpsX2(Rindex2),conpsY2(Rindex2),conpsZ2(Rindex2)];

    rou2 = 1;

    C2 = rou2*DL_Solid2;
    JS21 = Rs2*Cps2 + found_para(3)*v3s2'*(t2/2) + found_para(3)*v32'*(ts2/2);
    JS22 = Rt2*Cps2 + found_para(3)*v3t2'*(t2/2) + found_para(3)*v32'*(tt2/2);
    JS23 = v32'*(t2/2);
    JS2  = [JS21;JS22;JS23];
    invJS2 = inv(JS2);
    F2 = blkdiag(invJS2,invJS2,invJS2);
    G2 = zeros(9,5*ndof2);
    N2 = zeros(3,5*ndof2);
    gs2 = found_para(3);
    for i = 1:ndof2
        G2(1:9,(5*i-4):5*i)=[ Rs2(i)   0    0   gs2*(V1N2(1,i)*Rs2(i))*ti2/2   -gs2*(V2N2(1,i)*Rs2(i))*ti2/2
                              Rt2(i)   0    0   gs2*(V1N2(1,i)*Rt2(i))*ti2/2   -gs2*(V2N2(1,i)*Rt2(i))*ti2/2
                              0      0    0               V1N2(1,i)*R2(i)*ti2/2   -V2N2(1,i)*R2(i)*ti2/2
                              0    Rs2(i)  0   gs2*(V1N2(2,i)*Rs2(i))*ti2/2   -gs2*(V2N2(2,i)*Rs2(i))*ti2/2
                              0    Rt2(i)  0   gs2*(V1N2(2,i)*Rt2(i))*ti2/2   -gs2*(V2N2(2,i)*Rt2(i))*ti2/2
                              0      0    0               V1N2(2,i)*R2(i)*ti2/2   -V2N2(2,i)*R2(i)*ti2/2
                              0      0  Rs2(i) gs2*(V1N2(3,i)*Rs2(i))*ti2/2   -gs2*(V2N2(3,i)*Rs2(i))*ti2/2
                              0      0  Rt2(i) gs2*(V1N2(3,i)*Rt2(i))*ti2/2   -gs2*(V2N2(3,i)*Rt2(i))*ti2/2
                              0      0    0               V1N2(3,i)*R2(i)*ti2/2   -V2N2(3,i)*R2(i)*ti2/2];

        N2(1:3,(5*i-4):5*i)=R2(i)*[eye(3) gs2*ti2/2*V1N2(:,i) -gs2*ti2/2*V2N2(:,i)];
    end

    B2 = H2*F2*G2;

    elementDof1_son{loop_guass} = sort([5*Rindex1,5*Rindex1-4,5*Rindex1-3,5*Rindex1-2,5*Rindex1-1]);
    elementDof2_son{loop_guass} = sort([5*Rindex2,5*Rindex2-4,5*Rindex2-3,5*Rindex2-2,5*Rindex2-1]);
    Kn_upper_left_son{loop_guass}  = transpose(N1)*n*C1*B1*fac1;
    Kn_upper_right_son{loop_guass} = transpose(N1)*n*C2*B2*fac1;
    Kn_lower_left_son{loop_guass}  = transpose(N2)*n*C1*B1*fac1;
    Kn_lower_right_son{loop_guass} = transpose(N2)*n*C2*B2*fac1;
    Ks_upper_left_son{loop_guass}  = transpose(N1)*N1*fac1;
    Ks_upper_right_son{loop_guass} = transpose(N1)*N2*fac1;
    Ks_lower_left_son{loop_guass}  = transpose(N2)*N1*fac1;
    Ks_lower_right_son{loop_guass} = transpose(N2)*N2*fac1;
end

for loop_guass = 1:length(temp_coupling_points)
    elementDof1 = elementDof1_son{loop_guass};
    elementDof2 = elementDof2_son{loop_guass};
    Kn_upper_left(elementDof1,elementDof1)  = Kn_upper_left(elementDof1,elementDof1)+Kn_upper_left_son{loop_guass};
    Kn_upper_right(elementDof1,elementDof2) = Kn_upper_right(elementDof1,elementDof2)+Kn_upper_right_son{loop_guass};
    Kn_lower_left(elementDof2,elementDof1)  = Kn_lower_left(elementDof2,elementDof1)+Kn_lower_left_son{loop_guass};
    Kn_lower_right(elementDof2,elementDof2) = Kn_lower_right(elementDof2,elementDof2)+Kn_lower_right_son{loop_guass};
    Ks_upper_left(elementDof1,elementDof1)  = Ks_upper_left(elementDof1,elementDof1)+Ks_upper_left_son{loop_guass};
    Ks_upper_right(elementDof1,elementDof2) = Ks_upper_right(elementDof1,elementDof2)+Ks_upper_right_son{loop_guass};
    Ks_lower_left(elementDof2,elementDof1)  = Ks_lower_left(elementDof2,elementDof1)+Ks_lower_left_son{loop_guass};
    Ks_lower_right(elementDof2,elementDof2) = Ks_lower_right(elementDof2,elementDof2)+Ks_lower_right_son{loop_guass};
end
Kn_upper_left  = sparse(Kn_upper_left);
Kn_upper_right = sparse(Kn_upper_right);
Kn_lower_left  = sparse(Kn_lower_left);
Kn_lower_right = sparse(Kn_lower_right);
Ks_upper_left  = sparse(Ks_upper_left);
Ks_upper_right = sparse(Ks_upper_right);
Ks_lower_left  = sparse(Ks_lower_left);
Ks_lower_right = sparse(Ks_lower_right);

if(patch_index1<=patch_index2)
    Kn = [-gamma*Kn_upper_left -(1-gamma)*Kn_upper_right;
        gamma*Kn_lower_left   (1-gamma)*Kn_lower_right];
    Ks = [alpha*Ks_upper_left -alpha*Ks_upper_right;
        -alpha*Ks_lower_left  alpha*Ks_lower_right];
else
    Kn = [(1-gamma)*Kn_lower_right  gamma*Kn_lower_left ;
        -(1-gamma)*Kn_upper_right -gamma*Kn_upper_left];
    Ks = [alpha*Ks_lower_right -alpha*Ks_lower_left;
        -alpha*Ks_upper_right  alpha*Ks_upper_left];
end
interfacial_stiffness = Kn+transpose(Kn)+Ks;
end