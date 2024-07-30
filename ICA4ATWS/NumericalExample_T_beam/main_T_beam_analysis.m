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

%% Clear all
clear
clc
close all
format short e

%% Number of elements
NuandNv(1,:) = [100 10];
NuandNv(2,:) = [10 100];

%% Material property
E = 1e7;
miu = 0;

%% Import geometry
load('patches_T_beam.mat')

%% Parameter domain normalization([a, b]->[0, 1])
for loop_nbsf = 1:length(nurbssrfs)
    nurbssrfs{loop_nbsf,1}.knots{1,1} = (nurbssrfs{loop_nbsf,1}.knots{1,1}-min(nurbssrfs{loop_nbsf,1}.knots{1,1}))/(max(nurbssrfs{loop_nbsf,1}.knots{1,1})-min(nurbssrfs{loop_nbsf,1}.knots{1,1}));
    nurbssrfs{loop_nbsf,1}.knots{1,2} = (nurbssrfs{loop_nbsf,1}.knots{1,2}-min(nurbssrfs{loop_nbsf,1}.knots{1,2}))/(max(nurbssrfs{loop_nbsf,1}.knots{1,2})-min(nurbssrfs{loop_nbsf,1}.knots{1,2}));
end

%% Pre-processing
Total_Skin = cell(1,length(nurbssrfs));
for loop_nbsf = 1:length(nurbssrfs)
    Total_Skin{loop_nbsf} = skin_operation_TTO_multipatch(nurbssrfs{loop_nbsf},NuandNv(loop_nbsf,:),E,miu,loop_nbsf);
    fprintf('No.%d patch has been operated!\n', loop_nbsf)
end

%% Coupling_information
coupling_information_total

%% Boundary condition
activeDof = Boundary_condition_multi(Total_Skin,coupling_information);

%% Load condition
f = ForceVector_condition_multi(Total_Skin,coupling_information);

%% Coupling_type: 1 Nitsche_coupling 2 Penalty_coupling 3 Mortar_coupling
coupling_type = 1;
coupling_points = find_coupling_points(Total_Skin, coupling_information);
switch coupling_type
    case 1
        alpha_nitsche = 2e8;
        gamma_nitsche = 0.5;
        [multi_K, multi_K_b] = Nitsche_coupling_guass(Total_Skin, coupling_information, coupling_points, gamma_nitsche, alpha_nitsche);
    case 2
        alpha_penalty = 1e8;
        [multi_K, multi_K_b] = Penalty_coupling_guass(Total_Skin, coupling_information, coupling_points, alpha_penalty);
    case 3
        [multi_K, multi_K_b] = Mortar_coupling_guass(Total_Skin, coupling_information, coupling_points);
        f = [f;sparse(size(multi_K,1)-size(multi_K_b,1),1)];
        activeDof = [activeDof size(multi_K_b,1)+1:size(multi_K,1)];
end

%% Solve the displacement (static analysis)
Disp = multi_K(activeDof,activeDof)\f(activeDof);
Displacement = zeros(size(multi_K,1),1);
Displacement(activeDof) = Disp;
disp('Solve the displacement finish')
fprintf('coupled energy value = %f\n', Displacement'*multi_K*Displacement/2)
fprintf('non-coupled energy value = %f\n', (Displacement(1:size(multi_K_b,1)))'*multi_K_b*Displacement(1:size(multi_K_b,1))/2)
Total_Skin = Displacement_operation(Total_Skin, Displacement(1:size(multi_K_b,1)));

%% Displacement at the loaded point
[load_u,load_v,load_w] = Evaluate_Displacement([0 1 0],Total_Skin{2});
fprintf('Solve the displacement at the loaded point = %f\n', abs(load_v))

%% Linear buckling analysis
multi_K_G = K_G_multi(Total_Skin);
mode  = 5;
multi_K_G = blkdiag(multi_K_G,sparse(size(multi_K,1)-size(multi_K_b,1),size(multi_K,1)-size(multi_K_b,1)));
[active_fi, lambda] = eigs(multi_K(activeDof,activeDof),-multi_K_G(activeDof,activeDof),mode,'smallestabs');
fprintf('Solve the first mode eigen value = %f\n', lambda(1,1))