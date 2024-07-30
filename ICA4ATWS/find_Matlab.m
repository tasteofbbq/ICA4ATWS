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

function [found_para, val] = find_Matlab(X,Skin,local_surface_index) 
% X                      global coordinate on patch 1
% Skin                   information of Skin of patch 2
% local_surface_index    local index of patch 2
%% different initial value for patch 2
conpsX_temp  = Skin.conpsX;
conpsY_temp  = Skin.conpsY;
conpsZ_temp  = Skin.conpsZ;
for loop_conps = 1:length(conpsX_temp)
    pre_distance = norm(X-[conpsX_temp(loop_conps) conpsY_temp(loop_conps) conpsZ_temp(loop_conps)]);
    if(loop_conps == 1)
        pre_distance_temp = pre_distance;
        index_conps = loop_conps;
        continue
    end
    if(pre_distance<pre_distance_temp)
        index_conps = loop_conps;
        pre_distance_temp = pre_distance;
    end
end

%% set initial value according to pre_para
switch local_surface_index
    case 1
        pre_para = [0.5 0.5];
        para0 = [pre_para 0];
    case 2
        pre_para = [0.5 0.5];
        para0 = [pre_para 0];
    case 3
    case 4
        pre_para = [1 0.5];
        para0 = [pre_para 0];
    case 5
    case 6
        pre_para = [0 0.5];
        para0 = [pre_para 0];
    case 7
        pre_para = [0.5 0.5];
        para0 = [pre_para 0];
end

[t,ts,tt] = Thickness([0.5,0.5],0,Skin);

    function [funeval,grad] = gradfun(X,Skin, para0)
        ieplane = Skin.ieplane;
        conpsX  = Skin.conpsX;
        conpsY  = Skin.conpsY;
        conpsZ  = Skin.conpsZ;

        knots    = ieplane.knots;
        unique_knots = ieplane.unique_knots;
        if para0(1) == unique_knots{1}(end)
            index1 = length(unique_knots{1})-1;
        else
            index1 = find(unique_knots{1} > para0(1), 1)-1;
        end
        if para0(2) == unique_knots{2}(end)
            index2 = length(unique_knots{2})-1;
        else
            index2 = find(unique_knots{2} > para0(2), 1)-1;
        end
        index_ele = (length(unique_knots{1})-1)*(index2-1)+index1;
        Rindex = Skin.EL_nodeId(index_ele,:);
        degree   = ieplane.order-1;
        CP_wghts = Skin.CP_wghts;
        EN_wghts = CP_wghts(Rindex,1);
        Rders = Nurbs2DBasisFun2ndDer_Cpp(EN_wghts,knots{1},knots{2},degree(1),degree(2),[para0(1);para0(2)]);
        R     = Rders(1,:);
        Rs    = Rders(2,:);
        Rt    = Rders(3,:);
        Rss   = Rders(4,:);
        Rtt   = Rders(5,:);
        Rst   = Rders(6,:);
        conpsXRindex = conpsX(Rindex);
        SKL(:,:,1) = [R*conpsXRindex   Rt*conpsXRindex  Rtt*conpsXRindex;
                      Rs*conpsXRindex  Rst*conpsXRindex 0;
                      Rss*conpsXRindex 0                0;];
        conpsYRindex = conpsY(Rindex);
        SKL(:,:,2) = [R*conpsYRindex   Rt*conpsYRindex  Rtt*conpsYRindex;
                      Rs*conpsYRindex  Rst*conpsYRindex 0;
                      Rss*conpsYRindex 0                0;];
        conpsZRindex = conpsZ(Rindex);
        SKL(:,:,3) = [R*conpsZRindex   Rt*conpsZRindex  Rtt*conpsZRindex;
                      Rs*conpsZRindex  Rst*conpsZRindex 0;
                      Rss*conpsZRindex 0                0;];
        [~, ~ ,v3 ,~] = localCoordinate(SKL);
        [~,~,~,~,v3s,v3t] = dVdst(SKL);

        Cps = [conpsX(Rindex),conpsY(Rindex),conpsZ(Rindex)];
        JS(1,:) = Rs*Cps + para0(3)*v3s'*(t/2) + para0(3)*v3'*(ts/2);
        JS(2,:) = Rt*Cps + para0(3)*v3t'*(t/2) + para0(3)*v3'*(tt/2);
        JS(3,:) = v3'*(t/2);
        point_on_2 = R*Cps + para0(3)*v3'*(t/2);
        funeval = norm(point_on_2-X,2)^2;
        grad = 2*JS*((point_on_2-X)');
        grad = grad';
    end
    function [funeval] = fun(X,Skin, para0)
        ieplane = Skin.ieplane;
        conpsX  = Skin.conpsX;
        conpsY  = Skin.conpsY;
        conpsZ  = Skin.conpsZ;

        knots    = ieplane.knots;
        unique_knots = ieplane.unique_knots;
        if para0(1) == unique_knots{1}(end)
            index1 = length(unique_knots{1})-1;
        else
            index1 = find(unique_knots{1} > para0(1), 1)-1;
        end
        if para0(2) == unique_knots{2}(end)
            index2 = length(unique_knots{2})-1;
        else
            index2 = find(unique_knots{2} > para0(2), 1)-1;
        end
        index_ele = (length(unique_knots{1})-1)*(index2-1)+index1;
        Rindex = Skin.EL_nodeId(index_ele,:);
        degree   = ieplane.order-1;
        CP_wghts = Skin.CP_wghts;
        EN_wghts = CP_wghts(Rindex,1);
        Rders = Nurbs2DBasisFun2ndDer_Cpp(EN_wghts,knots{1},knots{2},degree(1),degree(2),[para0(1);para0(2)]);
        R     = Rders(1,:);
        Rs    = Rders(2,:);
        Rt    = Rders(3,:);
        Rss   = Rders(4,:);
        Rtt   = Rders(5,:);
        Rst   = Rders(6,:);
        conpsXRindex = conpsX(Rindex);
        SKL(:,:,1) = [R*conpsXRindex   Rt*conpsXRindex  Rtt*conpsXRindex;
                      Rs*conpsXRindex  Rst*conpsXRindex 0;
                      Rss*conpsXRindex 0                0;];
        conpsYRindex = conpsY(Rindex);
        SKL(:,:,2) = [R*conpsYRindex   Rt*conpsYRindex  Rtt*conpsYRindex;
                      Rs*conpsYRindex  Rst*conpsYRindex 0;
                      Rss*conpsYRindex 0                0;];
        conpsZRindex = conpsZ(Rindex);
        SKL(:,:,3) = [R*conpsZRindex   Rt*conpsZRindex  Rtt*conpsZRindex;
                      Rs*conpsZRindex  Rst*conpsZRindex 0;
                      Rss*conpsZRindex 0                0;];
        [~, ~ ,v3 ,~] = localCoordinate(SKL);

        Cps = [conpsX(Rindex),conpsY(Rindex),conpsZ(Rindex)];
        point_on_2 = R*Cps + para0(3)*v3'*(t/2);
        funeval = norm(point_on_2-X,2)^2;
    end
%% searching option
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 -1];
ub = [1 1 1];
nonlcon = [];

%% fmincon to optimize
options=optimset('GradObj','on','display','off', ...
    'TolFun',1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4);
[found_para, val] = fmincon(@(para)gradfun(X, Skin, para), para0, A, b, Aeq, beq, lb, ub, nonlcon, options);
val = sqrt(val);

end
