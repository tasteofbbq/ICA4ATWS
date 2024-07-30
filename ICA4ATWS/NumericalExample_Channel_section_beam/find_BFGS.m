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

function [found_para, val] = find_BFGS(X,Skin,local_surface_index)
% X                      global coordinate on patch 1
% Skin                   information of Skin of patch 2
% local_surface_index    local index of patch 2
%% different initial value for patch 2
switch local_surface_index
    case 1
    case 2
    case 3
        para0 = [0.5,0.5,1];
    case 4
        para0 = [1,0.5,0];
    case 5
    case 6
        para0 = [0,0.5,0];
    case 7
        para0 = [0.5,0.5,0];
end
    function [funeval,grad] = gradfun(X,Skin, para0)
        ieplane = Skin.ieplane;
        conpsX  = Skin.conpsX;
        conpsY  = Skin.conpsY;
        conpsZ  = Skin.conpsZ;
        [R, Rindex]   = nrbbasisfun ({ para0(1),para0(2) },ieplane);
        [Rs,Rt] = nrbbasisfunder ({ para0(1),para0(2) },ieplane);
        [SKL] = RatSurfaceDerivs( ieplane,para0(1),para0(2),2);
        [~, ~ ,v3 ,theta] = localCoordinate(SKL);
        [~,~,~,~,v3s,v3t] = dVdst(SKL);
        [t,ts,tt] = Thickness([para0(1),para0(2)],0,Skin);
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
        [R, Rindex]   = nrbbasisfun ({ para0(1),para0(2) },ieplane);
        [SKL] = RatSurfaceDerivs( ieplane,para0(1),para0(2),2);
        [~, ~ ,v3 ,~] = localCoordinate(SKL);
        [t,~,~] = Thickness([para0(1),para0(2)],0,Skin);
        Cps = [conpsX(Rindex),conpsY(Rindex),conpsZ(Rindex)];
        point_on_2 = R*Cps + para0(3)*v3'*(t/2);
        funeval = norm(point_on_2-X,2)^2;
    end
%% Armjio and BFGS
k = 0;
maxk = 1000;
rho = 0.55;
sigma = 4e-1;
e = 1e-5;
n = length(para0);
Hk = eye(n);
x0 = para0;
x0 = x0';
sk = 1000000;
    function x = trans01(x)
        for i = 1:size(x,2)
            if(i<3)
                if(x(1,i)<=0)
                    x(1,i) = 0;
                end
            end
            if(i==3)
                if(x(1,i)<=-1)
                    x(1,i) = -1;
                end
            end
            if(x(1,i)>=1)
                x(1,i) = 1;
            end
        end
    end
while(k<maxk)
    [~,gk] = gradfun(X,Skin, x0');
    gk = gk';

    if(max(abs(sk))<e),break;end
    dk = -Hk\gk;
    m = 0;
    mk = 0;

    [temp_a] = fun(X,Skin, x0');
    while(m<100)
        temp = x0'+rho^m*(dk');
        if(~((0<=temp(1)&&temp(1)<=1)||(0<=temp(2)&&temp(2)<=1)||(-1<=temp(3)&&temp(3)<=1)))
            m = m+1;
            continue;
        end
        [s] = fun(X,Skin, trans01(x0'+rho^m*(dk')));
        a = temp_a+sigma*rho^m*gk'*dk;
        if(s<a)
            mk=m;
            break;
        end
        m = m+1;
    end
    x = trans01((x0+dk*rho^mk)')';

    sk = x-x0;
    [~,yk] = gradfun(X,Skin, x');
    yk = yk';
    yk = yk-gk;
    if(yk'*sk>0)
        Hk = Hk-(Hk*(sk*sk')*Hk)/(sk'*Hk*sk)+(yk*yk')/(yk'*yk);
    end
    x0 = x;
    k = k+1;
end
found_para = x0';

val = fun(X,Skin, found_para);
val = sqrt(val);
end