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

function DL = Elastic_Solid(E,nu)
DL = E*(1-nu)/((1+nu)*(1-2*nu))*[  1         nu/(1-nu) nu/(1-nu) 0                   0                   0;
                                   nu/(1-nu) 1         nu/(1-nu) 0                   0                   0;
                                   nu/(1-nu) nu/(1-nu) 1         0                   0                   0;
                                   0         0         0         (1-2*nu)/(2*(1-nu)) 0                   0;
                                   0         0         0         0                   (1-2*nu)/(2*(1-nu)) 0;
                                   0         0         0         0                   0                   (1-2*nu)/(2*(1-nu))];
               
end
