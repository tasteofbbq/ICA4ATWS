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

function [ K ] = K_for_Skin(Parameter, DL, t,theta)
ieplane = Parameter.ieplane;
elementNode = Parameter.elementNode;
noElement = Parameter.noElement;
NoCP = ieplane.number(1,1)*ieplane.number(1,2);
dof = 5*NoCP; 
[ element_CP ] = elementCP( ieplane,elementNode.center,noElement );

%% initial global stiffness matrix
K = sparse(dof,dof);
elementType = zeros(1,noElement);

%% get the index of elements
fullEle = find(elementType == 0);
nofullEle = size(fullEle,2);

%% assemble global stiffness matrix
for i = 1:nofullEle
    ind = fullEle(i);
    ID = element_CP(ind,:);
    elementDof = sort([5*ID,5*ID-4,5*ID-3,5*ID-2,5*ID-1]);
    s1 = elementNode.vertex{1, ind}(1,1);
    s2 = elementNode.vertex{1, ind}(2,1);
    t1 = elementNode.vertex{1, ind}(2,2);
    t2 = elementNode.vertex{1, ind}(3,2);
    
    [ GaussPoint,Weight_st ] = GaussQ9( s1, s2, t1, t2 );

    [ stiffness ] = K_Element_Skin_STO(Parameter,GaussPoint,Weight_st,DL,t,theta,i);
    stiffness_sparse = sparse(dof,dof);
    stiffness_sparse(elementDof,elementDof) = stiffness;
    K = K + stiffness_sparse;
end

end