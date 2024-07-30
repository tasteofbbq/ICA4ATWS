function [ elementNode,noElement ] = GetParaMesh_Surface( kontU,kontV,p,q,F_parameter)

U = unique(kontU(p+1:end-p));
V = unique(kontV(q+1:end-q));
noElementU = length(U)-1;
noElementV = length(V)-1;
noElement = noElementU*noElementV;
elementNode.vertex = cell(1,noElement);
elementNode.center = cell(1,noElement);
elementNode.r_in = zeros(1,noElement);
elementNode.r_out = zeros(1,noElement);
for i = 1:noElementV
    for j = 1:noElementU
        node = zeros(4,2);
        indexElement = (i-1)*noElementU+j;
        node(1,1) = U(j);
        node(2,1) = U(j+1);
        node(3,1) = node(2,1);
        node(4,1) = node(1,1);
        
        node(1,2) = V(i);
        node(2,2) = node(1,2);
        node(3,2) = V(i+1);
        node(4,2) = node(3,2);
        
        node(:,3) = 0;
        
        center(1) = ( U(j)+U(j+1) )/2;
        center(2) = ( V(i)+V(i+1) )/2;
        center(3) = 0;
        
        elementNode.vertex{1,indexElement} = node;
        elementNode.center{1,indexElement} = center;
        
        b = (U(j+1)-U(j))/2; h=(V(i+1)-V(i))/2;
        elementNode.r_in (1,indexElement)= min( b,h);
        elementNode.r_out(1,indexElement) = sqrt(b^2+h^2);
    end
end

if F_parameter~=0
    figure(F_parameter)
    for k = 1:noElementV+1
        line( [U(1),U(end)], [V(k),V(k)], 'LineWidth', 1.2, 'color', 'k' )
    end

    for k = 1:noElementU+1
        line( [U(k),U(k)],[V(1),V(end)], 'LineWidth', 1.2, 'color', 'k' )
    end
    axis equal;
end

end

