function [ GaussPoint,weight ] = GaussQ4( a_U, b_U, a_V, b_V )
GaussPoint=zeros(4,2);
weight=zeros(4,1);
xi=[ -0.577350269189626 ,0.577350269189626 ];
GaussPoint(1,1) = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
GaussPoint(1,2) = (a_V+b_V)/2 + (b_V-a_V)*xi(1)/2;
GaussPoint(2,1) = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
GaussPoint(2,2) = (a_V+b_V)/2 + (b_V-a_V)*xi(1)/2;
GaussPoint(3,1) = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
GaussPoint(3,2) = (a_V+b_V)/2 + (b_V-a_V)*xi(2)/2;
GaussPoint(4,1) = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
GaussPoint(4,2) = (a_V+b_V)/2 + (b_V-a_V)*xi(2)/2;
weight(:,1)= (b_U-a_U)*(b_V-a_V)/4;
end

