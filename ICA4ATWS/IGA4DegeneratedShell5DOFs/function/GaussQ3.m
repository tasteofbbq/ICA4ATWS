function [ GaussPoint,weight ] = GaussQ3( a_U, b_U )
%GAUSSQ3 Summary of this function goes here
GaussPoint=zeros(3,1);
weight=zeros(3,1);
xi=[ -0.774596669241483,0, 0.774596669241483];
H = [0.555555555555556,0.888888888888889];

xi_1 = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
xi_2 = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
xi_3 = (a_U+b_U)/2 + (b_U-a_U)*xi(3)/2;

HU_1 = (b_U-a_U)*H(1)/2; HU_2 = (b_U-a_U)*H(2)/2;

GaussPoint = [xi_1;xi_2;xi_3];
weight = [HU_1;HU_2;HU_1];

end

