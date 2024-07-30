function [ GaussPoint,weight  ] = GaussQ6( a_U, b_U )
%GAUSSQ6 Summary of this function goes here
GaussPoint=zeros(6,1);
weight=zeros(6,1);
xi = [-0.932469514203152,-0.661209386466265,-0.238619186083197,0.238619186083197,0.661209386466265,0.932469514203152];
H = [0.171324492379170,0.360761573048139,0.467913934572691];
xi_1 = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
xi_2 = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
xi_3 = (a_U+b_U)/2 + (b_U-a_U)*xi(3)/2;
xi_4 = (a_U+b_U)/2 + (b_U-a_U)*xi(4)/2;
xi_5 = (a_U+b_U)/2 + (b_U-a_U)*xi(5)/2;
xi_6 = (a_U+b_U)/2 + (b_U-a_U)*xi(6)/2;

HU_1 = (b_U-a_U)*H(1)/2; HU_2 = (b_U-a_U)*H(2)/2; HU_3 = (b_U-a_U)*H(3)/2;

GaussPoint = [xi_1;xi_2;xi_3;xi_4;xi_5;xi_6];
weight = [HU_1;HU_2;HU_3;HU_3;HU_2;HU_1];
end

