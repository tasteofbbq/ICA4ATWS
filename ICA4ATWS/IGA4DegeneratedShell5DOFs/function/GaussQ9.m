function [ GaussPoint,weight  ] = GaussQ9( a_U, b_U, a_V, b_V )
%GAUSS9 Summary of this function goes here
GaussPoint=zeros(9,2);
weight=zeros(9,1);
xi=[ -0.774596669241483,0, 0.774596669241483];
H = [0.555555555555556,0.888888888888889];
xi_1 = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
xi_2 = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
xi_3 = (a_U+b_U)/2 + (b_U-a_U)*xi(3)/2;
eta_1 = (a_V+b_V)/2 + (b_V-a_V)*xi(1)/2;
eta_2 = (a_V+b_V)/2 + (b_V-a_V)*xi(2)/2;
eta_3 = (a_V+b_V)/2 + (b_V-a_V)*xi(3)/2;

HU_1 = (b_U-a_U)*H(1)/2; HU_2 = (b_U-a_U)*H(2)/2;
HV_1 = (b_V-a_V)*H(1)/2; HV_2 = (b_V-a_V)*H(2)/2;

GaussPoint(:,1) = [xi_1,xi_2,xi_3,xi_1,xi_2,xi_3,...
                   xi_1,xi_2,xi_3];
GaussPoint(:,2) = [eta_1,eta_1,eta_1,eta_2,eta_2,eta_2,...
                  eta_3,eta_3,eta_3];
              
weight(1:3) = HV_1*[ HU_1, HU_2, HU_1];
weight(4:6) = HV_2*[ HU_1, HU_2, HU_1];
weight(7:9) = HV_1*[ HU_1, HU_2, HU_1];

end

