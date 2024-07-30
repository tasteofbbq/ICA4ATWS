function [ GaussPoint,weight ] = GaussQ36( a_U, b_U, a_V, b_V )
GaussPoint=zeros(36,2);
weight=zeros(36,1);
xi = [-0.932469514203152,-0.661209386466265,-0.238619186083197,0.238619186083197,0.661209386466265,0.932469514203152];
H = [0.171324492379170,0.360761573048139,0.467913934572691];
xi_1 = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
xi_2 = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
xi_3 = (a_U+b_U)/2 + (b_U-a_U)*xi(3)/2;
xi_4 = (a_U+b_U)/2 + (b_U-a_U)*xi(4)/2;
xi_5 = (a_U+b_U)/2 + (b_U-a_U)*xi(5)/2;
xi_6 = (a_U+b_U)/2 + (b_U-a_U)*xi(6)/2;
eta_1 = (a_V+b_V)/2 + (b_V-a_V)*xi(1)/2;
eta_2 = (a_V+b_V)/2 + (b_V-a_V)*xi(2)/2;
eta_3 = (a_V+b_V)/2 + (b_V-a_V)*xi(3)/2;
eta_4 = (a_V+b_V)/2 + (b_V-a_V)*xi(4)/2;
eta_5 = (a_V+b_V)/2 + (b_V-a_V)*xi(5)/2;
eta_6 = (a_V+b_V)/2 + (b_V-a_V)*xi(6)/2;

HU_1 = (b_U-a_U)*H(1)/2; HU_2 = (b_U-a_U)*H(2)/2; HU_3 = (b_U-a_U)*H(3)/2;
HV_1 = (b_V-a_V)*H(1)/2; HV_2 = (b_V-a_V)*H(2)/2; HV_3 = (b_V-a_V)*H(3)/2;

GaussPoint(:,1) = [xi_1,xi_2,xi_3,xi_4,xi_5,xi_6,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6,...
    xi_1,xi_2,xi_3,xi_4,xi_5,xi_6,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6,...
    xi_1,xi_2,xi_3,xi_4,xi_5,xi_6,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6];
GaussPoint(:,2) = [eta_1,eta_1,eta_1,eta_1,eta_1,eta_1,eta_2,eta_2,eta_2,eta_2,eta_2,eta_2,...
    eta_3,eta_3,eta_3,eta_3,eta_3,eta_3,eta_4,eta_4,eta_4,eta_4,eta_4,eta_4,...
    eta_5,eta_5,eta_5,eta_5,eta_5,eta_5,eta_6,eta_6,eta_6,eta_6,eta_6,eta_6];

weight(1:6) = HV_1*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];
weight(7:12) = HV_2*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];
weight(13:18) = HV_3*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];
weight(19:24) = HV_3*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];
weight(25:30) = HV_2*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];
weight(31:36) = HV_1*[ HU_1, HU_2, HU_3, HU_3, HU_2, HU_1];




end

