function [ GaussPoint,weight ] = GaussQ16( a_U, b_U, a_V, b_V )
GaussPoint=zeros(16,2);
weight=zeros(16,1);
xi=[ -0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053 ];
H = [0.347854845137454,0.625145154862546];
xi_1 = (a_U+b_U)/2 + (b_U-a_U)*xi(1)/2;
xi_2 = (a_U+b_U)/2 + (b_U-a_U)*xi(2)/2;
xi_3 = (a_U+b_U)/2 + (b_U-a_U)*xi(3)/2;
xi_4 = (a_U+b_U)/2 + (b_U-a_U)*xi(4)/2;
eta_1 = (a_V+b_V)/2 + (b_V-a_V)*xi(1)/2;
eta_2 = (a_V+b_V)/2 + (b_V-a_V)*xi(2)/2;
eta_3 = (a_V+b_V)/2 + (b_V-a_V)*xi(3)/2;
eta_4 = (a_V+b_V)/2 + (b_V-a_V)*xi(4)/2;
HU_1 = (b_U-a_U)*H(1)/2; HU_2 = (b_U-a_U)*H(2)/2;
HV_1 = (b_V-a_V)*H(1)/2; HV_2 = (b_V-a_V)*H(2)/2;
GaussPoint(:,1) = [xi_1,xi_2,xi_3,xi_4,xi_1,xi_2,xi_3,xi_4,...
                   xi_1,xi_2,xi_3,xi_4,xi_1,xi_2,xi_3,xi_4];
GaussPoint(:,2) = [eta_1,eta_1,eta_1,eta_1,eta_2,eta_2,eta_2,eta_2,...
                   eta_3,eta_3,eta_3,eta_3,eta_4,eta_4,eta_4,eta_4];
weight(1:4) = HV_1*[ HU_1, HU_2, HU_2, HU_1];
weight(5:8) = HV_2*[ HU_1, HU_2, HU_2, HU_1];
weight(9:12) = HV_2*[ HU_1, HU_2, HU_2, HU_1];
weight(13:16) = HV_1*[ HU_1, HU_2, HU_2, HU_1];

end

