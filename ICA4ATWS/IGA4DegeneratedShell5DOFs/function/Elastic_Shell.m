function DL = Elastic_Shell(E,nu)
DL = E/(1-nu^2)*[  1   nu   0      0           0               0;
                   nu   1   0      0           0               0;
                   0    0   0      0           0               0;
                   0    0   0  (1-nu)/2        0               0;
                   0    0   0      0     5/6*(1-nu)/2         0;
                   0    0   0      0           0         5/6*(1-nu)/2];
               
end
