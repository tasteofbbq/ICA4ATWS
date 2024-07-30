function [CP] = ControlPoint(conpsX,conpsY,conpsZ,Ncops)
CP = cell(1,Ncops);

for i = 1:Ncops
    controlP(1,1) = conpsX(i);
    controlP(1,2) = conpsY(i);
    controlP(1,3) = conpsZ(i);
    CP{1,i} = controlP;
end