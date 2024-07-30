function [Skin] = Local_Coordinate_Surface(Skin)
nurbsSurf = Skin.ieplane;
CPw = Skin.CPw;
Nodenum = nurbsSurf.number(1,1)*nurbsSurf.number(1,2);
Skin.LocalCoor = cell(1,Nodenum);

for i = 1:Nodenum
    ID = i;
                           
    [SKL] = RatSurfaceDerivs( nurbsSurf,CPw{1,ID}(1),CPw{1,ID}(2),2);
    [v1, v2 ,v3 ,~] = localCoordinate(SKL);
    
    Skin.LocalCoor{1,i}(:,1) = v1;
    Skin.LocalCoor{1,i}(:,2) = v2;
    Skin.LocalCoor{1,i}(:,3) = v3;
    
    Skin.CPz(:,i) =  Skin.CP{i}';
    Skin.CPwz(:,i) = CPw{1,ID}';
    Skin.LCSv1(:,i) = v1;
    Skin.LCSv2(:,i) = v2;
    Skin.LCSv3(:,i) = v3;   
end