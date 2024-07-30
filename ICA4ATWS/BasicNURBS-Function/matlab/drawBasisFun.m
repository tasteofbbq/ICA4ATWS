U = [ 0 0 0.2 0.4 0.6 0.8 1 1 ];
n = 5; p =1; m = 7;
u = (0:0.01:1);
N = zeros(n+1,101);
for i = 1:n+1
    for j=1:101
        N(i,j) = OneBasisFun( p,m,U,i,u(j) );
    end
end

figure
for i = 1:n+1
    plot(u,N(i,:),'-')
    hold on
end