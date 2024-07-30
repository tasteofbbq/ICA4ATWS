function [ B1, B2, B3, B4 ] = NodesOnBoundary( n, m, noCP )

B1 = (1:n+1);
B4 = (noCP-n:noCP);
B2 = (1:n+1:noCP-n);
B3 = (n+1:n+1:noCP);
        
  

end

