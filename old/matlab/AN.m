
function A = AN(N)
   A = ones(N,N);
   for k = 1:N-1,
      A(1:k,k) = (N+1)-k;
      A(k,1:k) = (N+1)-k;
   end
 
