function A = nmatrix(n)
   for i = 1:n,
      for j = 1:n,
         A(i,j) = (n+1)-max([i j]);
      end,
   end
