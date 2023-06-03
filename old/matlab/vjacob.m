function df = vjacob(func,v,p)
   h = 1e-6;
   n = length(v);
   f0 = feval(func,v,p);
   df = zeros(n,n);
   for k = 1:n,
       dv = h*((1:n)'==k);
       df(1:n,k) = (feval(func,v+dv,p)-f0)/h;
   end
   
