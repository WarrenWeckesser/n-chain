function z = releqfunc(phi,omega)
   n = length(phi);
   
   z = AN(n)*sin(phi) - diag(n:-1:1)*tan(phi)/omega^2;
