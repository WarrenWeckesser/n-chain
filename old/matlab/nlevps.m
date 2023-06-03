function f = nlevps(x,D,B)
   N = size(x,1)-1;
   omegasq = x(N+1)^4;
   f = [omegasq*B*sin(x(1:N))-D*tan(x(1:N)); 0];
