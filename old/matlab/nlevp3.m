function f = nlevp3(x,D,B)
   N = size(x,1)-1;
   omegasq = exp(x(N+1));
   f = [omegasq*B*sin(x(1:N))-D*tan(x(1:N)); 0];
