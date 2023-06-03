function f = nlevp(x,omegasq,D,B)
   f = omegasq*B*sin(x)-D*tan(x);
