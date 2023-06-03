function plotnchain(phi)
   n = length(phi);
   clf
   x0 = 0;
   y0 = 0;
   r = 1/n;
   plot([-.1 .1],[0 0],'b');
   hold on
   plot([0 0],[-.1 .1],'b');
   for k = 1:n,
       x1 = x0 + r*sin(phi(k));
       y1 = y0 - r*cos(phi(k));
       plot([x0 x1],[y0 y1],'k');
       plot(x1,y1,'k.','Linewidth',[2]);
       x0 = x1;
       y0 = y1;
   end
   axis equal
   axis([-1 1 -1 1])
   axis off
   
