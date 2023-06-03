curve = ReadNChainSolution('autodata/s.n07phiM004');

N = 7;
k = 4;
m = size(curve,1);

omega = curve(:,2);
phi = curve(:,[1 3:(N+1)]);

cosphi = cos(phi(:,1:N));
sinphi = sin(phi(:,1:N));

x =  cumsum(sinphi,2);
y = -cumsum(cosphi,2);

yn = input('Mark the masses with circles in the plot? ','s');


figure(1)
clf
hold on
for i = 1:m,
   plot([0 x(i,:)],[0 y(i,:)],'k')
   if yn == 'y' 
      plot(x(i,:),y(i,:),'bo')   
   end
end
xmax = max(max(abs(x)));
axis equal
axis([-xmax-.1 xmax+.1 -N 0])
ts = sprintf('Whirling Mode %d of the %d-Chain',k,N);
title(ts)

yn = input('Make a movie? ','s');
if (yn == 'y')
   figure(2)
   clf
   hold off
   yn = input('Mark the masses with circles in the movie? ','s');
   xmax = max(max(abs(x)));
   hs = plot([0 x(1,:)],[0 y(1,:)],'k','Linewidth',2);
   set(hs,'EraseMode','xor');
   hold on
   if yn == 'y'
      hm = plot(x(1,:),y(1,:),'bo','Linewidth',2);
      set(hm,'EraseMode','xor');
   end
   axis equal
   axis([-xmax-.1 xmax+.1 -N 0])
   ts = sprintf('Whirling Mode %d of the %d-Chain',k,N);
   title(ts)
   ts = sprintf('omega =%6.2f',omega(1));
   th = text(0,-0.96*N,ts);
   set(th,'EraseMode','xor');
   set(th,'HorizontalAlignment','center');
   for replay = 1:1
      for i = 1:m,
         set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
         if yn == 'y'
            set(hm,'XData',x(i,:),'YData',y(i,:));
         end
         ts = sprintf('omega =%6.2f',omega(i));
         set(th,'String',ts);
         pause(0.02)
      end
      for i = m:-1:1,
         set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
         if yn == 'y'
            set(hm,'XData',x(i,:),'YData',y(i,:));
         end
         ts = sprintf('omega =%6.2f',omega(i));
         set(th,'String',ts);
         pause(0.02)
      end
   end
end
