N = 7;

curve = ReadNChainSolution('autodata/s.n07phiM001');
omega1 = curve(:,2);
phi1 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM002');
omega2 = curve(:,2);
phi2 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM003');
omega3 = curve(:,2);
phi3 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM004');
omega4 = curve(:,2);
phi4 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM005');
omega5 = curve(:,2);
phi5 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM006');
omega6 = curve(:,2);
phi6 = curve(:,[1 3:(N+1)]);

curve = ReadNChainSolution('autodata/s.n07phiM007');
omega7 = curve(:,2);
phi7 = curve(:,[1 3:(N+1)]);

figure(1)
clf
hold on
plot(omega1,phi1(:,1),'k','Linewidth',2)
plot(omega2,phi2(:,1),'k','Linewidth',2)
plot(omega3,phi3(:,1),'k','Linewidth',2)
plot(omega4,phi4(:,1),'k','Linewidth',2)
plot(omega5,phi5(:,1),'k','Linewidth',2)
plot(omega6,phi6(:,1),'k','Linewidth',2)
plot(omega7,phi7(:,1),'k','Linewidth',2)
xlabel('\omega')
ylabel('\phi_1')
axis([0 6 0 1.75])
plot([0 30],[pi/2 pi/2],'k--','Linewidth',2)
plot([0 30],[pi/6 pi/6],'k--','Linewidth',2)

figure(2)
clf
hold on
plot(omega1,sin(phi1(:,1)),'k','Linewidth',2)
plot(omega2,sin(phi2(:,1)),'k','Linewidth',2)
plot(omega3,sin(phi3(:,1)),'k','Linewidth',2)
plot(omega4,sin(phi4(:,1)),'k','Linewidth',2)
plot(omega5,sin(phi5(:,1)),'k','Linewidth',2)
plot(omega6,sin(phi6(:,1)),'k','Linewidth',2)
plot(omega7,sin(phi7(:,1)),'k','Linewidth',2)

xlabel('\omega','Fontsize',12)
ylabel('sin(\phi_1)','Fontsize',12)
axis([0 6 0 1.05])
plot([0 30],[1 1],'k--','Linewidth',2)
plot([0 30],[0.5 0.5],'k--','Linewidth',2)
set(gca,'YTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1])
set(gca,'Fontsize',12)
title('7-Chain Bifurcation Diagram')
