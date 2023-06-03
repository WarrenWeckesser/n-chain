%
%
%
N = input('Enter the number of segments in the chain: ');
B = nmatrix(N);
D = diag([N:-1:1]);
[evecs,mevals] = eig(D,B);
[evals,eind] = sort(diag(mevals)');
disp('The bifurcation values of omega^2 are:');
disp(evals)
%disp('The associated eigenfunctions are:');
%disp(evecs(:,eind));
k = input('Which mode do you want to follow? ');
omega2 = evals(k);
%
% The k-th eigenvector gives the initial search direction.
%
evec = evecs(:,eind(k));
dir = [evec' 0];
% (dir is a row)

delta = 0.1;
clear phi
%
% Starting values:
%
phi(1,:) = zeros(1,N+1);
phi(1,N+1) = sqrt(sqrt(omega2));
%
m = 75;
%
solve_opts = optimset('TolFun',1e-10,'TolX',1e-10,'Diagnostics','off','Display','off');

for i = 2:m,
   dir = dir/norm(dir);
   phiguess = phi(i-1,:)+delta*dir;
   phinext = fsolve('nlevps',phiguess',solve_opts,D,B);
   dstr = sprintf('%3d phinext(1)=%12.8f phinext(N)=%12.8f phinext(N+1)=%12.8f',i,phinext(1),phinext(N),phinext(N+1));
   disp(dstr);
   phi(i,:) = sign(phinext(1))*phinext';
   phi(i,N+1) = abs(phi(i,N+1));
   dir = (phi(i,:)-phi(i-1,:));
end

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

dfn = sprintf('N%03dM%03d.eps',N,k);
disp('Enter n in the next line to not print the plot to a file.');
promptstr = sprintf('File name for plot (default is %s): ',dfn);
fn = input(promptstr,'s');
if isempty(fn) fn = dfn; end
if (fn ~= 'n')
   disp('Printing to file...');
   estr = sprintf('print  -deps %s',fn);
   eval(estr);
end

yn = input('Make a movie? ','s');
if (yn == 'y')
   figure(2)
   clf
   hold off
   yn = input('Mark the masses with circles in the movie? ','s');
   xmax = max(max(abs(x)));
   hs = plot([0 x(1,:)],[0 y(1,:)],'k');
   set(hs,'EraseMode','xor');
   hold on
   if yn == 'y'
      hm = plot(x(1,:),y(1,:),'bo');
      set(hm,'EraseMode','xor');
   end
   axis equal
   axis([-xmax-.1 xmax+.1 -N 0])
   ts = sprintf('Whirling Mode %d of the %d-Chain',k,N);
   title(ts)
   ts = sprintf('omega =%6.2f',phi(1,N+1)^2);
   th = text(0,-0.96*N,ts);
   set(th,'EraseMode','xor');
   set(th,'HorizontalAlignment','center');
   for replay = 1:1
      for i = 1:m,
         set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
         if yn == 'y'
            set(hm,'XData',x(i,:),'YData',y(i,:));
         end
         ts = sprintf('omega =%6.2f',phi(i,N+1)^2);
         set(th,'String',ts);
         pause(0.02)
      end
      for i = m:-1:1,
         set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
         if yn == 'y'
            set(hm,'XData',x(i,:),'YData',y(i,:));
         end
         ts = sprintf('omega =%6.2f',phi(i,N+1)^2);
         set(th,'String',ts);
         pause(0.02)
      end
   end
end


yn = input('Find the eigenvalues for this family? ','s');

if yn == 'y'
   % Do stability check
   Z = zeros(size(B));
   clear evs mr
   st = 1;
   for j = 1:m,
      phij = phi(j,1:N);
      tanphij = tan(phij);
      cosphij = cos(phij);
      omega2  = phi(j,N+1);
      omega = sqrt(omega2);
      M = [B + B .* (tanphij'*tanphij) zeros(size(B)); Z B];
      C = 2*sqrt(omega2)*[ Z -B; B Z];
      w1 = [N:-1:1] ./ (cosphij.^3);
      w2 = [N:-1:1] ./ cosphij;
      W = [diag(w1)-omega2*B Z; Z diag(w2)-omega2*B];
      Minv = inv(M);
      A = [zeros(2*N,2*N) eye(2*N,2*N); -Minv*W -Minv*C];
      evs(j,:) = eig(A)';
      mr(j) = max(real(evs(j,:)));
      if st == 1
         if mr > 1e-6
            disp('Stability lost.');
            omega
            st = 0;
         end
      end
   end
   %
   % Replot
   %
   figure(1)
   clf
   hold on
   for i = 1:m,
      plot([0 x(i,:)],[0 y(i,:)],'k')
      sevs = sort(evs(i,:));
      nzevs = sevs(3:length(sevs));
      if max(real(nzevs)) > 1e-5
         plot(x(i,:),y(i,:),'ro')
      else 
         plot(x(i,:),y(i,:),'go')
      end
   end
   xmax = max(max(abs(x)));
   axis equal
   axis([-xmax-.1 xmax+.1 -N 0])
   ts = sprintf('Whirling Mode %d of the %d-Chain',k,N);
   title(ts)

   figure(3)
   clf
   clear sevs
   for i = 1:m,
      [evtmp,sind] = sort(imag(evs(i,:)));
      evtmp = evs(i,sind);
      sevs(i,:) = evtmp(2*N+1:4*N);
   end
   hold on
   for j = 1:2*N,
      plot(imag(sevs(:,j)))
   end
   hold off
   a = axis;
   a(1,1) = 1;
   axis(a);
   ts = sprintf('Imaginary Parts of the Eigenvalues in Mode %d of the %d-Chain',k,N);
   title(ts)

%   for j = 2*N:-1:1,
%      t1 = [imag(sevs(:,j))' 9999];
%      t2 = [-1 imag(sevs(:,j))'];
%      cmp = sum(t2 < t1);
%      if cmp ~= size(sevs,1)+1
%         st = sprintf('First non-monotonic eigenvalue at %d.',j);
%         disp(st);
%         break;
%      end
%   end

   yn = input('Make a movie? ','s');
   if (yn == 'y')
      yn = input('Mark the masses with circles in the movie? ','s');
      figure(4)
      clf
      hold off
      xmax = max(max(abs(x)));
      rmax = max(max(real(evs)));
      imax = max(max(imag(evs)));
      subplot(1,2,1);
      hs = plot([0 x(1,:)],[0 y(1,:)],'k');
      set(hs,'EraseMode','xor');
      hold on
      if yn == 'y'
         hm = plot(x(1,:),y(1,:),'bo')   ;
         set(hm,'EraseMode','xor');
      end
      axis equal
      axis([-xmax-.1 xmax+.1 -N 0])
      ts = sprintf('Whirling Mode %d of the %d-Chain',k,N);
      title(ts)
      ts = sprintf('omega =%6.2f',sqrt(phi(1,N+1)));
      th = text(0,-0.96*N,ts);
      set(th,'EraseMode','xor');
      set(th,'HorizontalAlignment','center');
      subplot(1,2,2);
      he = plot(real(evs(1,:)),imag(evs(1,:)),'mo');
      set(he,'EraseMode','xor');
      axis(1.1*[-rmax rmax 0 imax]);
      ts = sprintf('Eigenvalues in upper half of the complex plane');
      title(ts)
      xlabel('Real part')
      ylabel('Imaginary part')
      for replay = 1:5,
         for i = 1:m,
            set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
            if yn == 'y'
               set(hm,'XData',x(i,:),'YData',y(i,:));
            end
            ts = sprintf('omega =%6.2f',sqrt(phi(i,N+1)));
            set(th,'String',ts);
            set(he,'XData',real(evs(i,:)),'YData',imag(evs(i,:)));
            pause(0.05);
         end
         for i = m:-1:1,
            set(hs,'XData',[0 x(i,:)],'YData',[0 y(i,:)]);
            if yn == 'y'
               set(hm,'XData',x(i,:),'YData',y(i,:));
            end
            ts = sprintf('omega =%6.2f',sqrt(phi(i,N+1)));
            set(th,'String',ts);
            set(he,'XData',real(evs(i,:)),'YData',imag(evs(i,:)));
            pause(0.05);
         end
      end
   end


end


