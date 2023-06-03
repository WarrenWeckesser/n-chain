   phi2 = phi(1,:);
   k2 = 1;
   k1 = 2;
   while (k1 < m)
      while (k1 < m) & ((norm(phi(k1,1:N)) - norm(phi2(k2,1:N))) < 0.02)
         k1 = k1+1;
      end
      k2 = k2 + 1;
      phi2 = [phi2; phi(k1,:)];
      k1 = k1+1;
   end
size(phi2,1)