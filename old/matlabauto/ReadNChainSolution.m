function curve = ReadNChainSolution(filename)
   fid = fopen(filename,'rt');
   [header,count] = fscanf(fid,'%d',12);
   curve = [];
   while (count == 12),
       br = header(1);
       lbl = header(4);
       numdata = header(8);
       numpars = header(12);
       [data,count] = fscanf(fid,'%f',numdata);
       [pars,count] = fscanf(fid,'%f',numpars);
       data(1) = pars(1);
       curve = [curve; data'];
       [header,count] = fscanf(fid,'%d',12);
   end
