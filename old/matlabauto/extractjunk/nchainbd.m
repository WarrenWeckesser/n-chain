%
%
ylbls = {'L_2 Norm','\phi_1','\phi_2','\phi_3','\phi_4','\phi_5'};


yn = input('Reread the data? [n] ','s');
if (~isempty(yn) & yn == 'y'),
   fnroot = input('Enter root data file name: ','s');

   fn = sprintf('%s.dat',fnroot);
   fid = fopen(fn,'r');
   disp(['Reading ' fn]);
   dat = fscanf(fid,'%f',[7,inf]);
   fclose(fid);

   fn = sprintf('%s.blk',fnroot);
   fid = fopen(fn,'r');
   disp(['Reading ' fn]);
   blksize = fscanf(fid,'%f',[inf]);
   fclose(fid);
   numblks = length(blksize);

   fid = fopen([fnroot '.lbl'],'r');
   lbl = fscanf(fid,'%f',[9, inf]);
   fclose(fid);
end

ind = input('Enter index of field to plot: (1-6) ');
fld = ind+1;
figure(1);
clf;
hold on;
blkstart = 1;
for k = 1:numblks,
   blkend = blkstart+blksize(k)-1;
   plot(dat(1,blkstart:blkend),dat(fld,blkstart:blkend),'k-');
   blkstart = blkend+1;
end

yn = input('Include solution labels? [n] ','s');
if (~isempty(yn) & yn == 'y')
   numlbl = size(lbl,2);
   for k = 1:numlbl,
      text(lbl(3,k),lbl(fld+2,k),sprintf('%d',lbl(2,k)));
   end
end

xlbl = '\omega';
ylbl = ylbls(ind);

xlabel(xlbl,'FontSize',15,'FontWeight','bold')
ylabel(ylbl,'FontSize',15,'FontWeight','bold')
