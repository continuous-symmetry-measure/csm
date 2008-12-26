%    input     %
%--------------

clear all;

fid1 = fopen(['/result1/mark/danovich/Density-Cu7-d5h-dist1.bin'],'r','ieee-le');

fid2 = fopen(['Den-Cu7-d5h-dist1.bin'],'w','native');
fid3 = fopen(['Cu7-d5h-dist1.bin'],'w','native');

data = fread(fid1,[11],'int32');
fwrite(fid3,data,'int32');
fclose(fid3);

nn = 500;
for i=1:nn,
   for j=1:nn,
      data = fread(fid1,[nn],'int32');
      fwrite(fid2,data,'int32');
   end
end

fclose(fid2);
fclose(fid1);
