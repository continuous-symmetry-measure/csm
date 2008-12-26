function convert_file(inputFile, descFile, resFile)
% convert_file(inputFile, descFile, resFile)
% Converts the file to a binary file describing the information, 
% and a density file

fid1 = fopen(inputFile,'r','ieee-le');

fid2 = fopen(resFile,'w','native');
fid3 = fopen(descFile,'w','native');

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
