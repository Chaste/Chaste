function ret = LoadNonConstantLengthData(filename)
%
% LOADNONCONSTANTLENGTHDATA
%
%   Load a numeric text file with a varying number of entries
%   on each row. For example, if the file 'data.txt' is
%
%   1 2 3
%   4
%   5 6 7 8
%
%   x = LoadNonConstantLengthData('data.txt')
%
%   has x{1} = [1 2 3], x{2} = [4], x{3} = [5 6 7 8].
%
 
fid = fopen(filename);

if(fid<0)
   error('Unable to open file');
end;

i=1;
while 1
   line = fgetl(fid);
   if ~ischar(line), break, end
   ret{i} = str2num(line);
   i = i+1;
end;

fclose(fid);