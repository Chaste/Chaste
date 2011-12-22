function [vals, headers]=LoadFileWithHeader(filename)
%
% LOADFILEWITHHEADER    Load a CHASTE-type output file with one header line 
%   followed by n columns of data
% 
%   [DATA, HEADERS] = LOADFILEWITHHEADER(FILENAME) read a file and loads the 
%   data into the first return parameter and optionally the header into the 
%   second return param, as a vector of strings
%
%   Usage: 
%     [data, headers] = LoadFileWithHeader('data.txt')
%     PlotCols(1,3,data,headers)    % plots col 3 against col 1 & labels axes
%
%   See also PLOTCOLS
%


% filename for a new, temporary file
tempfile_name = [filename, '.matlabtemp'];

% call system command copying the file to a new file 
% without the first line
command = ['sed "1d" ', filename, ' > ' , tempfile_name];
system(command);

% load new file
vals = load(tempfile_name);

% remove temporary file
command = ['rm -f ', tempfile_name];
system(command);

if(nargout>1)
   % copy the header to a new file
   tempfile_name = [filename, '.matlabtemp2'];
   command = ['head -n 1 ', filename, ' > ' , tempfile_name];
   system(command);

   % read the labels
   fid = fopen(tempfile_name);
   headers_struct=textscan(fid,'%s');
   headers = headers_struct{1};
   fclose(fid);

   % remove temporary file
   command = ['rm -f ' , tempfile_name];
   system(command);
end;
