function PlotCols(x,y,data,labels,plottype)
% PLOTCOLS    Plot the columns of a matrix of data
%    PLOT(X,Y,DATA) plots column Y of data against column X. Y can be a vector if 
%    columns are to be plotted (in which case different colours are used (if 
%    length(Y)<8)).
%
%    PLOT(X,Y,DATA,LABELS), where LABELS is a vectors of strings, uses these strings
%    to label the axes (or create a legend). LABELS should be of size num columns in DATA
%
%    PLOT(X,Y,DATA,LABELS,PLOTTYPES) where PLOTTYPES is a struct of strings can be
%    used to provide to the type of plot, eg 'b-' for blue lines, 'r.' for red dots etc
%    Use LABELS=[] if no labels
%
%    Examples (here 'data' is a n-by-4 matrix, 'labels' a 4 by 1 array of strings):
%       PlotCols(1,2,data);                     % plot column 2 against column 1
%       PlotCols(1,[2,3],data);                 % plot columns 2 and 3 against column 1
%       PlotCols(1,[2,4],data,[],{'b-','k*'});  % plot columns 2 (blue) and 4 (black stars) 
%                                               % against column 1
%       PlotCols(1,[2,3],data,labels);          % labels the axes
%
%    Typical usage:
%       [data,headers] = LoadFileWithHeader('Lr91BackwardEuler.dat');
%       PlotCols(1,[2 4 3 7 8 9],data,headers);  % plots all the gating variables
%
%    See also LOADFILEWITHHEADER
%


if x<=0 || x>size(data,2)
   error('X too large or not positive, should be a index for a column of the matrix DATA');
end;

if nargin > 3 
   if (length(labels)>0) 
      if length(labels)~=size(data,2)
         error('Size of LABELS does not match number of columns in DATA')
      end;
   end;
end;

if nargin > 4
   if (length(plottypes)~=length(y))
      error('Size of plottypes does not match the number of columns to be plotted');
   end;
else
   plottypes = {'b','k','r','g','c','y','m'};
   if(length(y)>7)
      % too many grpahs to be drawm, but user didn't specify plottypes explicity, 
      % just draw them all blue
      for i=8:length(y)
         plottypes = {plottypes{1:end}, 'b'};
      end;
   end;
end;


figure
for i=1:length(y)
   if y(i)<=0 || y(i)>size(data,2)
      error('Y(i) too large or not positive, should be a index for a column of the matrix DATA');
   end;

   if x==y(i)
      error('Plotting a column against itself')
   end;

   hold on
   plot(data(:,x), data(:,y(i)), plottypes{i});
end;

% use the headers, if given, to label the axes (or create a legend if several columns were
% plotted)
if nargin > 3 
   if (length(labels)>0) 
      if (length(y)==1)
         xlabel(labels(x),'fontsize',16);
         ylabel(labels(y(1)),'fontsize',16);
      else
         xlabel(labels(x),'fontsize',16);
         ylabels = [];
         for i=1:length(y)
             ylabels = [ylabels, labels(y(i))];
         end;
         legend(ylabels,'fontsize',16);
      end;
   end;
end;

hold off
