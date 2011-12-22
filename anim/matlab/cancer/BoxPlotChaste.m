function hout=BoxPlotChaste(x,varargin)
%BOXPLOT Display boxplots of a data sample.
%   BOXPLOT(X) produces a box and whisker plot with one box for each column
%   of X.  The boxes have lines at the lower quartile, median, and upper
%   quartile values.  The whiskers are lines extending from each end of the
%   boxes to show the extent of the rest of the data.  Outliers are data
%   with values beyond the ends of the whiskers.
%
%   BOXPLOT(X,G) produces a box and whisker plot for the vector X grouped
%   by G.  G is a grouping variable defined as a categorical variable,
%   vector, string matrix, or cell array of strings.  G can also be a cell
%   array of several grouping variables (such as {G1 G2 G3}) to group the
%   values in X by each unique combination of grouping variable values.
%
%   BOXPLOT(...,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%      'notch'       'on' to include notches (default is 'off').
%      'symbol'      Symbol and color to use for all outliers (default is 'r+').
%      'orientation' Box orientation, 'vertical' (default) or 'horizontal'.
%      'whisker'     Maximum whisker length (default 1.5).
%      'labels'      Character array or cell array of strings containing
%                    labels for each column of X, or each group in G.
%      'colors'      A string or a three-column matrix of box colors.  Each
%                    box (outline, median line, and whiskers) is drawn in the
%                    corresponding color.  Default is to draw all boxes with
%                    blue outline, red median, and black whiskers.  Colors are
%                    recycled if necessary.
%      'widths'      A numeric vector or scalar of box widths.  Default is
%                    0.5, or slightly smaller for fewer than three boxes.
%                    Widths are recycled if necessary.
%      'grouporder'  When G is given, a character array or cell array of
%                    group names, in the order in which the groups in G
%                    should be plotted.  Ignored when G is not given.
%      'positions'   A numeric vector of box positions.  Default is 1:n.
%
%   In a notched box plot the notches represent a robust estimate of the
%   uncertainty about the medians for box-to-box comparison.  Boxes whose
%   notches do not overlap indicate that the medians of the two groups
%   differ at the 5% significance level.  Whiskers extend from the box
%   out to the most extreme data value within WHIS*IQR, where WHIS is the
%   value of the 'whisker' parameter and IQR is the interquartile range
%   of the sample.
%
%   BOXPLOT(AX,...) plots into the axes with handle AX.
%
%   H = BOXPLOT(...) returns the handle H to the lines in the box plot.
%   H has one column per box, consisting of the handles for the various
%   parts of the box.  Each column contains 7 handles for the upper
%   whisker, lower whisker, upper adjacent value, lower adjacent value,
%   box, median, and outliers.
%
%   Example:  Box plot of car gas mileage grouped by country
%      load carsmall
%      boxplot(MPG, Origin)
%      boxplot(MPG, Origin, 'sym','r*', 'colors',hsv(7))
%      boxplot(MPG, Origin, 'grouporder', ...
%                   {'France' 'Germany' 'Italy' 'Japan' 'Sweden' 'USA'})
%
%   Example: Plot by median gas mileage
%      [sortedMPG,sortedOrder] = sort(grpstats(MPG,Origin,@median));
%      pos(sortedOrder) = 1:6;
%      boxplot(MPG, Origin, 'position', pos)
%
%   See also ANOVA1, KRUSKALWALLIS, MULTCOMPARE.

%   Older syntax still supported:
%       BOXPLOT(X,NOTCH,SYM,VERT,WHIS)
% 
%   BOXPLOT calls BOXUTIL to do the actual plotting.

%   References
%   [1] McGill, R., Tukey, J.W., and Larsen, W.A. (1978) "Variations of
%       Boxplots", The American Statistician, 32:12-16.
%   [2] Velleman, P.F. and Hoaglin, D.C. (1981), Applications, Basics, and
%       Computing of Exploratory Data Analysis, Duxbury Press.
%   [3] Nelson, L.S., (1989) "Evaluating Overlapping Confidence
%       Intervals", Journal of Quality Technology, 21:140-141.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 2.15.4.17 $  $Date: 2006/10/02 16:33:58 $

whissw = 0; % don't plot whisker inside the box.

if isscalar(x) && ishandle(x)
    ax = x;
    x = varargin{1};
    varargin(1) = [];
else
    ax = [];
end
error(nargchk(1+length(ax),Inf,nargin));

if isvector(x)
   % Might have one box, or might have a grouping variable. n will be properly
   % set later for the latter.
   x = x(:);
   n = 1; %
else
   % Have a data matrix, use as many boxes as columns.
   n = size(x,2);
end

% Detect if there is a grouping variable by looking at the second input
nargs = nargin - length(ax);
if nargs < 2
   g = [];
else
   g = varargin{1};
   if isempty(g) || isequal(g,1) || isequal(g,0) || (ischar(g) && size(g,1)==1)
      % It's a NOTCH value or a parameter name
      g = [];
   else
      % It's a grouping variable
      if ~isvector(x)
         error('stats:boxplot:VectorRequired',...
               'X must be a vector when there is a grouping variable.');
      end
      varargin(1) = [];
      nargs = nargs - 1;
   end
end

% Set defaults
notch  = 0;
sym    = 'r+';
vert   = 1;
whis   = 1.5;
labels = {};
colors = []; % default is blue box, red median, black whiskers
posns  = []; % default is 1:n
widths = []; % default is 0.5, smaller for three or fewer boxes
grporder = []; % default is 1:n

% Determine if we have parameter names or the old syntax
if nargs > 1
   if ischar(varargin{1})
      okargs =   {'notch' 'symbol' 'orientation' 'whisker' 'labels' 'colors' 'positions' 'widths' 'grouporder'};
      defaults = { notch   sym      vert          whis      labels   colors   posns       widths   grporder};
      [eid,emsg,notch,sym,vert,whis,labels,colors,posns,widths,grporder] = ...
                                     statgetargs(okargs,defaults,varargin{:});
      if ~isempty(eid)
         error(sprintf('stats:boxplot:%s',eid),emsg);
      end
   else
      if (nargs>=2) && ~isempty(varargin{1}), notch = varargin{1}; end
      if (nargs>=3) && ~isempty(varargin{2}), sym   = varargin{2}; end
      if (nargs>=4) && ~isempty(varargin{3}), vert  = varargin{3}; end
      if (nargs>=5) && ~isempty(varargin{4}), whis  = varargin{4}; end
   end
end

% Convert wordy inputs to internal codes

if isequal(notch,'on')
   notch = 1;
elseif isempty(notch) || isequal(notch,'off')
   notch = 0;
elseif ~isscalar(notch) || ~ismember(notch,0:1)
   error('stats:boxplot:InvalidNotch','Invalid value for ''notch'' parameter');
end

if isempty(vert)
   vert = 1;
elseif ischar(vert)
   vert = strmatch(vert,{'horizontal' 'vertical'}) - 1;
end
if isempty(vert) || ~isscalar(vert) || ~ismember(vert,0:1)   
   error('stats:boxplot:InvalidOrientation',...
         'Invalid value for ''orientation'' parameter');
end

if ~isscalar(whis) || ~isnumeric(whis)
   error('stats:boxplot:BadWhisker',...
         'The ''whisker'' parameter value must be a numeric scalar.');
end

% Deal with grouping variable before processing more inputs
if ~isempty(g)
   if vert, sep = '\n'; else, sep = ','; end
   [g,glabel,gname,multiline] = mgrp2idx(g,size(x,1),sep);
   n = size(gname,1);
   if numel(g) ~= numel(x)
      error('stats:boxplot:InputSizeMismatch',...
            'X and G must have the same length.');
   end
else
    multiline = false;
end

% Reorder the groups if necessary
if ~isempty(g) && ~isempty(grporder)
   if iscellstr(grporder) || ischar(grporder)
      % If we have a grouping vector, grporder may be a list of group names.
      if ischar(grporder), grporder = cellstr(grporder); end
      [dum,grporder] = ismember(glabel,grporder(:));
      % Must be a permutation of the group names
      if ~isequal(sort(grporder),(1:n)')
         error('stats:boxplot:BadOrder', ...
               'The ''grouporder'' parameter value must contain all the unique group names in G.');
      end
   else
      error('stats:boxplot:BadOrder', ...
            'The ''grouporder'' parameter value must be a character array or a cell array of strings.');
   end
   g = grporder(g);
   glabel(grporder) = glabel;
   gname(grporder,:) = gname;
end

% Process the rest of the inputs

if isempty(labels)
   if ~isempty(g)
      labels = glabel;
   end
else
   if ~(iscellstr(labels) && numel(labels)==n) && ...
      ~(ischar(labels) && size(labels,1)==n)
      % Must have one label for each box
      error('stats:boxplot:BadLabels','Incorrect number of box labels.');
   end
   if ischar(labels), labels = cellstr(labels); end
   multiline = false;
end
dfltLabs = (isempty(labels) && isempty(g)); % box labels are just column numbers

if isempty(widths)
   widths = repmat(min(0.15*n,0.5),n,1);
elseif ~isvector(widths) || ~isnumeric(widths) || any(widths<=0)
   error('stats:boxplot:BadWidths', ...
         'The ''widths'' parameter value must be a numeric vector of positive values.');
elseif length(widths) < n
   % Recycle the widths if necessary.
   widths = repmat(widths(:),ceil(n/length(widths)),1);
end

if isempty(colors)
   % Empty colors tells boxutil to use defaults.
   colors = char(zeros(n,0));
elseif ischar(colors) && isvector(colors)
   colors = colors(:); % color spec string, make it a column
elseif isnumeric(colors) && (ndims(colors)==2) && (size(colors,2)==3)
   % RGB matrix, that's ok
else
   error('stats:boxplot:BadColors',...
         'The ''colors'' parameter value must be a string or a three-column numeric matrix.');
end
if size(colors,1) < n
   % Recycle the colors if necessary.
   colors = repmat(colors,ceil(n/size(colors,1)),1);
end

if isempty(posns)
   posns = 1:n;
elseif ~isvector(posns) || ~isnumeric(posns)
   error('stats:boxplot:BadPositions', ...
         'The ''positions'' parameter value must be a numeric vector.');
elseif length(posns) ~= n
   % Must have one position for each box
   error('stats:boxplot:BadPositions', ...
         'The ''positions'' parameter value must have one element for each box.');
else
   [dum,ord] = sort(posns);
   if isempty(labels) % labels never empty when grouping vector supplied
       % If we have matrix data with no labels, and the positions are not 1:n,
       % we need to force the default column number tick labels 1:n, but in
       % the right positions.
       labels = cellstr(num2str(ord(:)));
   else
       % Permute the labels to match the plot position order.
       labels = labels(ord);
   end
end

%
% Done processing inputs
%

% Put at least the widest box or half narrowest spacing at each margin
if n > 1
    wmax = max(max(widths), 0.5*min(diff(posns)));
else
    wmax = 0.5;
end
xlims = [min(posns)-wmax, max(posns)+wmax];

ymin = min(x(:));
ymax = max(x(:));
if ymax > ymin
   dy = (ymax-ymin)/20;
else
   dy = 0.5;  % no data range, just use a y axis range of 1
end
ylims = [(ymin-dy) (ymax+dy)];

% Scale axis for vertical or horizontal boxes.
ax = newplot(ax);
axes(ax);
oldstate = get(ax,'NextPlot');
set(ax,'NextPlot','add','Box','on');
if isempty(xlims)
    xlims = [0 1];
end
if isempty(ylims)
    ylims = [0 1];
end
if vert
    axis(ax,[xlims ylims]);
    set(ax,'XTick',sort(posns),'Units','normalized');
    drawnow;
    setappdata(ax,'NormalizedOuterPosition',get(ax,'OuterPosition'));
    ylabel(ax,'Values');
    if dfltLabs, xlabel(ax, 'Column Number'); end
else
    axis(ax,[ylims xlims]);
    set(ax,'YTick',sort(posns));
    xlabel(ax,'Values');
    if dfltLabs, ylabel(ax,'Column Number'); end
end
if nargout>0
   hout = [];
end

xvisible = NaN(size(x));
notnans = ~isnan(x);
for i= 1:n
   if ~isempty(g)
      thisgrp = find((g==i) & notnans);
   else
      thisgrp = find(notnans(:,i)) + (i-1)*size(x,1);
   end
   [outliers,hh] = boxutil(ax,x(thisgrp),notch,posns(i),widths(i), ...
                           colors(i,:),sym,vert,whis,whissw);
   outliers = thisgrp(outliers);
   xvisible(outliers) = x(outliers);
   if nargout>0
      hout = [hout, hh(:)];
   end
end

if ~isempty(labels)
   if multiline && vert
      % Turn off tick labels and axis label
      set(ax, 'XTickLabel','');
      setappdata(ax,'NLines',size(gname,2));
      xlabel(ax,'');
      ylim = get(ax, 'YLim');
      
      % Place multi-line text approximately where tick labels belong
      ypos = repmat(ylim(1),size(posns));
      text(posns,ypos,labels,'HorizontalAlignment','center', ...
                             'VerticalAlignment','top', 'UserData','xtick');

      % Resize function will position text more accurately
      f = ancestor(ax,'figure');
      set(f, 'ResizeFcn', @resizefcn, ...
               'Interruptible','off', 'PaperPositionMode','auto');
      resizefcn(f);
   elseif vert
      set(ax, 'XTickLabel',labels);
   else
      set(ax, 'YTickLabel',labels);
   end
end
set(ax,'NextPlot',oldstate);

% Store information for gname function
set(ax, 'UserData', {'boxplot' xvisible g vert});


%=============================================================================

function [outlier,hout] = boxutil(ax,x,notch,lb,lf,clr,sym,vert,whis,whissw)
%BOXUTIL Produces a single box plot.

% define the median and the quantiles
pctiles = prctile(x,[25;50;75]);
q1 = pctiles(1,:);
med = pctiles(2,:);
q3 = pctiles(3,:);

% find the extreme values (to determine where whiskers appear)
vhi = q3+whis*(q3-q1);
upadj = max(x(x<=vhi));
if (isempty(upadj)), upadj = q3; end

vlo = q1-whis*(q3-q1);
loadj = min(x(x>=vlo));
if (isempty(loadj)), loadj = q1; end

x1 = repmat(lb,1,2);
x2 = x1+[-0.25*lf,0.25*lf];
outlier = x<loadj | x > upadj;
yy = x(outlier);

xx = repmat(lb,1,length(yy));
lbp = lb + 0.5*lf;
lbm = lb - 0.5*lf;

if whissw == 0
   upadj = max(upadj,q3);
   loadj = min(loadj,q1);
end

% Set up (X,Y) data for notches if desired.
if ~notch
    xx2 = [lbm lbp lbp lbm lbm];
    yy2 = [q3 q3 q1 q1 q3];
    xx3 = [lbm lbp];
else
    n1 = med + 1.57*(q3-q1)/sqrt(length(x));
    n2 = med - 1.57*(q3-q1)/sqrt(length(x));
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
    lnm = lb-0.25*lf;
    lnp = lb+0.25*lf;
    xx2 = [lnm lbm lbm lbp lbp lnp lbp lbp lbm lbm lnm];
    yy2 = [med n1 q3 q3 n1 med n2 q1 q1 n2 med];
    xx3 = [lnm lnp];
end
yy3 = [med med];

    
% Determine if the boxes are vertical or horizontal.
% The difference is the choice of x and y in the plot command.
if vert
    hout = plot(ax,x1,[q3 upadj],'k--', x1,[loadj q1],'k--',...
                x2,[upadj upadj],'k-', x2,[loadj loadj],'k-', ...
                xx2,yy2,'b-', xx3,yy3,'r-', xx,yy,sym);
else
    hout = plot(ax,[q3 upadj],x1,'k--', [loadj q1],x1,'k--',...
                [upadj upadj],x2,'k-', [loadj loadj],x2,'k-', ...
                yy2,xx2,'b-', yy3,xx3,'r-', yy,xx,sym);
end

% If there's a color given, show everything in that color.  If the outlier
% symbol has a color in it, leave those alone.
if ~isempty(clr)
   if any(ismember(sym,'bgrcmykw'))
      set(hout(1:6),'Color',clr);
   else
      set(hout,'Color',clr);
   end
end

set(hout(1),'Tag','Upper Whisker');
set(hout(2),'Tag','Lower Whisker');
set(hout(3),'Tag','Upper Adjacent Value');
set(hout(4),'Tag','Lower Adjacent Value');
set(hout(5),'Tag','Box');
set(hout(6),'Tag','Median');
if length(hout)>=7
   set(hout(7),'Tag','Outliers');
else
   hout(7) = NaN;
end


%=============================================================================

function resizefcn(f,dum)

% Adjust figure layout to make sure labels remain visible
h = findobj(f, 'UserData','xtick');
if (isempty(h))
   set(f, 'ResizeFcn', '');
   return;
end

% Loop over all axes
allax = findall(f,'Type','Axes');
for j=1:length(allax)
    ax = allax(j);
    nlines = getappdata(ax, 'NLines');
    
    if ~isempty(nlines)
        % Try to retain the original normalized outer position
        nop = getappdata(ax,'NormalizedOuterPosition');
        set(ax,'Units','normalized');
        if isempty(nop)
            nop = get(ax,'OuterPosition');
        end
        set(ax,'OuterPosition',nop);
        
        % Adjust position so the fake X tick labels have room to display
        temp = hgconvertunits(f,[0 0 1 1],'character','normalized',f);
        charheight = temp(4);
        li = get(ax,'LooseInset');
        tickheight = min((nlines+1)*charheight, nop(4)/2);
        p = [nop(1)+li(1)*nop(3),    nop(2)+tickheight, ...
             nop(3)*(1-li(1)-li(3)), nop(4)*(1-li(4))-tickheight];
        p = max(p,0.0001);
        set(ax, 'Position', p);
        
        % The following lines do not change the position, but they leave
        % the axes in a state such that MATLAB will try to preserve the
        % outerposition rather than the ordinary position
        op = get(ax,'OuterPosition');
        set(ax, 'OuterPosition', op);
    end
end


% Position the labels at the proper place
xl = get(ax, 'XLabel');
set(xl, 'Units', 'data');
p = get(xl, 'Position');
ylim = get(ax, 'YLim');
p2 = (p(2)+ylim(1))/2;
for j=1:length(h)
   p = get(h(j), 'Position') ;
   p(2) = p2;
   set(h(j), 'Position', p);
end
