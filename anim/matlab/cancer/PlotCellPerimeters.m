function PlotCellPerimeters(y0, y1)
%
% PLOTCELLPERIMETERS
% 
%   PlotCellPerimeters(y0,y1) loads 'VoronoiAreaAndPerimeter.dat'
%   and plots the average perimeter for each cell who's y value 
%   is between y0 and y1.
 
data = LoadNonConstantLengthData('VoronoiAreaAndPerimeter.dat');

times = [];
average_perims = [];
for i=1:length(data)
    num_nodes = (length(data{i})-1)/5;
    
    av_perim = 0;
    num_cells_in_range = 0;
    for j=1:num_nodes
        y = data{i}(1+5*j-2);
        if(y>=y0 & y<=y1)
           av_perim = av_perim + data{i}(1+5*j);
           num_cells_in_range = num_cells_in_range+1;
        end;
    end;

    if(num_cells_in_range>0)
        av_perim = av_perim/num_cells_in_range;

        times = [times; data{i}(1)];
        average_perims = [average_perims; av_perim];
    end;
end;

plot(times, average_perims,'.');