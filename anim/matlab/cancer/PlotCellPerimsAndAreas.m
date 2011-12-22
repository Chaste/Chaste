% PlotCellPerimsAndAreas.m
% It reads in the positions of all cells at the beginning and end of a 
% Meineke-stlye labelling experiment and plots the percentages of cells 
% that are labelled in ranges 
%

%  close all
clear

% Experiment Setup
title_string = 'Simple Wnt Cells in Sunter i) Geometry';
crypt_height = 30;
%  path = '/tmp/pmxaw/testoutput/Noddy_WNT_No_Area_No_Length/results_from_time_500/';
path1 = '/tmp/pmxaw/testoutput/Noddy_WNT_Yes_Area_Yes_Length/';
colour1 = 'r-';
path2 = '/tmp/pmxaw/testoutput/Noddy_WNT_Yes_Area_No_Length/';
colour2 = 'b--';
path3 = '/tmp/pmxaw/testoutput/Noddy_WNT_No_Area_Yes_Length/';
colour3 = 'g.-';
path4 = '/tmp/pmxaw/testoutput/Noddy_WNT_No_Area_No_Length/';
colour4 = 'k.';
% End of setup

path = path3;
colour = colour3;

addpath('../');	% Adds the LoadNonConstantLengthData function.

Voronoi_data = LoadNonConstantLengthData([path 'All_results']);

X = [];
Y = [];
Area = [];
Perimeters = [];


buckets = 0:1:24;
for i = 1:length(buckets)-1
	Area_in_bucket{i} = [];
	Perim_in_bucket{i} = [];
	
	temp_Area_in_bucket{i} = [];
	temp_Perim_in_bucket{i} = [];
end

for i=1:length(Voronoi_data) % time loop
	num_cells = (length(Voronoi_data{i})-1)/5;
	
	for j = 1:num_cells
		Y = Voronoi_data{i}(5*j-1);
		for k = 1:length(buckets)-1
			if Y >= buckets(k) && Y < buckets(k+1)
				temp_Area_in_bucket{k} = [temp_Area_in_bucket{k} Voronoi_data{i}(5*j)];
				temp_Perim_in_bucket{k} = [temp_Perim_in_bucket{k} (Voronoi_data{i}(5*j+1))^2/Voronoi_data{i}(5*j)]; % shape index perim^2/ area
				break
			end
		end
	end
	if mod(i,50) == 0
		disp(i)
	end
	
	for k = 1:length(buckets)-1
		Area_in_bucket{k} = [Area_in_bucket{k} temp_Area_in_bucket{k}];
		Perim_in_bucket{k} = [Perim_in_bucket{k} temp_Perim_in_bucket{k}];
		
		temp_Area_in_bucket{k} = [];
		temp_Perim_in_bucket{k} = [];
	end
end

for j = 1:length(Area_in_bucket)
	if length(Area_in_bucket{j}) ~= 0
		average_area(j) = mean(Area_in_bucket{j});
		average_perim(j) = mean(Perim_in_bucket{j});
	else
		average_area(j) = 0;
		average_perim(j) = 0;
	end
end

subplot(3,1,1)
hold on
plot(buckets(1:end-2)+0.5*(buckets(2) - buckets(1)), average_area(1:end-1),colour);
subplot(3,1,2)
hold on
plot(buckets(1:end-2)+0.5*(buckets(2) - buckets(1)), 1-average_area(1:end-1),colour);
subplot(3,1,3)
hold on
plot(buckets(1:end-2)+0.5*(buckets(2) - buckets(1)), average_perim(1:end-1),colour);

%  figure;
%  BoxPlotChaste(Area);
%  figure;
%  BoxPlotChaste(Perimeters);

%  buckets = 0:1:24;
%  for i = 1:length(buckets)-1
%  	Area_in_bucket{i} = [];
%  	Perim_in_bucket{i} = [];
%  end
%  
%  for i = 1:length(X)
%  	for k = 1:length(buckets)-1
%  		if Y(i) >= buckets(k) && Y(i) < buckets(k+1)
%  			Area_in_bucket{k} = [Area_in_bucket{k} Area(i)];
%  			Perim_in_bucket{k} = [Perim_in_bucket{k} Perimeters(i)];
%  			break
%  		end
%  	end
%  end
%  
%  for j = 1:length(Area_in_bucket)
%  	if length(Area_in_bucket{j}) ~= 0
%  		average_area(j) = mean(Area_in_bucket{j});
%  		average_perim(j) = mean(Perim_in_bucket{j});
%  	else
%  		average_area(j) = 0;
%  		average_perim(j) = 0;
%  	end
%  end
%  
%  figure
%  plot(buckets(1:end-2)+0.5*(buckets(2) - buckets(1)), average_area(1:end-1),'r');
%  figure;
%  plot(buckets(1:end-2)+0.5*(buckets(2) - buckets(1)), average_perim(1:end-1),'b');