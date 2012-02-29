%%% 
% To remove Octave copyight, run with 
%   octave annulus.m | tail -n +19 > temp.poly
% Note that the area quality flag is based on the 
% last line of temp.poly (half-height-base)
%%%

%{ 
#length 0.025
octave annulus.m | tail -n +19 > temp.poly
triangle -e -a0.0003 -q temp.poly
mv temp.1.node ../annulus.node
mv temp.1.edge ../annulus.edge
mv temp.1.ele ../annulus.ele

#length 0.1
octave annulus.m | tail -n +19 > temp.poly
triangle -e -a0.005 -q temp.poly
mv temp.1.node ../annulus_coarse.node
mv temp.1.edge ../annulus_coarse.edge
mv temp.1.ele ../annulus_coarse.ele
%}

length=0.025;
%length=0.1;
rad_out=0.9;
rad_in=0.6;

out_nodes=ceil(2*pi*rad_out/length);
in_nodes=ceil(2*pi*rad_in/length);
all_nodes=out_nodes+in_nodes;
%Write total number of nodes
printf("%i\n", all_nodes)

%Write all outer nodes
for step=0:out_nodes-1
    theta=2*pi*step/out_nodes;
    printf("%i\t%f\t%f\n", step, rad_out*cos(theta), rad_out*sin(theta))
end
%Write all inner nodes
for step=0:in_nodes-1
    theta=-2*pi*step/in_nodes; % Winds clockwise, the other way
    printf("%i\t%f\t%f\n", step+out_nodes, rad_in*cos(theta), rad_in*sin(theta))
end

%Write total edges
printf("%i\t0\n", all_nodes)

%Write all outer edges
for step=0:out_nodes-2
   printf("%i\t%i\t%i\n", step, step, step+1)
end
%Close the loop
printf("%i\t%i\t%i #End of outer loop\n", out_nodes-1, out_nodes-1, 0)

%Write all inner edges
for step=out_nodes:all_nodes-2
   printf("%i\t%i\t%i\n", step, step, step+1)
end
%Close the loop
printf("%i\t%i\t%i #End of inner loop\n", all_nodes-1, all_nodes-1, out_nodes)

%Indicate the hole
printf("1 #One hole\n0 0.0 0.0\n")
printf("#Length is %f\n#Outer radius %f\n#Inner radius %f\n#Area %f\n", length, rad_out, rad_in, length^2/2)
