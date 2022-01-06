function plot_chaste_triangle_mesh(mesh)
% 
% Simple code to read in a triangle mesh and pot it in matlab, useful for debugging without showme 
%
% NOT TESTED OR ROBUST

pos = read_chaste_node_file([mesh,'.node']);

ele = read_chaste_ele_file([mesh,'.ele'])+1;

trimesh(ele,pos(:,1),pos(:,2))
axis([min(pos(:,1))-1,max(pos(:,1))+1,min(pos(:,2))-1,max(pos(:,2))+1])