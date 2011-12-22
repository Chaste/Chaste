function pos = read_chaste_node_file(file)
% READ_CHASTE_NODE_FILE
% 
% Very simple method for reading in a chaste node file 
% in triangles(tetgen) format, which just calls load
% after deleting comments and the first line (the header
% line).
%
% NOT TESTED OR ROBUST

% delete any comment lines
command = ['grep -v ^# ',file,' > /tmp/chaste_temp_node_file'];
system(command);
% remove first line (the header line)
command = ['tail -n +2 /tmp/chaste_temp_node_file > /tmp/chaste_temp_node_file2'];
system(command);
% read the new file, then delete
pos = load('/tmp/chaste_temp_node_file2');
system('rm /tmp/chaste_temp_node_file2');
% get rid of first column (the node indices)
pos = pos(:,2:end);
