# Created by Chaste version 3.0.15620 (with uncommitted modifications) on Mon, 21 May 2012 10:16:28 +0000.  Chaste was built on Mon, 21 May 2012 10:16:14 +0000 by machine (uname) 'Linux compphys11 3.0.0-19-generic #33-Ubuntu SMP Thu Apr 19 19:05:14 UTC 2012 x86_64' using settings: default, no Chaste libraries.
# Read the mesh 
gfx read node SimulationResults.exnode 
gfx read elem SimulationResults.exelem 
gfx define faces egroup SimulationResults
# Create a window 
gfx cre win 1 
# Modify the scene (obtained by gfx list g_element XXXX commands) to visualize first var on lines and nodes 
gfx modify g_element SimulationResults general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none; 
gfx modify g_element SimulationResults lines select_on material default data V spectrum default selected_material default_selected; 
gfx modify g_element SimulationResults node_points glyph point general size "1*1*1" centre 0,0,0 font default select_on material default data V spectrum default selected_material default_selected; 
# Load the data 
for ($i = 0; $i<9; $i++) { 
    gfx read node many_variables_$i.exnode time $i
}
