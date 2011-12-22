# Created by Chaste version 2.1.10890 on Mon, 22 Nov 2010 05:09:03 +0000.  Chaste was built on Mon, 22 Nov 2010 05:08:50 +0000 by machine (uname) 'Linux compphys11 2.6.32-25-generic #45-Ubuntu SMP Sat Oct 16 19:52:42 UTC 2010 x86_64' using settings: default_2, no Chaste libraries.
# Read the mesh 
gfx read node SimulationResults.exnode 
gfx read elem SimulationResults.exelem generate_faces_and_lines 
# Create a window 
gfx cre win 1 
# Modify the scene (obtained by gfx list g_element XXXX commands) to visualize first var on lines and nodes 
gfx modify g_element SimulationResults general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none; 
gfx modify g_element SimulationResults lines select_on material default data V spectrum default selected_material default_selected; 
gfx modify g_element SimulationResults node_points glyph point general size "1*1*1" centre 0,0,0 font default select_on material default data V spectrum default selected_material default_selected; 
# Load the data 
for ($i=0; $i<2; $i++) { 
    gfx read node cube_2mm_12_elements_$i.exnode time $i
}
