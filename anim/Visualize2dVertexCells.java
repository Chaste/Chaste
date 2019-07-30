/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.StringTokenizer;
import java.lang.Math;
import javax.imageio.ImageIO;

import javax.swing.JPanel;
import javax.swing.JLabel;

public class Visualize2dVertexCells implements ActionListener, AdjustmentListener, ItemListener, Runnable
{
    private Thread updateThread;

    static CustomVertexCanvas2D canvas;

    Button run;
    
    public Frame frame = new Frame();
    
    public static boolean parsed_all_files = false;
    public static boolean drawAncestors = false;
    public static boolean drawPotts = false;
    public static boolean drawAxes = true;
    public static boolean drawCells = true;
    public static boolean drawCircles = false;
    public static boolean writeFiles = false;
    public static boolean drawCylinder = false;
    public static boolean drawCylinderOverride = true;
    public static boolean ancestorsFilePresent = false;
    public static boolean elementFilePresent = true;
    // by default the last timestep isn't read or visualised; this
    // allows the visualiser to be run as a simulation is being run.
    // To visualise the last timestep, use "showlaststep" as an argument
    public static boolean showLastStep = false; 
    public static boolean axesEqual = false;
    
    public static int timeStep = 0;
    public static int delay = 50;
    public static int numTimeSteps = 0;
    public static int memory_factor = 2; // fudge factor
    public static int[] numNodes;
    public static int[] numCells;
    public static int[] numElements;
    public static int[][] ancestor_values;
    public static int[][] element_nodes; // stores the node indices associated with each element
    public static int[][] element_num_nodes; // stores the number of nodes associated with each element
    public static int[][] cell_types;
    public static int[][] image_nodes;
    public static int[][] is_boundary_nodes;
    
    
    public static double max_x = -1e12;
    public static double max_y = -1e12;
    public static double min_x =  1e12;
    public static double min_y =  1e12;
    public static int lattice_sites_up = 0;
    public static int lattice_sites_across = 0;
    public static double crypt_width = 0.0;
    public static double half_width = 0.0;
    public static double[] times;
    
    public static RealPoint[][] node_positions;
    public static RealPoint[][] element_centroid_positions;
     

    public static Scrollbar delay_slider = new Scrollbar(Scrollbar.HORIZONTAL, delay, 1, 1, 100);
    public static Scrollbar time_slider = new Scrollbar(Scrollbar.HORIZONTAL, timeStep, 1, 0, 2);
    
    public static Checkbox output = new Checkbox("Output");
    public static Checkbox cells = new Checkbox("Cells");
    public static Checkbox ancestors = new Checkbox("Clonal Populations");
    public static Checkbox potts = new Checkbox("Potts Simulation");
    public static Checkbox axes = new Checkbox("Axes");
    public static Checkbox axes_equal = new Checkbox("Axes Equal");
    
    public static JLabel nearest_node_label = new JLabel();
    public static JLabel nearest_element_centroid_label = new JLabel();
    
    public static File node_file;
    public static File boundary_node_file;
    public static File cell_type_file;
    public static File element_file;
    public static File ancestors_file;
    public static File setup_file;
    
    public static Button refresh;
    
    public Visualize2dVertexCells() 
    {
        frame.setSize(700, 700);
        frame.setLayout(new BorderLayout());
        
        canvas = new CustomVertexCanvas2D(this);
        canvas.setPreferredSize(new Dimension(frame.getWidth(),frame.getHeight()));
        canvas.addMouseMotionListener(canvas);
        
        addButtons(frame);
        addTimeSlider(frame);
        
        JPanel canvasPanel = new JPanel();
        canvasPanel.add(canvas);
        frame.add(canvasPanel, BorderLayout.CENTER);
                
        frame.addWindowListener(new WindowAdapter() 
        {
            public void windowClosing(WindowEvent e) 
            {
                System.exit(0);
            }
        });
        frame.pack();
        frame.setVisible(true);
    }

    public void actionPerformed(ActionEvent event) 
    {
        String pressed = event.getActionCommand();
                
        if (pressed == "Quit") 
        {
            frame.dispose();
        }
        if (pressed == "Run") 
        {
            if (timeStep == numTimeSteps - 1) 
            {
                timeStep = 0;
                time_slider.setValue(timeStep);
            }
            if (updateThread == null) 
            {
                run.setLabel("Pause");
                updateThread = new Thread(this);
                updateThread.start();
            }
        }
        if (pressed == "Reset") 
        {
            timeStep = 0;
            time_slider.setValue(timeStep);
            canvas.drawBufferedImage();
            canvas.repaint();
        }
        if (pressed == "Pause") 
        {
            if (updateThread != null) 
            {
                Thread fred = updateThread;
                updateThread = null;
                fred.interrupt();
                run.setLabel("Run");
            }
        }
        if (pressed == "Refresh")
        {
            refresh.setEnabled(false);
        	LoadAllFiles();
            canvas.drawBufferedImage();
            canvas.repaint();
            refresh.setEnabled(true);
        }
    }
    
    public void itemStateChanged(ItemEvent e) 
    {
        Object cb = e.getItemSelectable();
        boolean state = (e.getStateChange() == ItemEvent.SELECTED);
               
        if (cb == output) 
        {
            writeFiles = state;
            System.out.println("Writing output files = "+writeFiles);
        }
        else if (cb == cells)
        {
            drawCells = state;
            System.out.println("Drawing cells = "+drawCells);
        }
        else if (cb == axes)
        {
            drawAxes = state;
            System.out.println("Drawing axes = "+drawAxes); 
        }
        else if (cb == axes_equal)
        {
            axesEqual = state;
            
            if (axesEqual==true) // Make the axes equal
            {
	            if (min_y < min_x)
	        	{
	        		min_x = min_y;
	        	}
	            else
	            {
	            	min_y = min_x;
	            }
	            if (max_x < max_y)
	        	{
	        		max_x = max_y;
	        	}
	            else
	            {
	            	max_y = max_x;
	            }
            }
            else // reset the axes
            { 
            	CalculateCanvasDimensions();
            }
            System.out.println("Drawing axes equal = "+axesEqual); 
        }
        else if (cb == ancestors)
        {
            drawAncestors = state;
            System.out.println("Drawing clonal populations = " + drawAncestors); 
        }
        else if (cb == potts)
        {
            drawPotts = state;
            System.out.println("Drawing Potts Simulation = " + drawPotts); 
        }
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    public void adjustmentValueChanged(AdjustmentEvent e) 
    {
        delay = delay_slider.getValue();
        timeStep = time_slider.getValue();
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    public void run() 
    {
        while (updateThread != null) 
        {            
            if (timeStep < numTimeSteps - 1) 
            {
            	timeStep++;
            } 
            else 
            {
                if (updateThread != null) 
                {
                    Thread thread = updateThread;
                    updateThread = null;
                    thread.interrupt();
                    run.setLabel("Run");
                }
            }
            try 
            {
                Thread.sleep((100 - delay) * 1);
            } 
            catch (InterruptedException e) 
            {
                return;
            }
            canvas.drawBufferedImage();
            canvas.repaint();
        }
    }

    public void addButtons(Frame frame) 
    {
        JPanel buttonPanel = new JPanel(new GridLayout(0,4));
        Button quit = new Button("Quit");
        quit.addActionListener(this);

        run = new Button("Run");
        run.addActionListener(this);
        
        Button reset = new Button("Reset");
        reset.addActionListener(this);
        
        refresh = new Button("Refresh");
        refresh.setEnabled(true);
        refresh.addActionListener(this);
        
        buttonPanel.add(quit);
        buttonPanel.add(run);
        buttonPanel.add(reset);
        buttonPanel.add(refresh);
                
        JPanel scrollPanel = new JPanel();
        
        delay_slider.setPreferredSize(new Dimension(frame.getWidth(),20));
        delay_slider.addAdjustmentListener(this);
        
        Label slow = new Label("Slow");
        Label fast = new Label("Fast");
        
        scrollPanel.add(slow);
        scrollPanel.add(delay_slider);
        scrollPanel.add(fast);
                
        JPanel northPanel = new JPanel(new GridLayout(2,0));
        northPanel.add(buttonPanel);
        northPanel.add(scrollPanel);
        
        frame.add(northPanel,BorderLayout.NORTH);
    }
    
    public void addTimeSlider(Frame frame) 
    {
        JPanel scrollPanel_time = new JPanel();

        time_slider.setPreferredSize(new Dimension(frame.getWidth(),20));
        time_slider.addAdjustmentListener(this);

        Label start_time = new Label("Start");
        Label end_time = new Label("End");

        scrollPanel_time.add(start_time);
        scrollPanel_time.add(time_slider);
        scrollPanel_time.add(end_time);

        JPanel checkPanel = new JPanel(new GridLayout(0,5));
        output.addItemListener(this);
        cells.addItemListener(this);
        axes.addItemListener(this);
        axes_equal.addItemListener(this);
        ancestors.addItemListener(this);
        potts.addItemListener(this);
        
        checkPanel.add(output);
        checkPanel.add(cells);
        checkPanel.add(axes);
        checkPanel.add(axes_equal);
        checkPanel.add(ancestors);
        checkPanel.add(potts);
        checkPanel.add(nearest_node_label);
        checkPanel.add(nearest_element_centroid_label);

        JPanel southPanel = new JPanel(new GridLayout(2,0));
        southPanel.add(scrollPanel_time);
        southPanel.add(checkPanel);
        frame.add(southPanel,BorderLayout.SOUTH);
    }    
    
    public static void main(String args[]) 
    {
        System.out.println("Copyright The Chaste Project");

        // Set default states for options
        cells.setState(true);
        output.setState(false);
        axes.setState(true);
        axes_equal.setState(false);

        // Update states for options according to input args
        for (int i=1; i<args.length; i++)
        {
            if (args[i].equals("output"))
            {
                writeFiles = true;
                output.setState(true);                  
            }
            else if (args[i].equals("nocells"))
            {
                drawCells = false;
                cells.setState(false);
            }
            else if (args[i].equals("axes"))
            {
                drawAxes = true;
                axes.setState(true);
            }
            else if (args[i].equals("axesequal"))
            {
                axesEqual = true;
                axes_equal.setState(true);
            }
            else if (args[i].equals("ancestors"))
            {
                drawAncestors = true;
                ancestors.setState(true);
            }
            else if (args[i].equals("potts"))
            {
                drawPotts = true;
                potts.setState(true);
            }
            else
            {
                System.out.println("Input option not recognised");
            }
        }

        // Read in results files for visualization
        node_file = new File(args[0] + "/results.viznodes");
        if (!node_file.isFile())
        {
            System.out.println("The file " + args[0] + "/results.viznodes doesn't exist");
            return;
        }
        
        boundary_node_file = new File(args[0] + "/results.vizboundarynodes");
        if (!boundary_node_file.isFile())
        {
            System.out.println("The file " + args[0] + "/results.vizboundarynodes doesn't exist");
            return;
        }
        
        cell_type_file = new File(args[0] + "/results.vizcelltypes");
        if (!cell_type_file.isFile())
        {
            System.out.println("The file " + args[0] + "/results.vizcelltypes doesn't exist");
            return;
        }
        
        element_file = new File(args[0] + "/results.vizelements");
        if (!element_file.isFile())
        {
            System.out.println("The file " + args[0] + "/results.vizelements doesn't exist");
            return;
        }        
        
        ancestors_file = new File(args[0]+"/results.vizancestors");
        if (!ancestors_file.isFile())
        {
            ancestors.setVisible(false);
            ancestors.setState(false);
            drawAncestors = false;
        }
        else
        {
            ancestorsFilePresent = true;
            ancestors.setState(true);
            drawAncestors = true;
        }
        
        setup_file = new File(args[0] + "/results.vizsetup");
        if (!setup_file.isFile())
        {
            System.out.println("The file " + args[0] + "/results.vizsetup doesn't exist");
        }
        else 
        {
            try
            {
                BufferedReader in_setup_file = new BufferedReader(new FileReader(setup_file));
                String line_setup = in_setup_file.readLine();  
                
                // Read setup information
                while (line_setup != null)
                {
                    StringTokenizer st_setup = new StringTokenizer(line_setup);
                    String parameter = st_setup.nextToken();
                    if (parameter.equals("MeshWidth"))  
                    {
                        crypt_width = Double.valueOf(st_setup.nextToken()).doubleValue();
                        half_width = crypt_width/2.0;
                        System.out.println("Mesh Width = " + crypt_width);
                        drawCylinder = true && drawCylinderOverride;    // this is made true only if mesh width exists
                    }
                    if (parameter.equals("Complete")) 
                    {
                	    showLastStep = true;
                    }
                    if (parameter.equals("PottsSimulation")) 
                    {
                	    drawPotts = true;
                	    drawCells = false;
                        cells.setState(false);
                    }
                    line_setup = in_setup_file.readLine();
                }
            }
            catch (Exception e) 
        	{
            	System.out.println("Error occured. Exception message: "+e.getMessage());
        	}
        }             
        
        final Visualize2dVertexCells vis = new Visualize2dVertexCells();
        
        LoadAllFiles();
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    
    public static void LoadAllFiles()
    {
    	parsed_all_files = false;
        try 
        {                  	
        	BufferedReader skim_node_file = new BufferedReader(new FileReader(node_file));

            int num_lines = 0;
            while (skim_node_file.readLine() != null) 
            {
                num_lines++;
            }
            
            // By default we don't read or print the final line, so
            // the visualiser can be run as a simulation is being run,
            // and the visualiser will work on incomplete data.
            if (num_lines>1 && !showLastStep)
            {
            	num_lines -= 1;
            }

            numTimeSteps = num_lines;
            time_slider.setMaximum(numTimeSteps);
            times = new double[num_lines];
            node_positions = new RealPoint[num_lines][];
            element_centroid_positions = new RealPoint[num_lines][];
            cell_types = new int [num_lines][];
            
            numNodes = new int[num_lines];
            image_nodes = new int[num_lines][];
            is_boundary_nodes = new int[num_lines][];
            
            numCells = new int[num_lines];
            numElements = new int[num_lines];
            element_nodes = new int[num_lines][];

            element_num_nodes = new int[num_lines][];

            if (ancestorsFilePresent)
            {
            	ancestor_values = new int[num_lines][];
            }
            
            // Create line readers
            BufferedReader in_node_file = new BufferedReader(new FileReader(node_file));
            BufferedReader in_boundary_node_file = new BufferedReader(new FileReader(boundary_node_file));
            BufferedReader in_cell_type_file = new BufferedReader(new FileReader(cell_type_file));
            BufferedReader in_element_file = new BufferedReader(new FileReader(element_file));
            
            
            // Create strings for lines to be read
            String line_node = in_node_file.readLine();
            String line_boundary_node = in_boundary_node_file.readLine();
            String line_cell_type = in_cell_type_file.readLine();
            String line_element = in_element_file.readLine();
            
            // Create line readers and strings for optional files
            String line_ancestors = "";
            BufferedReader in_ancestors_file = null;
            if (drawAncestors)
            {
               	ancestor_values = new int[num_lines][]; 
               	in_ancestors_file = new BufferedReader(new FileReader(ancestors_file));
               	line_ancestors = in_ancestors_file.readLine();
            }
            
            // If line is not end of file continue
            int time_step = 0;
            while (line_node != null && time_step<num_lines) 
            {
            	// Create string tokenizers with a colon sign as a delimiter
                StringTokenizer st_node = new StringTokenizer(line_node);
                StringTokenizer st_boundary_node = new StringTokenizer(line_boundary_node);
                StringTokenizer st_cell_type = new StringTokenizer(line_cell_type);
                StringTokenizer st_element = new StringTokenizer(line_element);
                StringTokenizer st_ancestors = null;
                
                if (drawAncestors)
                {
                    st_ancestors = new StringTokenizer(line_ancestors);
                    Double current_ancestors_time = Double.valueOf(st_ancestors.nextToken());
                }
///////////////// READ IN NODE DATA... /////////////////
                
                // Get current time from node file
                double current_time_for_nodes = Double.valueOf(st_node.nextToken()).doubleValue();
                // Get current time from boundary node file
                double current_time_for_boundary_nodes = Double.valueOf(st_boundary_node.nextToken()).doubleValue();
              
                // Check times agree
                if (Math.abs(current_time_for_nodes - current_time_for_boundary_nodes) > 1e-6)
                {
                	throw new Exception("Error: The time corresponding to each line of the boundary nodes file must match that of the node file");
                }
                                
                // Add to vector of times
                times[time_step] = current_time_for_nodes;
                
                // Count the number of entries in the node file and check correct 
                int current_num_node_entries = st_node.countTokens();       
                if (current_num_node_entries%2 != 0)
                {
                	System.out.println("Warning: The node output at time " 
                			           + current_time_for_nodes + 
                			           " is not of the required form: time,x,y,x,y,...");
                	break;
                }
                int current_num_nodes = current_num_node_entries/2;

                // Allocate memory for the current time step
                node_positions[time_step] = new RealPoint[memory_factor*current_num_nodes];
                is_boundary_nodes[time_step] = new int[memory_factor*current_num_nodes];

                for (int i=0; i<current_num_nodes; i++) 
                {
                    double d1 = Double.valueOf(st_node.nextToken()).doubleValue();
                    double d2 = Double.valueOf(st_node.nextToken()).doubleValue();
                    node_positions[time_step][i] = new RealPoint(d1, d2);
                    
                    is_boundary_nodes[time_step][i] = Integer.parseInt(st_boundary_node.nextToken());
                }
                numNodes[time_step] = current_num_nodes;

///////////////// ...NODE DATA READ /////////////////
                
///////////////// READ IN ELEMENT DATA... /////////////////
                
                // Get current time from element file
                double current_time_for_elements = Double.valueOf(st_element.nextToken()).doubleValue();
                
                // Check times agree
                if (Math.abs(current_time_for_nodes - current_time_for_elements) > 1e-6)
                {
                	throw new Exception("Error: The time corresponding to each line of the element file must match that of the node file");
                }

                int num_element_entries = st_element.countTokens();
                int current_num_elements = 0;
                int entry_position = 0;
                int num_entries_covered = 0;
                
                // Allocate memory for the current time step
                element_nodes[time_step] = new int[memory_factor*num_element_entries];
                element_num_nodes[time_step] = new int[memory_factor*num_element_entries];
                element_centroid_positions[time_step] = new RealPoint[memory_factor*num_element_entries];

                while (entry_position < num_element_entries)
                {
                 	entry_position++;
                 	current_num_elements++;

                   	int num_nodes_in_this_element = Integer.parseInt(st_element.nextToken());
                   	element_num_nodes[time_step][current_num_elements-1] = num_nodes_in_this_element;

                   	double this_element_centroid_position_x = 0.0;
                   	double this_element_centroid_position_y = 0.0;
                   	int index[] = new int[num_nodes_in_this_element];

                   	for (int i=0; i<num_nodes_in_this_element; i++)
                   	{
                   		int node_index = Integer.parseInt(st_element.nextToken());
                   		element_nodes[time_step][num_entries_covered + i] = node_index;

                   		this_element_centroid_position_x += node_positions[time_step][node_index].x/Double.valueOf(num_nodes_in_this_element).doubleValue();;
                   		this_element_centroid_position_y += node_positions[time_step][node_index].y/Double.valueOf(num_nodes_in_this_element).doubleValue();;

                   		entry_position++;
                   	}

        	       	element_centroid_positions[time_step][current_num_elements-1] = new RealPoint(this_element_centroid_position_x, this_element_centroid_position_y);

                   	num_entries_covered = num_entries_covered + num_nodes_in_this_element;
                }
                numElements[time_step] = current_num_elements;
                
///////////////// ...ELEMENT DATA READ /////////////////                
                
///////////////// READ IN CELL TYPES AND ANCESTOR DATA... /////////////////
                // Get current time from cell types file
                double current_time_for_cell_types = Double.valueOf(st_cell_type.nextToken());
                
                // Check times agree
                if (Math.abs(current_time_for_nodes - current_time_for_cell_types) > 1e-6)
                {
                	throw new Exception("Error: The time corresponding to each line of the cell types file must match that of the node file");
                }
                
                int num_cell_type_entries = st_cell_type.countTokens();
                int current_num_cell_types = 0;
                entry_position = 0;
                num_entries_covered = 0;
                
                // Allocate memory for the current time step
                cell_types[time_step] = new int[memory_factor*num_cell_type_entries];
                if (ancestorsFilePresent)
                {
                	ancestor_values[time_step] = new int[memory_factor*num_cell_type_entries];
                }
                while (entry_position < num_cell_type_entries)
                {
                 	entry_position++;
                 	current_num_cell_types++;

                 	cell_types[time_step][current_num_cell_types-1] = Integer.parseInt(st_cell_type.nextToken());
                 	if (drawAncestors)
                    {	
                        ancestor_values[time_step][current_num_cell_types-1] = Integer.parseInt(st_ancestors.nextToken());
                    }	
                }

                numCells[time_step] = current_num_cell_types;
///////////////// ...CELL TYPES DATA READ /////////////////    
                
                
///////////////// CHECK CELL AND ELEMENT NUMBERS AGREE... /////////////////
                
                // Check times agree                    
                if (Math.abs(current_time_for_nodes - current_time_for_cell_types) > 1e-6)
                {
                	throw new Exception("Error: The time corresponding to each line of the cell type file must match that of the node file");
                }
                
                if (Math.abs(current_time_for_cell_types - current_time_for_elements) > 1e-6)
                {
                	throw new Exception("Error: The time corresponding to each line of the cell type file must match that of the element file");
                }
                if (current_num_cell_types != current_num_elements)
    			{
    				System.out.println("Warning: At time " + current_time_for_cell_types + 
    						           ", element file gives " + current_num_elements + 
    						           " cells, but cell type file gives " + current_num_cell_types + 
    						           " cells");
    				break;
    			}
///////////////// ...CELL AND ELEMENT NUMBERS DO AGREE /////////////////
                
                // Read next line of each file
                line_node = in_node_file.readLine();
                line_boundary_node = in_boundary_node_file.readLine();
                line_cell_type = in_cell_type_file.readLine();
                line_element = in_element_file.readLine();
                
                if (drawAncestors)
                {
                	line_ancestors = in_ancestors_file.readLine();
                }   

                // Increment time
                time_step++;

            } // end while not at end of file
            
            System.out.println("Writing output files = " + writeFiles);
            System.out.println("Drawing cells = " + drawCells);
            System.out.println("Drawing cylindrically = " + drawCylinder);
            System.out.println("Drawing axes = " + drawAxes);
            System.out.println("Drawing axes equal = "+ axesEqual);
            System.out.println("Drawing clonal populations = "+ drawAncestors);
            
            if (drawCylinder) 
            {
            	ConvertCylindricalDataToPlane();
            }
            
            
            CalculateCanvasDimensions();
            parsed_all_files = true;
        }
        catch (Exception e) 
        {
            System.out.println("Error occured. Exception message: "+e.getMessage());
        }
    }
    
    
    public static void CalculateCanvasDimensions()
    { 	
        max_x = -1e12;
        max_y = -1e12;
        min_x =  1e12;
        min_y =  1e12;
        for (int row=0; row<numTimeSteps; row++)
        {
           for (int j=0; j<numNodes[row]; j++) 
           {
               if (node_positions[row][j].x > max_x) 
               {
                   max_x = node_positions[row][j].x;
               }
               if (node_positions[row][j].y > max_y) 
               {
                   max_y = node_positions[row][j].y;
               }
               if (node_positions[row][j].x < min_x) 
               {
                  min_x = node_positions[row][j].x;
               }
               if (node_positions[row][j].y < min_y) 
               {
                  min_y = node_positions[row][j].y;
               }
           }
        }
        lattice_sites_up = (int)max_y + 1;
        lattice_sites_across = (int)max_x + 1;
    }
    
    
    public static void ConvertCylindricalDataToPlane()
    {
        // Scan through each element
        for (int time_index=0; time_index<numTimeSteps; time_index++)
    	{
            image_nodes[time_index] = new int[memory_factor*numNodes[time_index]]; // reserve plenty of memory
            
            // Fill image_nodes with an identity map (at each time step each node maps to itself)            
            for (int i=0; i<numNodes[time_index]; i++) 
            {
                image_nodes[time_index][i] = i;
            }
            
            if (elementFilePresent)
            {
            	int num_entries_covered = 0;
            	for (int i=0; i<numElements[time_index]; i++)
                {   
                	int num_nodes_in_this_element = element_num_nodes[time_index][i];	        	
        	       	int index[] = new int[num_nodes_in_this_element];
        	       	
        	       	RealPoint real_points[] = new RealPoint[num_nodes_in_this_element];
        	       	for (int j=0; j<num_nodes_in_this_element; j++)
        	        {
        	       		// Global index of each node
        	        	index[j] = element_nodes[time_index][num_entries_covered + j];
                        // Find the co-ords of each node
        	        	real_points[j] = node_positions[time_index][index[j]];      	
        	        }
        	       	
                    // Identify pairs of nodes that are a long way apart
        	       	for (int j=0; j<num_nodes_in_this_element; j++)
        	       	{
        	       		for (int k=j+1; k<num_nodes_in_this_element; k++)
            	       	{
//        	       			System.out.println("Element " + i + " nodes = "+ j  + " and = "+ k);
        	       			if ((Math.abs(real_points[j].x - real_points[k].x) > 0.5*crypt_width))
	        	       		{
        	       				// Make the left node an image node
	        	       			int local_node_index_to_create_image;
	        	       			if (real_points[j].x < real_points[k].x)
	        	       			{
	        	       				local_node_index_to_create_image = j;
	        	       			}
	        	       			else // (real_points[j].x > real_points[k].x)
	        	       			{
	        	       				local_node_index_to_create_image = k; 
	        	       			}
	        	       			int global_node_index_to_create_image = index[local_node_index_to_create_image];
	        	       			
//	        	       			System.out.println("Element " + i + " NeedsImage = "+ local_node_index_to_create_image + " GlobalIndex = "+ global_node_index_to_create_image);
	        	       			MakeNewImageNode(time_index,global_node_index_to_create_image);
	        	       			
               	       		    //Get elements to use newly created image nodes
	        	       			element_nodes[time_index][num_entries_covered + local_node_index_to_create_image] = image_nodes[time_index][index[local_node_index_to_create_image]];
	        	       			
	        	       		}
            	       	}
        	       	}
                    num_entries_covered = num_entries_covered + num_nodes_in_this_element;
                }
            }            
        }
    }
    
    public static void MakeNewImageNode(int time_index, int node_index)
    {   
    	// Only make a new node if one hasn't already been made
        if (image_nodes[time_index][node_index] == node_index)
        {	
        	// Make a new image of Node A       
            RealPoint new_point = node_positions[time_index][node_index];
            RealPoint new_point2 = new RealPoint(0.0, 0.0);
            new_point2.y = new_point.y;
            
            if (new_point.x < half_width)
            {
                new_point2.x = new_point.x + crypt_width;
            }
            if (new_point.x > half_width)
            {
                new_point2.x = new_point.x - crypt_width;
            }

            // New image node
            node_positions[time_index][numNodes[time_index]] = new_point2;
            is_boundary_nodes[time_index][numNodes[time_index]] = is_boundary_nodes[time_index][node_index]; 
            //cell_types[time_index][numCells[time_index]] = canvas.INVISIBLE_COLOUR;
            
            // Update the image record
            image_nodes[time_index][node_index] = numNodes[time_index]; 
            numNodes[time_index]++;
        }    
    }
    
}


class RealPoint
{
    public double x, y; 
    public RealPoint(double xs, double ys)
    {
        x = xs;
        y = ys;
    }
    public RealPoint(RealPoint p1, RealPoint p2) // returns the average of p1 and p2
    {
        x = (p1.x + p2.x)/2.0;
        y = (p1.y + p2.y)/2.0;
    }
}


class PlotPoint
{
    public int x, y;
    public PlotPoint(int xs, int ys)
    {        
        x = xs;
        y = ys;
    }
}


class CustomVertexCanvas2D extends Canvas implements MouseMotionListener 
{
    private static final long serialVersionUID = 6997195399856046957L;

    Visualize2dVertexCells vis;

    boolean imageReady = false;   
    boolean imageDrawing = false;
    
    int width;
    int height;
    int node_radius = 2;
    
    BufferedImage buffered_image = null;
    Graphics2D g2 = null;
    
    
    Color background_white = new Color(255,255,255);
    Color apoptotic_grey = new Color(80,80,80);
    Color purple = new Color(121,126,234);
    
    public CustomVertexCanvas2D(Visualize2dVertexCells v) 
    {
        vis = v;
        setBackground(background_white);
    }

    public void paint(Graphics graphics)
    {
        if (vis.parsed_all_files == false)
        {
            graphics.drawString("Still parsing input...", 10,10);
            return;
        }
                
        if (!imageReady)
        {
            repaint();
        }
        imageDrawing = true;
        graphics.drawImage(buffered_image,0,0,this);
        
        if (vis.writeFiles)
        {
        	//String filename = String.format("image%1$05d.png", vis.timeStep); //Pre Java-1.6
        	String filename = String.format("image%1$05d.png", new Object[] {new Integer(vis.timeStep)});
            System.out.println("Writing file : "+filename+".");
            File f = new File(filename);
            try 
            {
                ImageIO.write(buffered_image, "png", f);
            } 
            catch (Exception e)
            {
            }
            System.out.println("Written file : "+filename+".");
        }
        imageDrawing = false;
    }
    
    public void drawBufferedImage() 
    {
        int cycle = 0;
        while (imageDrawing)
        {
            System.out.print("");
            if (cycle == 100000)
            {
                System.out.print(".");
                cycle = 0;
            }
            cycle++;
        }
        imageReady = false;
        
        int tick_length = 10;
        
        vis.time_slider.setValue(vis.timeStep); 
        
        if (g2==null)
        {
            height = getHeight();
            width = getWidth();
            buffered_image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2 = (Graphics2D) buffered_image.getGraphics();
        }
        
        g2.setColor(background_white);
        g2.fillRect(0, 0, width, height);
        g2.setColor(Color.black);
        g2.drawString("Time = " + vis.times[vis.timeStep], 10, 10);
        
        g2.setColor(Color.black);
        Shape original_clip = g2.getClip();

        // Store the nodes containting element. Only used for potts sims where this is unique.
        int containing_elements[];
        containing_elements = new int[vis.numNodes[vis.timeStep]];
        
        if (vis.drawPotts)
        {
	        // For a Potts sim Draw All nodes before anything else 
        	// to get nodes that are not in any elements, i.e. medium
			for (int i=0; i<vis.numNodes[vis.timeStep]; i++)
			{
			    PlotPoint p = scale(vis.node_positions[vis.timeStep][i]);
			
			    g2.setColor(Color.black);
			    
			    g2.fillOval(p.x - node_radius, p.y - node_radius, 2 * node_radius, 2 * node_radius);
			    
			    containing_elements[i] = -1; // Default of medium
			}
			
	        // Calculate the containing elements for each node
			int num_entries_covered = 0;
			for (int i=0; i<vis.numElements[vis.timeStep]; i++)
		    {	
		      	int num_nodes_in_this_element = vis.element_num_nodes[vis.timeStep][i];	        	
		       	int global_index;
		       	for (int j=0; j<num_nodes_in_this_element; j++)
		        {
		        	global_index = vis.element_nodes[vis.timeStep][num_entries_covered + j];
		        	containing_elements[global_index] = i;
		        }
		       	num_entries_covered = num_entries_covered + num_nodes_in_this_element;
			}
			
        }
                
	    // Draw elements first
        int num_entries_covered = 0;
	    for (int i=0; i < vis.numElements[vis.timeStep]; i++)
	    {	
	      	int num_nodes_in_this_element = vis.element_num_nodes[vis.timeStep][i];	        	
	       	int index[] = new int[num_nodes_in_this_element];
	       	
	       	PlotPoint vertices[] = new PlotPoint[num_nodes_in_this_element];
	       	RealPoint real_points[] = new RealPoint[num_nodes_in_this_element];
	       	for (int j=0; j<num_nodes_in_this_element; j++)
	        {
	        	// What nodes are we joining up?
	        	index[j] = vis.element_nodes[vis.timeStep][num_entries_covered + j];
				RealPoint this_node_location = vis.node_positions[vis.timeStep][index[j]];
	        	real_points[j] = vis.node_positions[vis.timeStep][index[j]];
	        }
	       	
	       	num_entries_covered = num_entries_covered + num_nodes_in_this_element;
	       	
	        // Where are they? Convert to integer pixels	           
	        for (int j=0; j<num_nodes_in_this_element; j++)
	        {
	        	vertices[j] = scale(real_points[j]);
	        }
	        
//	        if (vis.cell_types[vis.timeStep][i]== LATE_CANCER_COLOUR)
//	        {
		        if (vis.drawCells && !vis.drawPotts)
	            {
	                int xpoints[] = new int[num_nodes_in_this_element];
	                int ypoints[] = new int[num_nodes_in_this_element];
	                for (int node=0; node<num_nodes_in_this_element; node++)
	                {
	                	xpoints[node] = vertices[node].x;
	                	ypoints[node] = vertices[node].y;
	                }
	                SetCellColour(i);
	    	        
	                g2.fillPolygon(xpoints, ypoints, num_nodes_in_this_element);
	
	                // Plot cell boundary lines
	                g2.setColor(Color.black);
	                
	                for (int node=0; node<num_nodes_in_this_element; node++)
	                {
	                	g2.drawLine(xpoints[node], 
	                			    ypoints[node], 
	                			    xpoints[(node+1)%num_nodes_in_this_element],
	                			    ypoints[(node+1)%num_nodes_in_this_element]);
	                }
	            }	
//	        }
	    }    

        // Draw nodes second so that dots are on top of lines Loop over elements so dont plot image nodes
	    num_entries_covered = 0;
	    for (int i=0; i < vis.numElements[vis.timeStep]; i++)
	    {	
	      	int num_nodes_in_this_element = vis.element_num_nodes[vis.timeStep][i];	        	
	       	int global_index;
	       	PlotPoint point;
	       	RealPoint real_point;
	       	for (int j=0; j<num_nodes_in_this_element; j++)
	        {
	        	global_index = vis.element_nodes[vis.timeStep][num_entries_covered + j];
				real_point = vis.node_positions[vis.timeStep][global_index];      	
	        
	        	// Where are they? Convert to integer pixels	           
		        point = scale(real_point);

	        	SetNodeColour(global_index);

		        // \todo: Larger simulations would be clearer with smaller nodes
	        	g2.fillOval(point.x - node_radius, point.y - node_radius, 2 * node_radius, 2 * node_radius);
	        	
	        	if (vis.drawPotts)
	        	{
	        		int nodes_across = vis.lattice_sites_across;
		            int nodes_up = vis.lattice_sites_up;
	        		
	        		PlotPoint square_vertices[] = new PlotPoint[4];
	    	       	
	    	       	square_vertices[0] = scale(real_point.x - 0.5, real_point.y - 0.5);
    	        	square_vertices[1] = scale(real_point.x + 0.5, real_point.y - 0.5);
    	        	square_vertices[2] = scale(real_point.x + 0.5, real_point.y + 0.5);
    	        	square_vertices[3] = scale(real_point.x - 0.5, real_point.y + 0.5);
	    	        
	    	        int xpoints[] = new int[4];
	                int ypoints[] = new int[4];
	                for (int node=0; node<4; node++)
	                {
	                	xpoints[node] = square_vertices[node].x;
	                	ypoints[node] = square_vertices[node].y;
	                }
	                SetCellColour(i);
	    	        
	                g2.fillPolygon(xpoints, ypoints, 4);
	
	                // Plot cell boundary lines with dashed lines
	                g2.setColor(Color.black);
 
	                Stroke dashed_stroke = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[] { 1, 0 }, 0);
	                g2.setStroke(dashed_stroke);
	                
	                for (int node=0; node<4; node++)
	                {
	                	g2.drawLine(xpoints[node], 
	                			    ypoints[node], 
	                			    xpoints[(node+1)%4],
	                			    ypoints[(node+1)%4]);
	                }
	                
	                // Find cell edges and plot them with solid lines
	                Stroke solid_stroke = new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[] { 1 , 0}, 0);
	                g2.setStroke(solid_stroke);

	                if (global_index >= nodes_across) // not on bottom row 
	                {
	                	int bottom_neighbour = global_index - nodes_across;
	                	if (containing_elements[bottom_neighbour] != i) //i is current element
	                	{
	                		g2.drawLine(xpoints[0], 
		                			    ypoints[0], 
		                			    xpoints[1],
		                			    ypoints[1]);
		                }
	                }
	                else
	                {
	                	g2.drawLine(xpoints[0], 
                			        ypoints[0], 
                			        xpoints[1],
                			    	ypoints[1]);
	                }
	        		if (global_index < nodes_across*(nodes_up-1) ) // not on top row 
	        		{
	        			int top_neighbour = global_index + nodes_across;
	                	if (containing_elements[top_neighbour] != i) //i is current element
	                	{
	                		g2.drawLine(xpoints[2], 
		                			    ypoints[2], 
		                			    xpoints[3],
		                			    ypoints[3]);
		                }
	        		}
	        		else
	        		{
	                	g2.drawLine(xpoints[2], 
                			        ypoints[2], 
                			        xpoints[3],
                			        ypoints[3]);
	                }
	        		if (global_index%nodes_across != 0 ) // Not on left hand side
	        		{
	        			int left_neighbour = global_index - 1;
	                	if (containing_elements[left_neighbour] != i) //i is current element
	                	{
	                		g2.drawLine(xpoints[0], 
		                			    ypoints[0], 
		                			    xpoints[3],
		                			    ypoints[3]);
		                }
	        		}
	        		else
	        		{
	                	g2.drawLine(xpoints[0], 
                			        ypoints[0], 
                			        xpoints[3],
                			        ypoints[3]);
	                }
	        		if (global_index%nodes_across != (nodes_across - 1)) // Not on right hand side
	                {
	                	int right_neighbour = global_index + 1;
	                	if (containing_elements[right_neighbour] != i) //i is current element
		                {
	                		g2.drawLine(xpoints[1], 
		                			    ypoints[1], 
		                			    xpoints[2],
		                			    ypoints[2]);
		                }
	                }
	        		else
	        		{
	                	g2.drawLine(xpoints[1], 
                			        ypoints[1], 
                			        xpoints[2],
                			        ypoints[2]);
	                }
	        	}
	        }
	        num_entries_covered = num_entries_covered + num_nodes_in_this_element;
        }
	    
	    g2.setColor(Color.black);
        
        if (vis.drawAxes)
        {
        	drawXAxis(tick_length);
        	drawYAxis(tick_length);
        }
        
        imageReady = true;
    }
    
    private void drawXAxis(int tick_length) 
    {
        int min_x = (int) vis.min_x;
        if (vis.min_x<0)
        {
        	min_x -= 1;
        }
//Hack To get matching axes        
//min_x=-1;
        int max_x = (int) vis.max_x;
        if (vis.max_x>0)
        {
        	max_x += 1;
        }
//Hack To get matching axes
//max_x=8;
        // work out the number of ticks to use - if it's too big (>15, say)
        // this would make the axis look crowded, so halve num_ticks  
        if ((max_x-min_x)%2 != 0)
        {
        	max_x++;
        }
        int num_ticks = max_x - min_x;        
        int tick_spacing = 1;
        if (num_ticks > 11) 
        {        	
        	num_ticks = num_ticks/2;
        	tick_spacing = 2;
        }
        
        PlotPoint start = scale(min_x, 0);
        PlotPoint end = scale(max_x, 0);
        g2.drawLine(start.x, start.y, end.x, end.y);
          
                
        for (int i = 0; i <= num_ticks; i++) 
        {
            double x = (double) (min_x + tick_spacing*i);
            DecimalFormat df = new DecimalFormat("0.0");
            String x_1dp = df.format(x);
              
            // Tick lines
            PlotPoint posn =  scale(x,0);
            g2.drawLine(posn.x, posn.y, posn.x, posn.y+tick_length);
            g2.drawString(x_1dp, posn.x, posn.y + 2*tick_length);
        }
    }
        
    private void drawYAxis(int tick_length) 
    {
        int min_y = (int) vis.min_y;
        if (vis.min_y<0)
        {
        	min_y -= 1;
        }
// Hack To get matching axes        
//min_y=-1;
        int max_y = (int) vis.max_y;
        if (vis.max_y>0)
        {
        	max_y += 1;
        }
// Hack To get matching axes
//max_y=9;     
        
        // work out the number of ticks to use - if it's too big (>15, say)
        // this would make the axis look crowded, so halve num_ticks  
        if ((max_y-min_y)%2 != 0)
        {
        	max_y++;
        }
        int num_ticks = max_y - min_y;        
        int tick_spacing = 1;
        if (num_ticks > 11) 
        {        	
        	num_ticks = num_ticks/2;
        	tick_spacing = 2;
        }
                
        PlotPoint start = scale(0, min_y);
        PlotPoint end = scale(0, max_y);        
        g2.drawLine(start.x, start.y, end.x, end.y);
        
        for (int i = 0; i <= num_ticks; i++) 
        {
            double y = (double) (min_y + tick_spacing*i);
            DecimalFormat df = new DecimalFormat("0.0");
            String y_1dp = df.format(y);

            // Tick lines
            PlotPoint posn = scale(0,y);
            g2.drawLine(posn.x-tick_length, posn.y, posn.x, posn.y);
            g2.drawString(y_1dp, posn.x - 4*tick_length, posn.y );
        }
    }
    
    PlotPoint scale(double x, double y)
    {
        // Map min_x to eps and max_x to width-eps (to allow a border)
        int eps = 100;
        int ix = (int) ((x - vis.min_x) * (width-2*eps) /(vis.max_x - vis.min_x) +eps);
        int iy = (int) ((y - vis.min_y) * (height-2*eps) /(vis.max_y - vis.min_y) +eps);
        iy = height - iy; // this is because java is silly and has the y axis going down the screen
        return (new PlotPoint(ix,iy));
    }
    
    PlotPoint scale(RealPoint p) 
    {
        return (scale(p.x,p.y));    
    }
    
    RealPoint unscale(PlotPoint p)
    {
        int ix = p.x;
        int iy = height - p.y;
        int eps = 20;
        
        double x = (ix-eps)*(vis.max_x - vis.min_x) / (width-2.0*eps) + vis.min_x;
        double y = (iy-eps)*(vis.max_y - vis.min_y) / (height-2.0*eps) + vis.min_y;
        
        return (new RealPoint(x,y));
    }

    public static double RoundedDouble(double value, int numDecimalPlaces)
    {
    	double temp = (double)Math.pow(10, numDecimalPlaces);
    	value = value * temp;
    	double temp2 = Math.round(value);
    	return (double)temp2/temp;
    }
    
    public void mouseMoved(MouseEvent e) 
    {
        PlotPoint mouse_position = new PlotPoint(e.getX(), e.getY());
        RealPoint real_position = unscale(mouse_position);

        int nearest_node_index = -1;
        for (int i=0; i<vis.numNodes[vis.timeStep]; i++) 
        {
            int sq_dist = SquaredDistance(scale(vis.node_positions[vis.timeStep][i]), mouse_position);
            if (sq_dist < node_radius*node_radius)
            {
            	nearest_node_index = i;
                break;
            }
        }
        if (nearest_node_index >= 0)
        {
            RealPoint node_position = vis.node_positions[vis.timeStep][nearest_node_index];
            vis.nearest_node_label.setText("Node " + nearest_node_index + " is at " 
            		                               + RoundedDouble(node_position.x, 2) + "  " 
            		                               + RoundedDouble(node_position.y, 2));
        }
        else
        {
            vis.nearest_node_label.setText("");
        }

        int nearest_element_centroid_index = -1;
        for (int i=0; i<vis.numElements[vis.timeStep]; i++) 
        {
            int sq_dist = SquaredDistance(scale(vis.element_centroid_positions[vis.timeStep][i]), mouse_position);
            if (sq_dist < node_radius*node_radius)
            {
            	nearest_element_centroid_index = i;
                break;
            }
        }
        if (nearest_element_centroid_index >= 0)
        {
            RealPoint element_centroid_position = vis.element_centroid_positions[vis.timeStep][nearest_element_centroid_index];
            vis.nearest_element_centroid_label.setText("Element " + nearest_element_centroid_index);
        }
        else
        {
            vis.nearest_element_centroid_label.setText("");
        }
    }
    
    public void mouseDragged(MouseEvent e) 
    {
        // Not used
    }
    
    int SquaredDistance(PlotPoint p0, PlotPoint p1)
    {
        int diffx = p0.x-p1.x;
        int diffy = p0.y-p1.y;
        return diffx*diffx + diffy*diffy;
    }
     
    void SetNodeColour(int index)
    {
    	if(vis.drawPotts)
        {
            Color cell_colour = ancestorColourMap(index);
            int new_r = 0;
            int new_g = 0;
            int new_b = 0;
            if (cell_colour.getRed() - 40 > new_r) new_r = cell_colour.getRed() - 40;
            if (cell_colour.getGreen() - 40 > new_g) new_g = cell_colour.getGreen() - 40;
            if (cell_colour.getBlue() - 40 > new_b) new_b = cell_colour.getBlue() - 40;
            g2.setColor(new Color(new_r, new_g, new_b));
        }
    	else
    	{
	    	if(vis.is_boundary_nodes[vis.timeStep][index]==1) 
	    	{
	    		// Boundary nodes are Red
	    		g2.setColor(Color.red);	
	    	}
	    	else
	    	{
	    		// All other nodes are Black
	        	g2.setColor(Color.black);
	    	}
    	}
    }
 
     
    /* From C++ 
    enum cell_colours
    {
        STEM_COLOUR,//0
        TRANSIT_COLOUR,//1
        DIFFERENTIATED_COLOUR,//2
        EARLY_CANCER_COLOUR,//3
        LATE_CANCER_COLOUR,//4
        LABELLED_COLOUR,//5
        APOPTOSIS_COLOUR,//6
        INVISIBLE_COLOUR, // visualizer treats '7' as invisible
    };
    */

    public static final int STEM_COLOUR = 0;
    public static final int TRANSIT_COLOUR = 1;
    public static final int DIFFERENTIATED_COLOUR = 2;
    public static final int EARLY_CANCER_COLOUR = 3;
    public static final int LATE_CANCER_COLOUR = 4;
    public static final int LABELLED_COLOUR = 5;
    public static final int APOPTOSIS_COLOUR = 6;
    public static final int INVISIBLE_COLOUR = 7;
    
    void SetCellColour(int index)
    {
        if (vis.drawAncestors && (vis.ancestor_values[vis.timeStep][index]!=-1))
      	{	// If we are drawing ancestors and this cell's value has been set in simulation.
    		Color ancestor_colour = ancestorColourMap(vis.ancestor_values[vis.timeStep][index]);
    		g2.setColor(ancestor_colour);
      	}
    	else
        {
    		switch (vis.cell_types[vis.timeStep][index]) 
            {
	        	case STEM_COLOUR: // stem cell
	        		g2.setColor(Color.cyan); 
	        		break;
	        	case TRANSIT_COLOUR: // transit cell
	        		g2.setColor(Color.yellow); 
	        		break;
	        	case DIFFERENTIATED_COLOUR: // differentiated cell
	        		g2.setColor(Color.pink); 
	        		break;
	        	case EARLY_CANCER_COLOUR: // early cancer
	        		g2.setColor(Color.lightGray); 
	        		break;
	        	case LATE_CANCER_COLOUR:  // late cancer
	        		g2.setColor(Color.gray);
	        		break;
	        	case LABELLED_COLOUR: // labelled cell
	        		g2.setColor(purple); 
	        		break;
	        	case APOPTOSIS_COLOUR: // apoptotic cell
	        		g2.setColor(apoptotic_grey); 
	        		break;
	        	case INVISIBLE_COLOUR: // sloughed cell
	        		g2.setColor(background_white);                    
	        		break;
	        	default: 
	        		g2.setColor(background_white);                    
	        		break;
	        }
        }
    }

    public Color ancestorColourMap(int ancestor)
    {
        //Map the colour uniquely into [0, 255]
        int r=hash32shiftmult(ancestor, 256);
        int g=hash32shiftmult(ancestor+1, 256);
        int b=hash32shiftmult(ancestor*2, 256);
        return new Color(r,g,b);
    }
    
    public int hash32shiftmult(int key, int range)
    {
      //Mostly copied from http://www.acme.com/resources/classes/Acme/IntHashtable.java
      int c2=0x27d4eb2d; // a prime or an odd constant
      key = (key ^ 61) ^ (key >>> 16);
      key = key + (key << 3);
      key = key ^ (key >>> 4);
      key = key * c2;
      key = key ^ (key >>> 15);
      //We added the last two lines
      key = key & 0x7FFFFFFF; //Make positive unsigned
      return (key%range);//In 0<=key<range
    }
    
 }
