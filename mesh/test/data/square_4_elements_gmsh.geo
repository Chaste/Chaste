// Gmsh project created on Thu Dec 13 15:23:35 2012
Point(1) = {0, 0, 0, 1.5};
Point(2) = {1, 0, 0, 1.5};
Point(3) = {1, 1, 0, 1.5};
Point(4) = {0, 1, 0, 1.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(2) = {5};
Physical Surface(1) = {2};
Physical Line(2) = {3,4,1,2};

//Uncomment for quadratic elements
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 1;

//To regenerate the mesh file, use: gmsh -2 -o quad_square_4_elements_gmsh.msh square_4_elements_gmsh.geo