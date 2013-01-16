// Gmsh project created on Thu Dec 13 15:23:35 2012
Point(1) = {0, 0, 0, 2.5};
Point(2) = {1, 0, 0, 2.5};
Point(3) = {1, 1, 0, 2.5};
Point(4) = {0, 1, 0, 2.5};
Point(5) = {0, 0, 1, 2.5};
Point(6) = {1, 0, 1, 2.5};
Point(7) = {1, 1, 1, 2.5};
Point(8) = {0, 1, 1, 2.5};

// z=0, points 1,2,3,4
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};

// z=1, points 5,6,7,8
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};
Line Loop(10) = {6, 7, 8, 9};

// x=0, points 1,4,5,8
Line(11) = {1, 4};
Line(12) = {4, 8};
Line(13) = {8, 5};
Line(14) = {5, 1};
Line Loop(15) = {11, 12, 13, 14};

// x=1, points 2,3,6,7
Line(16) = {2,3};
Line(17) = {3,7};
Line(18) = {7,6};
Line(19) = {6,2};
Line Loop(20) = {16, 17, 18, 19};

// y=0, points 1,2,5,6
Line(21) = {1,2};
Line(22) = {2,6};
Line(23) = {6,5};
Line(24) = {5,1};
Line Loop(25) = {21, 22, 23, 24};

// y=1, points 3,4,7,8
Line(26) = {3,4};
Line(27) = {4,8};
Line(28) = {8,7};
Line(29) = {7,3};
Line Loop(30) = {26, 27, 28, 29};

Plane Surface(1) = {5};
Plane Surface(2) = {10};
Plane Surface(3) = {15};
Plane Surface(4) = {20};
Plane Surface(5) = {25};
Plane Surface(6) = {30};

Surface Loop(7) = {1,2,3,4,5,6};
Volume(8) = {7};

Physical Volume(1) = {8};
Physical Surface(2) = {1,2,3,4,5,6};

//Uncomment for quadratic elements
//Mesh.ElementOrder = 2;
//Mesh.SecondOrderLinear = 1;

//To regenerate the mesh file, use: gmsh -3 -o quad_cube_gmsh.msh simple_cube_gmsh.geo



