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
Line Loop(5) = {1, 2, 3, 4};

// z=1, points 5,6,7,8
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};
Line Loop(10) = {6, 7, 8, 9};

// x=0, points 1,4,5,8
Line(12) = {4, 8};
Line(14) = {5, 1};
Line Loop(15) = {-4, 12, 9, 14};

// x=1, points 2,3,6,7
Line(17) = {3,7};
Line(19) = {6,2};
Line Loop(20) = {2, 17, -7, 19};

// y=0, points 1,2,5,6
Line Loop(25) = {1,-19, -6, 14};

// y=1, points 3,4,7,8
Line Loop(30) = {3, 12, -8, -17};

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



