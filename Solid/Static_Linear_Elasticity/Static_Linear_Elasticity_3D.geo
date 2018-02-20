/*
Static_Linear_Elasticity_3D.geo

Author: Simon Garnotel
Date: 11/12/2017
*/

Mesh.Optimize = 1;

////////////////
///PARAMETERS///
////////////////
h = 1./5.;			//Mesh quality
L = 20.;			//Beam length
D = 1.;				//Beam height
Fixed = 1;			//Beam fixed label
Free = 2;			//Beam free label

////////////////
///ELEMENTARY///
////////////////
//Points
p = newp;
Point(p+0) = { 0.,   0.,   0.};
Point(p+1) = { D/2., 0.,   0., h};
Point(p+2) = { 0.,   D/2., 0., h};
Point(p+3) = {-D/2., 0.,   0., h};
Point(p+4) = { 0.,  -D/2., 0., h};

Point(p+5) = { 0.,   0.,   L};
Point(p+6) = { D/2., 0.,   L, h};
Point(p+7) = { 0.,   D/2., L, h};
Point(p+8) = {-D/2., 0.,   L, h};
Point(p+9) = { 0.,  -D/2., L, h};

//Lines
l = newl;
Circle(l+0) = {p+1, p+0, p+2};
Circle(l+1) = {p+2, p+0, p+3};
Circle(l+2) = {p+3, p+0, p+4};
Circle(l+3) = {p+4, p+0, p+1};

Circle(l+4) = {p+6, p+5, p+7};
Circle(l+5) = {p+7, p+5, p+8};
Circle(l+6) = {p+8, p+5, p+9};
Circle(l+7) = {p+9, p+5, p+6};

Line(l+10) = {p+1, p+6};
Line(l+11) = {p+2, p+7};
Line(l+12) = {p+3, p+8};
Line(l+13) = {p+4, p+9};

//Line Loops
ll = newll;
Line Loop(ll+0) = {l+0, l+1, l+2, l+3};
Line Loop(ll+1) = {l+4, l+5, l+6, l+7};
Line Loop(ll+2) = {l+0, l+11, -(l+4), -(l+10)};
Line Loop(ll+3) = {l+1, l+12, -(l+5), -(l+11)};
Line Loop(ll+4) = {l+2, l+13, -(l+6), -(l+12)};
Line Loop(ll+5) = {l+3, l+10, -(l+7), -(l+13)};

//Surfaces
s = news;
Plane Surface(s+0) = {ll+0};
Plane Surface(s+1) = {ll+1};
Ruled Surface(s+2) = {ll+2};
Ruled Surface(s+3) = {ll+3};
Ruled Surface(s+4) = {ll+4};
Ruled Surface(s+5) = {ll+5};

//Surface loops
sl = newsl;
Surface Loop(sl+0) = {s+0, s+1, s+2, s+3, s+4, s+5};

//Volumes
v = newv;
Volume(v+0) = {sl+0};

//////////////
///PHYSICAL///
//////////////
//Surfaces
Physical Surface("Fixed", Fixed) = {s+0, s+1};
Physical Surface("Free", Free) = {s+2, s+3, s+4, s+5};

//Volumes
Physical Volume("Volume", 1) = {v+0};



