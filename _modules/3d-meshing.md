---
name: 3D meshing example
category: academic, electromagnetism
layout: module
---

# Coil

Build a 3D coil mesh and solve the Poisson's problem.

## Algorithms

### 3D

3D coil mesh and Poisson's problem solve.

{% highlight cpp %}
load "medit"
load "msh3"
load "tetgen"
load "iovtk"

// COIL PARAMETERS
real R = 0.50;  // tube center to coil center distance
real r = 0.10;  // tube radius
real a = 0.05;  // stretching factor (depends on r)
int  l = 6   ;  // number of rotations

// TRANSFORMATION FROM SQUARE TO COIL
func f1 = ( R + r*cos(x) )*cos(y);
func f2 = ( R + r*cos(x) )*sin(y);
func f3 = a*y + r*sin(x);

// COIL SURFACE MESH
int   resx = 10*l, resy = 60*l; // coil surface mesh resolution
mesh  T2coil = square( resx, resy, [2.0*pi*x, 2.0*pi*y*l] );
mesh3 T3coil = movemesh23( T2coil, transfo = [f1, f2, f3] );

medit( "Open coil surface mesh", T3coil );

// CAPS TO CLOSE THE SURFACE
border C(t = 0, 2*pi){x = r*cos(t); y = r*sin(t);}
mesh  T2cap1 = buildmesh( C(resx) );
mesh3 T3cap1 = movemesh23( T2cap1, transfo = [R + x, 0, y] );
mesh3 T3cap2 = movemesh23( T2cap1, transfo = [R + x, 0, 0.5*pi+y] );

// GLUE SURFACE MESHES
mesh3 Thcoil = T3coil + T3cap1 + T3cap2;

medit( "Closed coil surface mesh", Thcoil );

// SURFACE TO VOLUME MESH
mesh3 Tcoil = tetg( Thcoil, switch = "pq1AAY" );

medit( "Coil volume mesh", Tcoil );

// SPHERE SURFACE MESH
real  rho   = 0.5*pi*l;
mesh  T2sph = square( 40*R, 80*R, [x*pi-0.5*pi,2*y*pi] );
mesh3 T3sph = movemesh23(T2sph, transfo = [rho*R*cos(x)*cos(y), rho*R*cos(x)*sin(y), (a*l)*pi + rho*R*sin(x)]);

// GLUE COIL AND SPHERE SURFACE MESHES
mesh3 Th = Thcoil + T3sph;

medit( "Combined surface mesh",Th);

// MAKE COMPLETE MESH
real[int] domain = [0, 0, 0, 1, 0.5*r, 0, R, 0, 2, 0.1*r];
mesh3 T = tetg(Th, switch = "pqaAAYYQ", nbofregions = 2, regionlist = domain);

// -----------------------------------------------------------------------------
// SAMPLE POISSON PROBLEM div( p.grad(u) ) = -z; u = 0 ON SPHERE
// -----------------------------------------------------------------------------

// FEM SPACES AND FUNCTIONS
fespace V0(T, P0); fespace V1(T, P1);
V0 mat, h = hTriangle; V1 u, v;

// DOMAIN INDICATORS
int reg1 = T(0,R,0).region; // coil is 2, see real[int] domain array
int reg2 = T(0,0,0).region; // void is 1, see real[int] domain array
cout << "Domain indicators" << endl;
cout << "Coil: " << reg1 << endl;
cout << "Void: " << reg2 << endl;

// MATERIAL FUNCTION
V0 p = 1.0*( region == reg1 ) +  1e-2*( region == reg2 );

// SET LABEL FOR OUTER SPHERE TO 100
func newlabel = ( x^2+y^2+(z-a*l*pi)^2 > rho*R-h[].min ) ? 100 : 0;
T = change(T, flabel= newlabel);

// SOLVE POISSON
solve poisson(u, v) = int3d(T)( p*(dx(u)*dx(v) + dy(u)*dy(v) + dz(u)*dz(v)) )
                    - int3d(T)( z*v )
                    + on(100, u=0);

medit("Solution", T, u);
{% endhighlight %}

|Geometry|
|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/coilGeometry.jpg)|

|Result|
|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/coilSolution.jpg)|

## Authors

Author: [Fotios Kasolis](mailto:fotios.kasolis@gmail.com)
