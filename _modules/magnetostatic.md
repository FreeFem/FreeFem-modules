---
name: Magnetostatic
category: electromagnetism
layout: module
author: Simon Garnotel
authorURL: https://github.com/sgarnotel
---

# Stokes

Algorithms for solving the 2D and 3D linear magnetostatic equations.

## Problem

Solve:

$
\displaystyle{todo}
$


## Variational form

$
\displaystyle{
  \int_{\Omega}{\frac{1}{\mu}\nabla\times\mathbf{A}\cdot\nabla\times\mathbf{AA}} - \int_{\Omega}{\mathbf{J}\cdot\mathbf{AA}} = 0
}
$

## Algorithms

### 2D

{% highlight cpp %}
//Parameters
real Mu0 = 4.*pi*1.e-7;	//Vacuum magnetic permeability
real MuC = 1.25e-6;		//Copper magnetic permeability
real J0 = 500.e4;		//Current density

//Mesh
int nn = 25;						//Mesh quality
real BoxRadius = 2.;				//"Inifite" boundary
real MagnetHeight = 1.;				//Magnet height
real MagnetInternalRadius = 0.1;	//Magnet internal radius
real MagnetExternalRadius = 0.2;	//Magnet external radius

int BoxWall = 10;

int nnL = max(2., BoxRadius*nn/2.);
border b1(t=0., 2.*pi){x=BoxRadius*cos(t); y=BoxRadius*sin(t); label=BoxWall;};

border b2(t=0., 1.){x=MagnetInternalRadius+(MagnetExternalRadius-MagnetInternalRadius)*t; y=-MagnetHeight/2.;};
border b3(t=0., 1.){x=MagnetExternalRadius; y=-MagnetHeight/2.+MagnetHeight*t;};
border b4(t=0., 1.){x=MagnetExternalRadius-(MagnetExternalRadius-MagnetInternalRadius)*t; y=MagnetHeight/2.;};
border b5(t=0., 1.){x=MagnetInternalRadius; y=MagnetHeight/2.-MagnetHeight*t;};

border b6(t=0., 1.){x=-MagnetInternalRadius-(MagnetExternalRadius-MagnetInternalRadius)*t; y=-MagnetHeight/2.;};
border b7(t=0., 1.){x=-MagnetExternalRadius; y=-MagnetHeight/2.+MagnetHeight*t;};
border b8(t=0., 1.){x=-MagnetExternalRadius+(MagnetExternalRadius-MagnetInternalRadius)*t; y=MagnetHeight/2.;};
border b9(t=0., 1.){x=-MagnetInternalRadius; y=MagnetHeight/2.-MagnetHeight*t;};

int nnH = max(2., 10.*nn*MagnetHeight);
int nnR = max(2., 10.*nn*(MagnetExternalRadius-MagnetInternalRadius));
mesh Th = buildmesh(
	  b1(nnL)
	+ b2(nnR) + b3(nnH) + b4(nnR) + b5(nnH)
	+ b6(-nnR) + b7(-nnH) + b8(-nnR) + b9(-nnH)
);

int Coil1 = Th((MagnetExternalRadius+MagnetInternalRadius)/2., 0.).region;
int Coil2 = Th(-(MagnetExternalRadius+MagnetInternalRadius)/2., 0.).region;
int Box = Th(0., 0.).region;

//Fespace
func Pk = P2;
fespace Ah(Th, Pk);
Ah Az;

//Macro
macro Curl(A) [dy(A), -dx(A)] //

//Current
Ah Jz = J0*(-1.*(Th(x,y).region == Coil1) + 1.*(Th(x,y).region == Coil2));
plot(Jz, nbiso=30, fill=true, value=true, cmm="J");

//Problem
func Mu = Mu0 + (MuC-Mu0)*((region==Coil1)+(region==Coil2));
func Nu = 1./Mu;
varf vLaplacian (Az, AAz)
	= int2d(Th)(
		  Nu * Curl(Az)' * Curl(AAz)
	)
	- int2d(Th)(
		  Jz * AAz
	)
	+ on(BoxWall, Az=0)
	;
matrix<real> Laplacian = vLaplacian(Ah, Ah, solver=sparsesolver);
real[int] LaplacianBoundary = vLaplacian(0, Ah);
Az[] = Laplacian^-1 * LaplacianBoundary;

//Magnetic induction
Ah Bx, By;
Bx = dy(Az);
By = -dx(Az);
Ah B = sqrt(Bx^2 + By^2);

//Magnetic field
Ah Hx, Hy;
Hx = Nu * Bx;
Hy = Nu * By;
Ah H = sqrt(Hx^2 + Hy^2);

//Current
Ah J = dx(Hy) - dy(Hx);

//Lorentz force
Ah Lx, Ly;
Lx = -J * By;
Ly = J * Bx;
Ah L = sqrt(Lx^2 + Ly^2);

//Plot
plot(Az, nbiso=30, fill=true, value=true, cmm="A");
plot(B, nbiso=30, fill=true, value=true, cmm="B");
plot(H, nbiso=30, fill=true, value=true, cmm="H");
plot(L, nbiso=30, fill=true, value=true, cmm="L");
{% endhighlight %}

|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/Magnetostatic2D.png)|
|Result|

### 3D

{% highlight cpp %}
load "gmsh"

//Parameters
real Mu0 = 4.*pi*1.e-7;	//Vacuum magnetic permeability
real MuC = 1.25e-6;		//Copper magnetic permeability
real J0 = 500.e4;		//Current density

//Mesh
real R = 2.;
real H = 1.;
real Rint = 0.1;
real Rext = 0.2;

int Coil = 10;
int BoxWall = 1;

mesh3 Th = gmshload3("Linear_Magnetostatic_3D.msh");

//FESpaces
func Pk = P2;
fespace Ah(Th, [Pk, Pk, Pk]);
Ah [Ax, Ay, Az] = [0, 0, 0];

//Macro
macro Curl(Ax, Ay, Az) [dy(Az)-dz(Ay), dz(Ax)-dx(Az), dx(Ay)-dy(Ax)] //
macro Divergence(Ax, Ay, Az) (dx(Ax) + dy(Ay) + dz(Az)) //

//Current
func r = sqrt(x^2+z^2);
Ah [Jx, Jy, Jz] = J0*[
	cos(atan2(z, x)+pi/2.) * (r >= Rint)*(r <= Rext)*(y >= -H/2.)*(y <= H/2.),
	0,
	sin(atan2(z, x)+pi/2.) * (r >= Rint)*(r <= Rext)*(y >= -H/2.)*(y <= H/2.)
	];

//Problem
func Mu = Mu0 + (MuC-Mu0)*(region==Coil);
func Nu = 1. / Mu;
varf vLaplacian ([Ax, Ay, Az], [AAx, AAy, AAz])
	= int3d(Th)(
		  Nu * Curl(Ax, Ay, Az)' * Curl(AAx, AAy, AAz)
		+ (1./Mu) * Divergence(Ax, Ay, Az) * Divergence(AAx, AAy, AAz)
	)
	- int3d(Th)(
		  [Jx, Jy, Jz]' * [AAx, AAy, AAz]
	)
	+ on(BoxWall, Ax=0, Ay=0, Az=0)
	;

matrix<real> Laplacian = vLaplacian(Ah, Ah, solver=sparsesolver);
real[int] LaplacianBoundary = vLaplacian(0, Ah);
Ax[] = Laplacian^-1 * LaplacianBoundary;
{% endhighlight %}

|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/Magnetostatic3D.png)|
|Result|