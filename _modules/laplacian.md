---
name: Laplacian
category: academic
layout: module
---

# Laplacian

Algorithms for solving the 2D and 3D Laplace equation

## Problem

Let $u\in\mathbb{R}$, solve:

$
-\displaystyle{\Delta u = f}
$

With $f$ an external force.

## Variational form

Let $\Omega\in\mathbb{R}^n$, $2\leq n\leq3$. Let $u, v\in H^1(\Omega)$. The variational form reads as follow:

$
\displaystyle{-\int_{\Omega}{\Delta u v} - \int_{\Omega}{f v} = 0}
$

Using the Green formula:

$
\displaystyle{\int_{\Omega}{\nabla u \cdot \nabla v} - \int_{\Omega}{f v} = 0}
$


## Algorithms

### 2D

Laplacian equation on a square.

{% highlight cpp %}
// Parameters
func f = 1.;

// Mesh
int nn = 25;	//Mesh quality
mesh Th = square(nn, nn);

// Fespace
func Pk = P2;
fespace Uh(Th, Pk);
Uh u;

// Macro
macro grad(A) [dx(A), dy(A)] //

// Problem
varf vLaplacian (u, uh)
	= int2d(Th)(
		  grad(u)' * grad(uh)
	)
	- int2d(Th)(
		  f * uh
	)
	+ on(1, 2, 3, 4, u=0)
	;
matrix<real> Laplacian = vLaplacian(Uh, Uh, solver=sparsesolver);
real[int] LaplacianBoundary = vLaplacian(0, Uh);
u[] = Laplacian^-1 * LaplacianBoundary;

// Plot
plot(u, nbiso=30, fill=true, value=true, cmm="A");
{% endhighlight %}

|Result|
|--|
![Result]({{ site.url }}{{ site.baseurl }}/assets/Laplacian2D.png)|


### 3D

Laplacian equation on a cube.

{% highlight cpp %}
include "cube.idp"

// Parameters
func f = 1.;

// Mesh
int nn = 10;
mesh3 Th = cube(nn, nn, nn);

// Fespace
func Pk = P2;
fespace Uh(Th, Pk);
Uh u;

// Macro
macro grad(A) [dx(A), dy(A), dz(A)] //

// Problem
varf vLaplacian (u, uh)
	= int3d(Th)(
		  grad(u)' * grad(uh)
	)
	- int3d(Th)(
		  f * uh
	)
	+ on(1, 2, 3, 4, 5, 6, u=0)
	;
matrix<real> Laplacian = vLaplacian(Uh, Uh, solver=sparsesolver);
real[int] LaplacianBoundary = vLaplacian(0, Uh);
u[] = Laplacian^-1 * LaplacianBoundary;

// Plot
plot(u, nbiso=30, fill=true, value=true, cmm="A");
{% endhighlight %}

|Result|
|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/Laplacian3D.png)|
