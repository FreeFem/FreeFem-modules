---
name: Wentzel
category: academic
layout: module
---

# Wentzel

Computation of Steklov or Wentzell eigenvalues.

Uses radial representation for the shape in terms of Fourier coefficients.

Parametrization uses a vector of the form.

`vec = [a0, as, bs]` where `as` are coefficients of `cos` and `bs` coefficients of `sin`

## Variational form

$
\displaystyle{
  \int_{\Omega}{\nabla\mathbf{u}:\nabla\mathbf{v} + \beta} + \beta\int_{\partial\Omega}{\nabla_{\tau}\mathbf{u}:\nabla_{\tau}\mathbf{v}} = \lambda\int_{\Omega}{\mathbf{u}\cdot\mathbf{v}}
}
$

## Algorithm

{% highlight cpp %}
// Parameters
real beta = 0.3; //parameter for the Laplace-Beltrami coefficient
				 //For Steklov eigenvalue choose 0
int M = 300;
int eigCount = 11; //the number of eigenvalues to compute

real[int] vec = [1.0,0,0,0,0.3];  //vector of Fourier coefficients for the
								  //boundary parametrization

// Function
//evaluation of a trigonometric polynomial
func real ptrig(real t, real[int] VV){
	real[int] vect = VV;
	int n = (vect.n-1)/2;
	real[int] as = vect(1:n) ;    //coeffs of cos
	real[int] bs = vect(n+1:2*n); //coeffs of sin
	real a0 =vect(0);
	real Sum = a0;   //initialize sum
	for(int i=0;i<n;i++){Sum+=(as[i]*cos((i+1)*t)+bs[i]*sin((i+1)*t));}
	return Sum;
}

// Mesh
border C(t=0,2*pi){x=cos(t)*ptrig(t,vec);
				  y=sin(t)*ptrig(t,vec);
				  label=1;}
mesh Th = buildmesh(C(200));
plot(Th);

// Fespace
fespace Vh(Th, P2);
Vh uh,vh;

// Problem
//weak form of the Wentzell problem
varf va (uh, vh)
	= int2d(Th)(
		  dx(uh)*dx(vh)
		+ dy(uh)*dy(vh)
	)
	+ int1d(Th, 1)(
		  beta*(
			  dx(uh)*dx(vh)
			- dx(uh)*N.x*(N.x*dx(vh)+N.y*dy(vh))
			- dx(vh)*N.x*(N.x*dx(uh)+N.y*dy(uh))
			+ N.x*(dx(vh)*N.x+dy(vh)*N.y)*N.x*(dx(uh)*N.x+dy(uh)*N.y)
			+ dy(uh)*dy(vh)
			- dy(uh)*N.y*(dx(vh)*N.x+dy(vh)*N.y)
			- dy(vh)*N.y*(dx(uh)*N.x+dy(uh)*N.y)
			+ (N.y)^2*(dx(vh)*N.x+dy(vh)*N.y)*(dx(uh)*N.x+dy(uh)*N.y)
		)
	)
	;

varf vb(uh, vh) = int1d(Th,1)( uh * vh); // linear form
matrix A = va(Vh, Vh ,solver = sparsesolver); // Matrix A on left side
matrix B = vb(Vh, Vh);						  // Matrix B on right side
real cpu = clock();  // get the clock in seconds

// Solve
real[int] ev(eigCount); // Holds Eigenvalues
Vh[int] eV(eigCount);   // Holds Eigenfunctions

int numEigs = EigenValue(A, B, sym=true, sigma=0, value=ev, vector=eV);
for(int i = 0; i < eigCount; i++) { //Plot the spectrum and show the eigenvalues
	// Plot
	plot(eV[i], fill=true, value=true, cmm= ev[i]);

	// Display
	cout << "Eigenvalue " << i << " :"  << ev[i] << endl;
}

cout << " Total CPU time = " << clock()-cpu << endl;
{% endhighlight %}

## References

[B. Bogosel, The method of fundamental solutions applied to boundary eigenvalue problems, 2016, Journal of Computational and Applied Mathematics](https://www.researchgate.net/publication/282658142_The_method_of_fundamental_solutions_applied_to_boundary_eigenvalue_problems)

### Authors

Author: [Beniamin Bogosel](http://www.cmap.polytechnique.fr/~beniamin.bogosel/)
