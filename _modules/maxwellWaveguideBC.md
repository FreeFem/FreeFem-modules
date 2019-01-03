---
name: "maxwellWaveguideBC"
category: "electromagnetism, academic"
layout: module
---

# maxwellWaveguideBC

This module illustrates how to implement waveguide boundary conditions in Maxwell equations. So far it is done for 2D planar case in Cartezian and Axial (m = 0) cases. Waveguides modes and waveguides port are hardcoded. 

## Problem
The module uses Maxwell equations written for $\mathbf{E]$ field in $\exp(-i \omega t)$ presentation: $curl curl \mathbf{E} = \frac{\omega^2}{c^2} \epsilon(x,y) \mathbf{E}$ (here: $\omega$ - is the frequency, $c$ - the speed of light, $\epsilon$ - dielectric permittivity; the CGS system of units is used) 
This module illustrates how to set boundary condition so that there is two waveguide ports through one of which an input signal enters. That is it solves a simple waveguide scattering problem

## Variational form
Variational form is 
$$\int_V dV ( curl(\mathbf{E}) curl(\mathbf{V}) - \frac{\omega^2}{c^2} \epsilon(x,y) \mathbf{E} \mathbf{V}) + \int_{dV} dS curl (\mathbf{E}) \mathbf{V} = 0$$
For waveguide ports the last term in the left hand side is written as
$$\int_{dS_{port}} curl \mathbf{E} \mathbf{V} =  \sum_{s=1,N_{modes}} C_s \int_{dS} curl \mathbf{E_s} \mathbf{V}$$
Here C_s is the amplitude of mode $s$ and $\mathbf(E_s)$ - the s-th mode's electric field. The $C_s$ in its turn is written as 
$$C_s = \frac{1}{N_s} \int_{S_{port}} dS ( \mathbf{E} \times \mathbf{H_s} - \mathbf{E_s} \times \mathbf{H})$$, here $\mathbf{H}$ and $\mathbf{H_s}$ are current and mode's magnetic field which are found from Maxwell equations, $N_s$ is the norm. 

Finally, the last term of the LHS side can be written in term of the direct product as
$$
\int_{dS_{port}} curl \mathbf{E} \mathbf{V} =  \sum_{s=1,N_{modes}} \frac{1}{N_s} \int_{S_{port}} dS ( \mathbf{E} \times \mathbf{H_s} - \mathbf{E_s} \times \mathbf{H}) \bigotimes \int_{dS} curl \mathbf{E_s} \mathbf{V},
$$
that allows one to implement it by means of FreeFem++ tools without low-level matrix assembly (for the price of performance, of course).

## Algorithms
### 2D - Cartesian

### 2D - axial-symmetric
```C
load "Element_Mixte"
real height = 3;
real length = 3;
func mesh generateMesh()
{
    return square(30,30, [x*length, y*height]);
}
func complex epsilon(real xx, real yy)
{
    return 1. + 0*1i;
//  if(xx > 0.25 && xx < 0.75/* && yy > 0.25 && yy < 0.75*/) return 5; else return  1.;
}

mesh Th = generateMesh();
//	plot(Th);
include "axial2D_modes.idp"

fespace Vh(Th,[RT1Ortho, P2]); // finite element space
	
func matrix<complex> MakswellEqs(real sigma, int m, real h)
{
    macro odx(F) (dx(F) +  1i*h*F) // dx -> odx 
    varf MKSWL([Hz,Hr, Hphi],[vHz,vHr, vHphi]) = 
         int2d(Th)(        
        	   (
		       (dy(Hphi) - 1i*m*Hr )*conj(dy(vHphi)- 1i*m*vHr )/y      //curl_z // vHphi => (Hphi*y)
        	     + (odx(Hr) - dy(Hz))*conj(odx(vHr) - dy(vHz))*y	       //curl_phi
        	     + (1i*m*Hz  - odx(Hphi))*conj(1i*m*vHz  - odx(vHphi))/y   //curl_r
        	   )	
        	   - sigma*y*epsilon(x,y)*(conj(vHr)*Hr+conj(vHphi)*Hphi/(y*y)+conj(vHz)*Hz)  // vHphi => (Hphi*y)
        	  )        
        	  + on(3,  Hz = 0, Hr = 0, Hphi = 0)
        	  + on(1,  Hphi = 0)        //Hr - at the axis does'not have meaning due to used RT1Ortho
                  ;
    matrix<complex> OP = MKSWL(Vh, Vh, solver=UMFPACK); 
    return OP;
}

func complex[int] vectorBCh(int modeNum, int boundaryLabel, real k)
{
    int m = modeM [abs(modeNum)];
    varf BoundCondRot([Ex,Ey,Ez],[vEx,vEy,vEz]) =   
        int1d(Th,boundaryLabel)
            (/*2.*pi*//*y*/(                          // x -> z, y ->  r, z -> phi
                anH(x, y, (N.x)*modeNum, 1, height, 0, k)*conj(/*0.*/ - N.y*vEz/*/y*/) // 0 => N.z
              + anH(x, y, (N.x)*modeNum, 2, height, 0, k)*conj(N.x*vEz/*/y*/ /*- 0.*/)
              + anH(x, y, (N.x)*modeNum, 3, height, 0, k)*conj(N.y*vEx - N.x*vEy)*y
	     )*(1i*k));                  // [n v*] curl E = [n v*] (i k H) 
     complex[int] BCr = BoundCondRot(0,Vh);
     return BCr;
}	
	
func complex[int] vectorBCe(int modeNum, int boundaryLabel, real k) 
{
    int m = modeM [abs(modeNum)];
    cout<<"modeNum = "<<modeNum<<endl;
    varf BoundCondField([Ex,Ey,Ez],[vEx,vEy,vEz]) =  
        int1d(Th,boundaryLabel)
            ((// x(1) -> z, y(2) ->  r, z(3) -> phi;    n [Es Hv]
              2.*pi*/*y*/
	          conj( - (vEz/*/y*/)*anH(x, y, (N.x)*modeNum, 2, height, 0, k)
			+ (anE(x, y,  (N.x)*modeNum, 3, height, 0, k))*((-dx(vEz)/*/y*/ +  1i*m*vEx/*/y*/) / (1i*k) ) //Ez = anE(..,3,..)*y
			+ (vEy) *anH(x, y,  (N.x)*modeNum, 3, height, 0, k) * y  
			- anE(x, y,  (N.x)*modeNum, 2, height, 0, k)*((dx(vEy) - dy(vEx))/(1i*k) * y)
		      )*(N.x)
	     )/fNorm(modeNum, height, k)
	    ); 
    complex[int] BCrf = BoundCondField(0,Vh);
    return BCrf;
}

func matrix<complex> matrixBC(int modeNum, int boundaryLabel, real k)
{		
   complex[int] BCr  = vectorBCh( modeNum, boundaryLabel,k);
   complex[int] BCrf = vectorBCe( -modeNum, boundaryLabel,k);
   matrix<complex> Br = BCr*BCrf';                           // cartesian product
   return Br;		
}	
	
func matrix<complex> generateBCMatrix(int numberOfModes, int boundaryLabel, int sign, real k)
{
    matrix<complex> Br1; 
    for(int i = 1; i <= numberOfModes; i++)
    {
         matrix<complex> Br2 = matrixBC(sign*i, boundaryLabel, k);
         Br1 = Br1 + Br2;
    }
    return Br1;
}	

//complimentary function for calculating scattering
func complex[int] modeAmplitudes  
(complex[int]& eVec, int numOfModes, int referenceMode, int inputBoundaryLabel, int otherBoundaryLabel, real k,  bool needdB)
{
    complex[int] result(numOfModes);
    Vh<complex> [ex,ey,ez]; 
    ex[] = eVec;
    if(needdB) for(int i = 0; i < numOfModes; i++) result(i) = log10(0.);
    complex[int] probe = vectorBCe(-referenceMode, inputBoundaryLabel, k);
    Vh<complex> [pex, pey, pez];
    pex[] = probe;	
    int m = modeM [abs(referenceMode)-1];
    complex power =  // eVec'*probe;
        int1d(Th,otherBoundaryLabel)((// x -> z, y ->  r, z -> phi;    n [Es Hv]
	  2.*pi*y*
	  conj(-(pez/y)*(-(dx(ex)/y + 1i*m*ex/y) / (1i*k) )       
	  + ((ez/y))*(-(dx(pex)/y + 1i*m*pex/y) / (1i*k) ) //Ez = anE(..,3,..)*y
          + (pey) *((dx(ey) - dy(ex))/(1i*k))
          - (ey)*((dx(pey) - dy(pex))/(1i*k))
          )*(N.x)
	 )
     /* / fNorm(referenceMode, height, k)*/);	
    cout<<"modeAmplitudes: mode = "<<referenceMode<<"   power = "<<abs(power)<<"(abs),   "<<power<<endl;
    for(int testmode = 0; testmode < numOfModes; testmode++)	
    {
      // complex[int] probe2 = vectorBCe( (testmode+1), otherBoundaryLabel);
      // complex res = eVec'*probe2/power;
         int m = modeM [testmode];
         int modeNum = -(testmode + 1);
	 complex res = int1d(Th,otherBoundaryLabel)((// x -> z, y ->  r, z -> phi;    n [Es Hv]
			     2.*pi*y*
			      conj(
				  - (ez/y)*anH(x, y, (N.x)*modeNum, 2, height, 0, k)  
                                  + (anE(x, y,  (N.x)*modeNum, 3, height, 0, k))*(-(dx(ez)/y + 1i*m*ex/y) / (1i*k) ) //Ez = anE(..,3,..)*y
				  + (ey) *anH(x, y,  (N.x)*modeNum, 3, height, 0, k)
				  - anE(x, y,  (N.x)*modeNum, 2, height, 0, k)*((dx(ey) - dy(ex))/(1i*k))
														     )*(N.x)
			          )
														 /(fNorm(modeNum, height, k)/*power*/));	
            cout << "inner loop: res=" << res <<",  abs(res) = "<<abs(res)<<endl; 														 
        if(needdB) result(testmode) =  20*log10(abs(res)); 
        else       result(testmode) =  res;
    }
	
    return result;
}
		
	int m = 0;                           //not used in this implementation, keep it for consistency 
	real sigma = 9.*pi*pi-75;//70.;
	real k = sqrt(sigma);
	cout<<"h1 = "<<H(k, Kappa(1,height,1))<<",  kappa1 = "<<Kappa(1,height,1)<<",  k = "<<k<<",  beta_ph = "<<k/H(k, Kappa(1,height,1))<<endl; 
	cout<<"h2 = "<<H(k, Kappa(2,height,1))<<",  kappa2 = "<<Kappa(2,height,1)<<",  k = "<<k<<",  beta_ph = "<<k/H(k, Kappa(1,height,1))<<endl; 
	cout<<"h3 = "<<H(k, Kappa(3,height,1))<<",  kappa3 = "<<Kappa(3,height,1)<<",  k = "<<k<<",  beta_ph = "<<k/H(k, Kappa(1,height,1))<<endl; 
	cout<<"h3 = "<<H(k, Kappa(4,height,1))<<",  kappa3 = "<<Kappa(4,height,1)<<",  k = "<<k<<",  beta_ph = "<<k/H(k, Kappa(1,height,1))<<endl; 	

int rightBorderLabel = 2;
int  leftBorderLabel = 4;
int numberOfModesToAccount = 3;
int targetMode = 3;
matrix<complex> Br = generateBCMatrix(numberOfModesToAccount, rightBorderLabel,1, k);
matrix<complex> Bl = generateBCMatrix(numberOfModesToAccount, leftBorderLabel, 1, k);
matrix<complex> OP = MakswellEqs(sigma, m, 0) ;
	
matrix<complex> A; //A =  OP + Bl + Br;   //Here is an odd thing happens. Initially, I wrote  "A = OP + generateBCMatrix(...) + generateBCMatrix(...);"
A = OP;                                   //but I've found out that the result depends on the order of calling the constructors (and,generally, is wrong),
A = A + Br;                               //so now I'm not sure now even if  "matrix A= OP + Bl +Br"  will work correctly and add them step by step
A = A + Bl; 
complex[int] u(Br.n);
set(A, solver = UMFPACK);
complex[int] RHS =  vectorBCh(targetMode , leftBorderLabel, k);	
    
u = A^-1*(RHS);
	

Vh<complex> [Ex,Ey,Ez]; 
Ex[] = u;

plot(Ex, fill = true, cmm = "Ex");
plot(Ey, fill = true, cmm = "Ey");
//fespace Uh (Th, P1); Uh<complex> Ezreal = Ez/(y+1e-10);
//plot(Ezreal, fill = true, cmm = "Ez");
	
cout<<	modeAmplitudes(u, numberOfModesToAccount, -targetMode,  leftBorderLabel,  rightBorderLabel, k, false)<<endl;
cout<<	modeAmplitudes(u, numberOfModesToAccount,  targetMode,  leftBorderLabel,  leftBorderLabel, k, false)<<endl; 
```

_[Optional]_
### 2D

[Case description]

### 3D

[Case description]

### Optional

[e.g., Gmsh script]

_[End optional]_

## Validation

[Validation algorithms and results, e.g. convergence curves]

## References

[Links to references, in open-access if possible]

## Authors

Petr Makhalov
