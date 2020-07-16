---
name: 3D meshing example
category: academic, electromagnetism
layout: module
---

# TEAM 7: Asymmetrical Conductor with a Hole

Build a mesh for the TEAM 7 problem and compute a div-free current in the coil.

## Algorithms

### 3D

TEAM 7 model and a Dirichlet-Laplace problem to get a div-free current, without geometry modifications.

{% highlight cpp %}
load "msh3"
load "tetgen"
load "medit"
load "iovtk"

// 2D INITIAL MESH RESOLUTION --------------------------------------------------
int n = 50;

// ALUMINIUM PLATE -------------------------------------------------------------
// Outer contour
border a01(t = 0.000, 0.294){x = t    ; y = 0.000; label = 1;}
border a02(t = 0.000, 0.294){x = 0.294; y = t    ; label = 1;}
border a03(t = 0.294, 0.000){x = t    ; y = 0.294; label = 1;}
border a04(t = 0.294, 0.000){x = 0.000; y = t    ; label = 1;}

// Interior contour
border a05(t = 0.018, 0.126){x = t    ; y = 0.018; label = 1;}
border a06(t = 0.018, 0.126){x = 0.126; y = t    ; label = 1;}
border a07(t = 0.126, 0.018){x = t    ; y = 0.126; label = 1;}
border a08(t = 0.126, 0.018){x = 0.018; y = t    ; label = 1;}

// Two-dimensional mesh and refine
mesh Tal = buildmesh(a01(n) + a02(n) + a03(n) + a04(n) + a05(-n/2) + a06(-n/2) + a07(-n/2) + a08(-n/2));
for (int cnt = 1; cnt < 5; cnt++){
  Tal = adaptmesh(Tal, 0.0025, IsMetric=true, nbjacoby=10, nbsmooth=10, nbvx=1000000);
}

// Three dimensional volume and surface mesh
mesh3 TalV = buildlayers(Tal, 5, zbound=[0, 0.019], transfo=[x, y, z]);
TalV = buildBdMesh(TalV);
meshS TalS = TalV.Gamma;
medit("Aluminium plate", TalV);

// COPPER COIL -----------------------------------------------------------------
// Outer contour
border a09(t = 0.144 , 0.244 ){x = t                  ; y = 0.000             ; label = 2;}
border a10(t = 1.5*pi, 2.0*pi){x = 0.244 + 0.05*cos(t); y = 0.05 + 0.05*sin(t); label = 2;}
border a11(t = 0.05  , 0.15  ){x = 0.294              ; y = t                 ; label = 2;}
border a12(t = 0.0*pi, 0.5*pi){x = 0.244 + 0.05*cos(t); y = 0.15 + 0.05*sin(t); label = 2;}
border a13(t = 0.244 , 0.144 ){x = t                  ; y = 0.200             ; label = 2;}
border a14(t = 0.5*pi, 1.0*pi){x = 0.144 + 0.05*cos(t); y = 0.15 + 0.05*sin(t); label = 2;}
border a15(t = 0.150 , 0.050 ){x = 0.094              ; y = t                 ; label = 2;}
border a16(t = 1.0*pi, 1.5*pi){x = 0.144 + 0.05*cos(t); y = 0.05 + 0.05*sin(t); label = 2;}

// Interior contour
border a17(t = 0.144 , 0.244 ){x = t                   ; y = 0.025              ; label = 3;}
border a18(t = 1.5*pi, 2.0*pi){x = 0.244 + 0.025*cos(t); y = 0.05 + 0.025*sin(t); label = 3;}
border a19(t = 0.05  , 0.15  ){x = 0.269               ; y = t                  ; label = 3;}
border a20(t = 0.0*pi, 0.5*pi){x = 0.244 + 0.025*cos(t); y = 0.15 + 0.025*sin(t); label = 3;}
border a21(t = 0.244 , 0.144 ){x = t                   ; y = 0.175              ; label = 3;}
border a22(t = 0.5*pi, 1.0*pi){x = 0.144 + 0.025*cos(t); y = 0.15 + 0.025*sin(t); label = 3;}
border a23(t = 0.150 , 0.05 ){x = 0.119                ; y = t                  ; label = 3;}
border a24(t = 1.0*pi, 1.5*pi){x = 0.144 + 0.025*cos(t); y = 0.05 + 0.025*sin(t); label = 3;}

// Two-dimensional mesh and refine
mesh Tcu = buildmesh(a09(n) + a10(n) + a11(n) + a12(n) + a13(n) + a14(n) + a15(n) + a16(n) +
                     a17(-n/2) + a18(-n/2) + a19(-n/2) + a20(-n/2) + a21(-n/2) + a22(-n/2) + a23(-n/2) + a24(-n/2));
for (int cnt = 1; cnt < 5; cnt++){
  Tcu = adaptmesh(Tcu, 0.005, IsMetric=true, nbjacoby=10, nbsmooth=10, nbvx=1000000);
}

// Three dimensional volume and surface mesh
int[int] top1 = [0, 4], down1 = [0, 5];
mesh3 TcuV = buildlayers(Tcu, 20, zbound=[0.049, 0.149], transfo=[x, y, z], labelup=top1, labeldown=down1);
TcuV = buildBdMesh(TcuV);
meshS TcuS = TcuV.Gamma;
medit("Copper coil", TcuV);

// BOUNDING BOX ----------------------------------------------------------------
int[int] labs = [6, 6, 6, 6];
mesh Tvd = square(n/2, n/2, [-1.353 + 3.0*x, -1.353 + 3.0*y], label=labs);
int[int] top2 = [0, 6], down2 = [0, 6];
mesh3 TvdV = buildlayers(Tvd, 3, zbound=[-0.3, 0.449], transfo=[x, y, z], labelup=top2, labeldown=down2);
TvdV = buildBdMesh(TvdV);
meshS TvdS = TvdV.Gamma;
medit("Bounding box", TvdV);

// GLUE SURFACE MESHES, DEFINE DOMAIN INDICATORS, AND VOLUME CONSTRAINTS -------
meshS TS = TalS + TcuS + TvdS;
int alp = 1; // aluminium plate domain indicator
int cuc = 2; // copper coil domain indicator
int air = 3; // air/void space domain indicator
real[int] domain = [0.15, 0.15, 0.01, alp, 0.0005, 0.15, 0.01, 0.1, cuc, 0.001, -1.0, -1.0, 0.0, air, 0.5];
mesh3 TV = tetg(TS, switch = "pqaAAYYQ", nbofregions = 3, regionlist = domain);
medit("Global mesh", TV);

// FINITE ELEMENT SPACES -------------------------------------------------------
fespace S0(TV, P0);
fespace S1(TV, P1);
fespace V0(TV, [P0, P0, P0]);

// DOMAIN INDICATOR FEM FUNCTION -----------------------------------------------
S0 inCoil = 1e-12*(region == alp) + 1.0*(region == cuc) + 1e-12*(region == air);

// PROBLEM TO GET DIV FREE CURRENT ---------------------------------------------
S1 u, v;
solve Current(u, v)
      = int3d(TV)( inCoil * [dx(u), dy(u), dz(u)]' * [dx(v), dy(v), dz(v)] )
      + on(2, u=1.0)  // outer coil surface
      + on(3, u=0.0); // inner coil surface

// CURRENT ---------------------------------------------------------------------
V0 [Jx, Jy, Jz] = inCoil*[dy(u), -dx(u), dz(u)];
medit("Current", TV, [Jx, Jy, Jz]);
savevtk("currentCoilVector.vtu", TV, [Jx, Jy, Jz]);
{% endhighlight %}

|Geometry|
|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/coilGeometry.jpg)|

|Result|
|--|
|![Result]({{ site.url }}{{ site.baseurl }}/assets/coilSolution.jpg)|

## Authors

Author: [Fotios Kasolis](mailto:fotios.kasolis@gmail.com)
