---
name: Piezoelectricity
category: solid
layout: module
---

# Elasticity

Algorithms for solving the axisymmetric linear piezoelectricity equations

## Problem

Solve:

$
\displaystyle{

}
$

With:

$
\displaystyle{

}
$

and

$
\displaystyle{
}
$

## Variational form

Let . The variational form reads as follows:

$
\displaystyle{

}
$

## Algorithms

### 2D

Piezoelectricity equation on a circular disc.

{% highlight cpp %}

{% endhighlight %}

|Result warped by a factor 1000|
|--|
|![Result warped by a factor 1000]({{ site.url }}{{ site.baseurl }}/assets/Piezoelectricity2D.png)|

### 3D

Piezoelectricity equation on a circular disc.

{% highlight cpp %}


{% endhighlight %}

|Result warped by a factor 1000|
|--|
|![Result warped by a factor 1000]({{ site.url }}{{ site.baseurl }}/assets/Elasticity3D.png)|

### Optional

Gmsh script:

{% highlight cpp %}



{% endhighlight %}

|Mesh|
|--|
|![Mesh]({{ site.url }}{{ site.baseurl }}/assets/Elasticity3DMesh.png)|

## Validation

TODO

## Authors

Author: [Marek Moszyński](https://github.com/marmoszy)
