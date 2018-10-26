# Add a module

- Fork the project on your Github account

- Create a new file in `_modules/`, e.g `myModule.js`

- Fill the front matter:

```
---
name: {required: module name}
category: {required: module category}
layout: module
author: {required: your name}
authorURL: {your personnal page/Github profil}
---
```

- Write your module page

Page architecture :

```
# {Module name}

{Short description}

## Problem

{Continuous equations, with parameters description, ...}

## Variational form

{Variational form and spaces}

## Algorithms

[Optionaly 2D/3D]

{Description of the case}

{FreeFem++ algorithm}

## Validation

[Optional]

## References

[Optional]

{Publications related to the algorithm, prefer researchgate if possible}
```

- Open a pull request