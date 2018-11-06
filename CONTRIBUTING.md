# How to add a module

- Fork the project on your Github account

- Create a new file in `_modules/` in Markdown format, e.g `myModule.md`. Please use the [MODULE_TEMPLATE.md](MODULE_TEMPLATE.md) file.

- Fill the header with your module name and the corresponding tag :

```
---
name: "module name"
category: "module category (solid, electromagnetism, academic, fluid, etc.)"
layout: module
---
```
Note: we do not support multiple tags yet. Let us know if you require this functionality.

- Write your module page

- Open a pull request

## Markdown tips

To create a link in markdown:
```markdown
[title][www.example.com]
```

To include an image, add the compressed PNG or JPG to the [assets](/assets) folder and use the following:
```markdown
|Title|
|--|
|![alt attribute - not displayed, usually the same as the title]({{ site.url }}{{ site.baseurl }}/assets/myimage.png)|
```
