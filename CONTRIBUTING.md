# Add a module

- Fork the project on your Github account

- Create a new file in `_modules/` in Markdown format, e.g `myModule.md` (you could use the [TEMPLATE](.github/TEMPLATE.md))

- Fill the front matter:

```
---
name: {required: module name}
category: {required: module category}
layout: module
---
```

- Write your module page

- Open a pull request

## Some tips

Create a link in markdown:
```markdown
[title][www.example.com]
```

Include an image:
```markdown
|Title|
|--|
|![not displayed title]({{ site.url }}{{ site.baseurl }}/assets/myimage.png)|
```
