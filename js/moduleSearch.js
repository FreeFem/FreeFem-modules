onSearch = (event) => {
  let text = event.target.value

/*
  if (text.length === 0) {
    toggleCategories()
    return
  }
*/
  text = text.toLowerCase()

  for (let i = 0; i < menu.children.length; i++) {
    const module = menu.children[i]
    let name = module.children[0].innerHTML
    name = name.toLowerCase()

    if (!name.includes(text))
      module.style.display = 'none'
    else
      module.style.display = null;
  }
}
