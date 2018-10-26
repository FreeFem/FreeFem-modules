onSearch = (event) => {
  let text = event.target.value

  if (text.length === 0) {
    toogleCategories()
    return
  }

  text = text.toLowerCase()

  for (let i = 0; i < menu.children.length; i++) {
    const module = menu.children[i]
    let name = module.children[0].innerHTML
    name = name.toLowerCase()

    if (name.includes(text) || name.includes("Home"))
      module.style.display = 'flex'
    else
      module.style.display = 'none'
  }
}
