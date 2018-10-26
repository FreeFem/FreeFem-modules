onSearch = (event) => {
  let text = event.target.value

  if (text.length === 0) {
    toogleCategories()
    return
  }

  text = text.toLowerCase()

  for (let i = 0; i < modulesList.children.length; i++) {
    const module = modulesList.children[i].children[0]
    let name = module.children[0].innerHTML
    name = name.toLowerCase()

    if (name.includes(text))
      module.style.display = 'flex'
    else
      module.style.display = 'none'
  }
}
