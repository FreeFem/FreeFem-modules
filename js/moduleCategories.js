categories = []

addToCategories = (item) => {
  const index = categories.indexOf(item)
  if (index > -1)
    return

  addCategory(item)
  categories.push(item)
}

addCategory = (item) => {
  const div = document.createElement('div')
  div.className = 'custom-control custom-checkbox light-' + categories.length

  const label = document.createElement('label')
  label.className = 'custom-control-label'
  label.htmlFor = item
  label.innerHTML = item

  const input = document.createElement('input')
  input.className = 'custom-control-input'
  input.type = 'checkbox'
  input.checked = true
  input.id = item
  input.onchange = function() { onCategoryChange(item, input.checked) }

  div.appendChild(input)
  div.appendChild(label)

  categoriesDiv.appendChild(div)
}

onCategoryChange = (category, checked) => {
  if (checked) {
    categories.push(category)
  } else {
    const index = categories.indexOf(category)
    if (index > -1)
      categories.splice(index, 1)
  }

  toogleCategories()
}

toogleCategories = () => {
  for (let i = 0; i < menu.children.length; i++) {
    const module = menu.children[i]
    if (module.children[0].innerHTML.includes("Home")) {
      module.style.display = 'flex'
      continue
    }

    let category = module.children[1].innerHTML
    category = category.toLowerCase()

    let show = false
    for (let j = 0; j < categories.length; j++) {
      if (category.includes(categories[j].toLowerCase()))
        show = true
    }

    if (show)
      module.style.display = 'flex'
    else
      module.style.display = 'none'
  }
}
