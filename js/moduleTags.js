tags = []

addToTags = (item) => {
  const items = item.split(', ')

  items.forEach(item => {
    const index = tags.find(x => x.value === item)
    if (index)
    	return

    tags.push({
		value: item,
		label: item.charAt(0).toUpperCase() + item.substring(1),
		selected: false,
		disabled: false
	})
  })
}

includesOnceOf = (string, array) => {
	const length = array.length

	for (let i = 0; i < length; i++) {
		if (string.includes(array[i])) {
			return true
		}
	}

	return false
}

toggleTags = (value) => {
	const usedTags = []

	for (let i = 0; i < value.target.children.length; i++) {
		usedTags.push(value.target.children[i].value)
	}

	const menu = document.getElementById('menu')

	if (usedTags.length === 0) {
		for (let i = 0; i < menu.children.length; i++) {
			const child = menu.children[i]
			child.classList.remove('hidden')
		}
		return
	}

	for (let i = 0; i < menu.children.length; i++) {
		const child = menu.children[i]
		const moduleTags = child.children[1].innerHTML

		if (includesOnceOf(moduleTags, usedTags)) {
			child.classList.remove('hidden')
		} else {
			child.classList.add('hidden')
		}
	}
}
