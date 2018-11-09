showNav = () => {
  if (nav.style.display === 'block') {
    nav.style.display = 'none'
    if (navImage.src.includes('menu-close.svg')) {
      navImage.src = navImage.src.replace('menu-close.svg', 'menu.svg')
    }
  }
  else {
    nav.style.display = 'block'
    if (navImage.src.includes('menu.svg')) {
      navImage.src = navImage.src.replace('menu.svg', 'menu-close.svg')
    }
  }
}
