showNav = () => {
  if (nav.style.display === 'block') {
    nav.style.display = 'none'
    content.style.display = 'block'
    navImage.style.display = 'block'
    navImageClose.style.display = 'none'
  }
  else {
    nav.style.display = 'block'
    content.style.display = 'none'
    navImage.style.display = 'none'
    navImageClose.style.display = 'block'
  }
}
