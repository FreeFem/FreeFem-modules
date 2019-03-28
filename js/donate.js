const donateLink = document.getElementById('donate')
const nav = donateLink.parentElement.querySelectorAll('nav > a')

donateLink.onclick = function() {
  donateLink.innerHTML = 'COMING SOON'
  if (!window.matchMedia("(max-width: 1200px)").matches) {
    donateLink.parentElement.querySelectorAll('nav > a').forEach(el => {
      el.style.marginLeft = '2rem'
    })
  }
  setTimeout(function() {
    donateLink.innerHTML = 'DONATE'
    if (!window.matchMedia("(max-width: 1200px)").matches) {
      donateLink.parentElement.querySelectorAll('nav > a').forEach(el => {
        el.style.marginLeft = ''
      })
    }
  }, 3000)
}
