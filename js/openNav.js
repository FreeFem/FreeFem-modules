let opened = false;

function openNav() {
  var nav = document.getElementById('topNav');
  if (opened) {
    nav.style.width = 0
    opened = false
  } else {
    nav.style.width = "100%"
    opened = true
  }
}
