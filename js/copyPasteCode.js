let codeBlocks = document.getElementsByClassName('highlight')

function copyPasteCode(code, button) {
  // Get text content (without html char and span)
  const content = code.children[0].children[0]
  let contentString = content.textContent
  contentString = contentString.replace(/<span class="[A-Za-z0-9]*">/g, '')
  contentString = contentString.replace(/<\/span>/g, '')

  // Put text in textarea (select undefined in code tag)
  let textArea = document.createElement('textarea')
  textArea.innerHTML = contentString
  content.appendChild(textArea)

  // Select and copy
  textArea.select()
  document.execCommand('copy')

  // Clear
  content.removeChild(textArea)

  // Color button
  button.style.backgroundColor = 'rgba(0, 255, 0, 0.5)'
  button.style.color = 'black'
  button.innerText = 'copied!'
  button.onmouseleave = function () {
    button.innerText = 'copy'
    button.style.backgroundColor = '#084887'
    button.style.color = '#f9ab55'
  }
  // setTimeout(function() {
  //   button.style.background = 'none'
  // }, 1000)
}

for (let i = 0; i < codeBlocks.length; i++) {
  // Copy / Paste button
  let button = document.createElement('button')
  button.className = 'copyPasteButton'
  button.innerText ='copy'
  button.onclick = function() { copyPasteCode(codeBlocks[i], button) }

  // Add button
  codeBlocks[i].appendChild(button)
}
