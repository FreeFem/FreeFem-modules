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
  setTimeout(function() {
    button.style.background = 'none'
  }, 1000)
}

for (let i = 0; i < codeBlocks.length; i++) {
  // Copy / Paste button
  var p = document.createElement('p')
  var node = document.createTextNode('copy')
  p.appendChild(node);

  let button = document.createElement('button')
  button.className = 'copyPasteButton'
  button.appendChild(p)
  button.onclick = function() { copyPasteCode(codeBlocks[i], button) }

  // Add button
  codeBlocks[i].appendChild(button)
}
