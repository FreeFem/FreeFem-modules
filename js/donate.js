const donateLink = document.getElementById('donate')
donateLink.position = 'relative'

const donateMessage = document.createElement('div')
donateMessage.innerHTML = '<p>Coming soon</p>'
donateMessage.style.display = 'none'
donateMessage.style.backgroundColor = 'white'
donateMessage.style.color = 'black'
donateMessage.style.padding = '1rem'
donateMessage.style.position = 'absolute'
donateMessage.style.top = 0
donateMessage.style.right = 0

donateLink.appendChild(donateMessage)

donateLink.onclick = function() {
	donateMessage.style.display = 'block'
	setTimeout(function() {
		donateMessage.style.display = 'none'
	}, 3000)
}
