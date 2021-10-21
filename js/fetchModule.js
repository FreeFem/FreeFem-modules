window.addEventListener("popstate", function (event) {
  if (history.state) {
    fetchModule(event, history.state.url, history.state.title, true);
    event.preventDefault();
    event.stopPropagation();
  }
});

function fetchModule(e, url, title, disableHistory) {
  e.preventDefault();

  let previousModule = document.getElementsByClassName("current");
  if (previousModule.length > 0) previousModule[0].classList.remove("current");

  try {
    e.target.classList.add("current");
  } catch (error) {
    try {
      const links = document.getElementsByTagName("a");
      for (let link of links)
        if (link.href === url) link.classList.add("current");
    } catch (err) {
      console.error("fatality!");
      console.error(err);
    }
  }

  fetch(url)
    .then(function (response) {
      return response.text();
    })
    .then(function (html) {
      // Get content (not surrounding header, footer, nav, ...)
      var parser = new DOMParser();
      var content = parser.parseFromString(html, "text/html");
      var algoDisplay = content.getElementById("algoDisplay");

      // Replace in the page
      document.querySelector("#algoDisplay").innerHTML = algoDisplay.innerHTML;

      // Reset scroll
      document.body.scrollTop = 0;
      document.documentElement.scrollTop = 0;

      // Set history
      if (!disableHistory)
        history.pushState({ url: url, title: title }, title, url);
      document.title = "FreeFEM - " + title;

      // Relaunch MathJax
      MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
    })
    .catch(function (err) {
      console.log("Failed to fetch page: ", err);
    });
}
