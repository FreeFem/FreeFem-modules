const onSearch = (event) => {
  let text = event.target.value;

  text = text.toLowerCase();

  for (let module of menu.children) {
    let name = module.children[0].innerHTML;
    name = name.toLowerCase();

    if (!name.includes(text)) module.style.display = "none";
    else module.style.display = null;
  }
};
