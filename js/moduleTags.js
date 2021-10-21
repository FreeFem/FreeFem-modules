const tags = [];

const addToTags = (itemsList) => {
  const items = itemsList.split(", ");

  items.forEach((item) => {
    const index = tags.find((x) => x.value === item);
    if (index) return;

    tags.push({
      value: item,
      label: item.charAt(0).toUpperCase() + item.substring(1),
      selected: false,
      disabled: false,
    });
  });
};

const includesOnceOf = (string, array) => {
  const length = array.length;

  for (let i = 0; i < length; i++) {
    if (string.includes(array[i])) {
      return true;
    }
  }

  return false;
};

const toggleTags = (value) => {
  const usedTags = [];

  for (let child of value.target.children) {
    usedTags.push(child.value);
  }

  const menu = document.getElementById("menu");

  if (usedTags.length === 0) {
    for (let child of menu.children) {
      child.classList.remove("hidden");
    }
    return;
  }

  for (let child of menu.children) {
    const moduleTags = child.children[1].innerHTML;

    if (includesOnceOf(moduleTags, usedTags)) {
      child.classList.remove("hidden");
    } else {
      child.classList.add("hidden");
    }
  }
};
