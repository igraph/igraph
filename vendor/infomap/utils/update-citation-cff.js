function readVersion(contents) {
  const matches = contents.match(/^version: (.*)$/im);
  if (matches && matches[1]) {
    return matches[1];
  }
}

function writeVersion(contents, version) {
  contents = contents.replace(/version: (.*)$/gim, `version: ${version}`);
  contents = contents.replace(
    /^date-released: (.*)$/gim,
    "date-released: " + new Date().toLocaleDateString("sv-SE")
  );
  return contents;
}

module.exports = { readVersion, writeVersion };
