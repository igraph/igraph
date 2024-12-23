function readVersion(contents) {
  const matches = contents.match(/INFOMAP\_VERSION = "(.*)";$/im);
  if (matches && matches[1]) {
    return matches[1];
  }
}

function writeVersion(contents, version) {
  return contents.replace(
    /INFOMAP\_VERSION = "(.*)";$/gim,
    `INFOMAP_VERSION = \"${version}\";`
  );
}

module.exports = { readVersion, writeVersion };
