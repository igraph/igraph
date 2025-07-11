const { exec } = require("child_process");

function getParameters(infomapBin) {
  return new Promise((resolve, reject) => {
    exec(`${infomapBin} --print-json-parameters`, (err, stdout, stderr) => {
      if (err) reject(err);

      if (stdout) resolve(JSON.parse(stdout));
      if (stderr) reject(stderr);
    });
  });
}

module.exports = getParameters;
