"use strict";

const gitRawCommits = require("git-raw-commits");
const conventionalCommitsParser = require("conventional-commits-parser");

function getCommits(from = "", to = "HEAD") {
  const commits = [];

  const gitOpts = { from, to, format: "%B%n-date-%n%aI" };

  return new Promise((resolve, reject) =>
    gitRawCommits(gitOpts)
      .pipe(conventionalCommitsParser())
      .on("data", (data) => commits.push(data))
      .on("finish", () => resolve(commits))
      .on("error", reject)
  );
}

module.exports = getCommits;
