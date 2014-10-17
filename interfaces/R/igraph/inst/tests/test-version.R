
context("igraph_version")

test_that("igraph_version returns a version string", {

  ## This is essentially a semver regex, we do not allow a
  ## leading 'v' and space after
  regex <- paste0(
    "\\b",                                             # word boundary
    "(?:0|[1-9][0-9]*)\\.",                           # major
    "(?:0|[1-9][0-9]*)\\.",                           # minor
    "(?:0|[1-9][0-9]*)",                              # patch
    "(?:-[\\da-zA-Z\\-]+(?:\\.[\\da-zA-Z\\-]+)*)?",   # prerelease
    "(?:\\+[\\da-zA-Z\\-]+(?:\\.[\\da-zA-Z\\-]+)*)?", # word boundary
    "\\b"
  )

  expect_true(grepl(regex, igraph_version()))

})
