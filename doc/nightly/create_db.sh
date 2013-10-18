#! /bin/sh

rm -f nightly-test.db
echo "
CREATE TABLE 'files' (
       filename TEXT,
       type     TEXT,
       version  TEXT,
       branch   TEXT,
       hash     TEXT,
       date     TEXT,
       size     INTEGER,
       count    INTEGER DEFAULT 0,
       PRIMARY KEY(filename)
);
CREATE TABLE 'tests' (
       filename TEXT,
       result   TEXT,
       PRIMARY KEY(filename, result),
       FOREIGN KEY(filename) REFERENCES files(filename)
);
" | sqlite3 nightly-test.db 
