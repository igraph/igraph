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
       filename   TEXT,
       platform   TEXT,
       test       TEXT,
       result     TEXT,
       resultcode INTEGER,
       url        TEXT,
       PRIMARY KEY(filename, platform, test),
       FOREIGN KEY(filename) REFERENCES files(filename)
);
" | sqlite3 nightly-test.db 
