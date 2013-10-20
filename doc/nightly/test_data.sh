#! /bin/sh

echo "
INSERT INTO 'files'
       VALUES ('c/igraph-0.6.6-pre+41.849c191.tar.gz', 'c', '0.6.6-pre',
               'release-0.6', '849c191', '2013-10-09', 2886425, 0);
INSERT INTO 'files'
       VALUES ('r/igraph_0.6.5.999-41.tar.gz', 'r', '0.6.6-pre',
               'release-0.6', '849c191', '2013-10-08', 2316191, 0);
INSERT INTO 'files'
       VALUES ('r/igraph_0.6.999-523.tar.gz', 'r', '0.7.0-pre',
               'develop', '50c1ec0', '2013-10-08', 2367155, 0);

INSERT INTO 'tests'
       VALUES ('r/igraph_0.6.999-523.tar.gz', 'r-3.0.2--osx-10.8.5--gcc',
               'R CMD check', 'WARNING', 2,
               'http://pave.igraph.org/job/igraph-r-test/5/console');
INSERT INTO 'tests'
       VALUES ('c/igraph-0.6.6-pre+41.849c191.tar.gz', 'osx-10.8.5--gcc',
               'make check', 'ERROR', 3, NULL);
" | sqlite3 nightly-test.db
