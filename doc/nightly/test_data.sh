#! /bin/sh

echo "
INSERT INTO 'files'
       VALUES ('c/igraph-0.7-pre+41.badcafe.tar.gz', 'c', '0.7-pre', 'develop',
               'f5115b6', '2013-09-04', '1200000', 0);
INSERT INTO 'files'
       VALUES ('c/igraph-0.6.5.tar.gz.', 'r', '0.6.5', 'master',
               'badcafe', '2013-09-03', '1200000', 0);
INSERT INTO 'files'
       VALUES ('python/python-igraph-0.6.5.tar.gz', 'python', '0.6.5', 'master',
               'f5115bb', '2013-08-05', '1200000', 0);
INSERT INTO 'files'
       VALUES ('c/igraph-0.7-pre+25.badcafe.tar.gz', 'c', '0.7', 'develop',
               'badcafe', '2013-09-02', '1200000', 0);
INSERT INTO 'files'
       VALUES ('r/igraph-0.7-pre+14.badcafe.tar.gz', 'r', '0.7', 'develop',
               'badcafe', '2013-09-01', '1200000', 0);
" | sqlite3 nightly.db
