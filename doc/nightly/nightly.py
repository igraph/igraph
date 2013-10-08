#! /usr/bin/env/python

#
#
# TO create the DB, run:
"""
rm -f nightly.db
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
" | sqlite3 nightly.db 
"""

# Some example data:
"""
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

"""

import bottle
import bottle_sqlite
import socket

plugin = bottle_sqlite.Plugin(dbfile='nightly.db')
nightly=bottle.Bottle()
nightly.install(plugin)

urlmap={ 'c': 'C library', 'r': 'R package', 'python': 'Python extension',
         'msvc': 'C library for MSVC' }

# This is the main web page, you can choose your download here
# It can be filtered
@nightly.route("/")
@nightly.route("/list/")
@nightly.route("/list")
@nightly.route("/list/<dtype>")
@nightly.route("/list/<dtype>/<version>")
@nightly.route("/list/<dtype>/<version>/<branch>")
def list_files(db, dtype="all", version="all", branch="all"):

    files = db.execute("SELECT * FROM files    \
                        WHERE type LIKE ? AND version LIKE ? AND branch LIKE ? \
                        ORDER BY version DESC, date DESC",
                       ("%" if dtype=="all" else dtype,
                        "%" if version=="all" else version,
                        "%" if branch=="all" else branch)).fetchall()

    versions=db.execute("SELECT DISTINCT version FROM files \
                         ORDER BY version DESC").fetchall()
    versions=[ e[0] for e in versions ]

    types=db.execute("SELECT DISTINCT type FROM files \
                      ORDER BY type").fetchall()
    types=[ e[0] for e in types ]

    branches=db.execute("SELECT DISTINCT branch FROM files \
                         ORDER BY branch").fetchall()
    branches=[ e[0] for e in branches ]

    return bottle.template('main', files=files, versions=versions,
                           types=types, branches=branches, urlmap=urlmap,
                           dtype=dtype, branch=branch, version=version)

@nightly.route("/get/<filename>")
def get_file(db, filename):
    db.execute("UPDATE files SET count=count+1 WHERE filename=?", \
               (filename,))
    return bottle.static_file(filename, ".")

@nightly.error(404)
def error404(error):
    return "Page does not exist, maybe your syntax is wrong"

myname=socket.gethostname()
if myname[-4] == ".com":
    bottle.run(nightly, server="cgi")
else:
    bottle.run(nightly, host="localhost", port=8080, debug=True, 
               reloader=True)
