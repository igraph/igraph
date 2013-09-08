#! /usr/bin/env/python

#
#
# TO create the DB, run:
"""
echo "
DROP TABLE 'downloads';
CREATE TABLE 'downloads' (
       type     TEXT,
       version  TEXT DEFAULT "",
       branch   TEXT,
       hash     TEXT,
       date     TEXT,
       size     INTEGER,
       count    INTEGER DEFAULT 0,
       PRIMARY KEY(type,hash)
);
" | sqlite3 nightly.db 
"""

# Some example data:
"""
echo "
INSERT INTO 'downloads'
       VALUES ('C library', '0.7', 'develop', 'f5115b6', '2013-09-04', 
               '1200000', 0);
INSERT INTO 'downloads'
       VALUES ('R package', '0.6.5', 'master', 'f5115ba', '2013-09-03', 
               '1200000', 0);
INSERT INTO 'downloads'
       VALUES ('Python extension', '0.6.5', 'master', 'f5115bb', '2013-08-05', 
               '1200000', 0);
INSERT INTO 'downloads'
       VALUES ('C library', '0.7', 'develop', 'f5115bc', '2013-09-02', 
               '1200000', 0);
INSERT INTO 'downloads'
       VALUES ('R package', '0.7', 'develop', 'f5115bd', '2013-09-01', 
               '1200000', 0);
" | sqlite3 nightly.db

"""

import bottle
import bottle_sqlite
import socket

plugin = bottle_sqlite.Plugin(dbfile='nightly.db')
nightly=bottle.Bottle()
nightly.install(plugin)

urlmap={ 'C library': 'c', 'R package': 'r', 'Python extension': 'python',
         'C library for MSVC': 'msvc' }
revurlmap=dict((v,k) for k, v in urlmap.iteritems())

# This is the main web page, you can choose your download here
# It can be filtered
@nightly.route("/")
@nightly.route("/list/")
@nightly.route("/list")
@nightly.route("/list/<dtype>")
@nightly.route("/list/<dtype>/<version>")
@nightly.route("/list/<dtype>/<version>/<branch>")
def list_files(db, dtype="all", version="all", branch="all"):

    if dtype=="all":
        dtype="%"
    else:
        dtype=revurlmap[dtype]
    if version=="all":
        version="%"
    if branch=="all":
        branch="%"

    files = db.execute("SELECT type, version, branch, hash, \
                               date, size FROM downloads    \
                        WHERE type LIKE ? AND version LIKE ? AND branch LIKE ? \
                        ORDER BY version DESC, date DESC",
                       (dtype, version, branch)).fetchall()

    versions=db.execute("SELECT DISTINCT version FROM downloads \
                         ORDER BY version DESC").fetchall()
    versions=[ e[0] for e in versions ]

    types=db.execute("SELECT DISTINCT type FROM downloads \
                      ORDER BY type").fetchall()
    types=[ e[0] for e in types ]

    branches=db.execute("SELECT DISTINCT branch FROM downloads \
                         ORDER BY branch").fetchall()
    branches=[ e[0] for e in branches ]

    return bottle.template('main', files=files, versions=versions,
                           types=types, branches=branches, urlmap=urlmap)

# This is the one that serves the files. 
# <type> can be 'c', 'python', 'r' or 'msvc' and 
# <hash> is a git hash for the commit to get. If <hash> is 
# None, then the latest version is served
@nightly.route("/get/<dtype>/<hash>")
@nightly.route("/get/<dtype>")
def get_file(db, dtype, hash=None):
    ltype=revurlmap[dtype]
    if hash is None:
        hash=db.execute("SELECT hash, max(date) FROM downloads WHERE type=?",
                        (ltype,)).fetchone()[0]
    filename=dtype + "/" + hash
    db.execute("UPDATE downloads SET count=count+1 WHERE type=? AND hash=?", \
               (ltype, hash))
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
