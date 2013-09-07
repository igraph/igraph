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
@nightly.route("/")
def list_files(db):
    files = db.execute("SELECT type, version, branch, hash, \
                               date, size FROM downloads    \
                        ORDER BY version DESC, date DESC").fetchall()
    versions = sorted(list(set([ f['version'] for f in files ])), reverse=True)
    types = sorted(list(set([ f['type'] for f in files ])))
    branches = sorted(list(set([ f['branch'] for f in files ])))
    return bottle.template('main', files=files, versions=versions,
                           types=types, branches=branches, urlmap=urlmap)

# This is the one that serves the files. 
# <type> can be 'c', 'python', 'r' or 'msvc' and 
# <hash> is a git hash for the commit to get. If <hash> is 
# None, then the latest version is served
@nightly.route("/get/<type>/<hash>")
def get_file(db, type, hash=None):
    filename=type + "/" + hash
    db.execute("UPDATE downloads SET count=count+1 WHERE type=? AND hash=?", \
               (revurlmap[type], hash))
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
