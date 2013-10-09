#! /usr/bin/env python

import cgi
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
if myname[-6] == ".local":
    bottle.run(nightly, host="localhost", port=8080, debug=True, 
               reloader=True)
else:
    bottle.run(nightly, server="cgi")
