#! /bin/env /home/csardi/software/bin/mypython

import sqlite3
import sys

if len(sys.argv) <= 1:
    print 'Give database file as argument'
    sys.exit(1)

con=sqlite3.connect(sys.argv[1])
cur=con.cursor()

con.execute('ALTER TABLE network ADD COLUMN directed BOOLEAN;')
con.execute('ALTER TABLE network ADD COLUMN bipartite BOOLEAN;')
con.execute('ALTER TABLE network ADD COLUMN weighted BOOLEAN;')
con.execute('ALTER TABLE network ADD COLUMN sid TEXT;')

con.execute('UPDATE network SET directed=0;')
con.execute('UPDATE network SET bipartite=0;')
con.execute('UPDATE network SET weighted=0;')
