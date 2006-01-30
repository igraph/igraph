#! /usr/bin/env python

cmdfile=open("atlas.py", "r")
cmd=cmdfile.read()
cmdfile.close()

exec(cmd.strip())

for g in descr_list:
    f=open("atlas-"+g[1]+".txt", "w")
    f.write(str(g[2]))
    f.write("\n")
    for edge in g[3]:
        f.write(str(edge[0]-1))
        f.write(" ")
        f.write(str(edge[1]-1))
        f.write("\n")
