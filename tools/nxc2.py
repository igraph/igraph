#! /usr/bin/env python

cmdfile=open("atlas.py", "r")
cmd=cmdfile.read()
cmdfile.close()

exec(cmd.strip())

p=0
plist=[]
f=open("atlas-edges.h", "w")
f.write("real_t igraph_i_atlas_edges[]={\n")
for g in descr_list:
    f.write("\t"+str(g[2])+","+str(len(g[3])))
    plist.append(p)
    p+=2
    for edge in g[3]:
        f.write(",")
        f.write(str(edge[0]-1))
        f.write(",")
        f.write(str(edge[1]-1))
        p+=2
    f.write(",\n")

f.write("};\n\nlong int igraph_i_atlas_edges_pos[]={")
f.write(", ".join(map(str,plist)))
f.write("};\n\n")

