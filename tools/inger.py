#! /usr/bin/env python

#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

import sys
import getopt
import re
import string
import shutil

lineno=1
verbose=False

##########################

def main():

    global lineno, verbose
    
    # command line arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'f:c:hv', ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    function_filename=config_filename=""

    for o,a in optlist:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o == "-f":
            function_filename = a
        if o == "-c":
            config_filename = a
        if o == "-v":
            verbose = True

    if config_filename=="" or function_filename=="":
        print("Error, config file or function file not given")
        usage()
        sys.exit(2)

    # read the config file
    templates=dict()
    types=dict()
    lineno=1
    try:
        config_file = open(config_filename, "r")
        config_str = config_file.read()
        config_file.close()
    except IOError:
        print("Error reading config file "+config_filename)
        sys.exit(3)
    config_str=remove_comments(config_str)
    if verbose:
        print("Parsing config file: " + config_filename)
    parse_configfile(config_str, templates, types)

    # read the function file
    functions=dict()
    lineno=1
    try:
        function_file = open(function_filename, "r")
        function_str = function_file.read()
        function_file.close()
    except IOError:
        print("Error reading function file: "+function_filename)
        sys.exit(4)
    function_str=remove_comments(function_str)
    if verbose:
        print("Parsing function file " + function_filename)
    parse_functionfile(function_str, functions)

    # Ok, for each function and each template, we do the job
    for tempname in templates:
        if verbose:
            print("Applying template " + tempname)
        apply_template(templates[tempname], functions)

########################

def apply_template(temp, functions):

    # check the output mode, the default is "append"
    if not temp.has_key("OMODE"):
        temp["OMODE"] = "append"        
    if temp["OMODE"] == "append":
        if not temp.has_key("OUTPUT"):
            print("Error: output file not defined in template")
            sys.exit()
        shutil.copyfile(temp["OUTPUT"]+".in", temp["OUTPUT"])
        output=open(temp["OUTPUT"], "a")
    else:
        print("Unkwown output mode")
        sys.exit()

    for funcname in functions:
        entries=dict()

        if verbose:
            print("  Function " + funcname)

        ind=0
        text=temp["TEMPLATE"]
        while True:
            so=re.compile("\%([A-Z]+)\%").search(text, ind)
            if so==None: break
            ind=so.end(1)+2
            entry=text[so.start(1):so.end(1)]
            value=get_entry(entry, temp, functions[funcname])
            entries[entry]=value

        for entryname in entries:
            text=text.replace("%"+entryname+"%", entries[entryname])

        # what to do with this?
        if temp["OMODE"] == "append":
            output.write(text)
            output.write("\n")

    # clean up
    if temp["OMODE"] == "append":
        output.close()

def get_entry(entry, temp, func):
    res="%"+entry+"%"
    ############### CNAME #####################
    if entry=="CNAME":
        res=func["NAME"]
    ############### INAME #####################
    elif entry=="INAME":
        lnames=func["INAME"]
        for ln in lnames:
            if ln[0]==None or ln[0]==temp["NAME"]:
                res=ln[1]
    ############### IRETURN ###################
    elif entry=="IRETURN":
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        res=ipref + func["RETURN"][0]
    ############### DECL ######################
    elif entry=="DECL":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            typerec=temp["types"][arg[1]]
            if typerec.has_key("DECL"):
                decl=typerec["DECL"]
                decl=decl.replace("%C%", cpref+arg[2])
                res=res+decl
            else:
                decl=typerec["CTYPE"] + " " + cpref+arg[2] + ";"
                res=res+decl
            if arg[0]=="OUT":
                decl=typerec["ITYPE"] + " " + ipref+arg[2] + ";"
                res=res+decl
    ############### INCONV ####################
    elif entry=="INCONV":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            typerec=temp["types"][arg[1]]
            if typerec.has_key("INCONV") and (arg[0]=="IN" or arg[0]=="INOUT"):
                inconv=typerec["INCONV"]
                inconv=inconv.replace("%C%", cpref+arg[2])
                inconv=inconv.replace("%I%", ipref+arg[2])
                if res != "" and inconv != "":
                    res=res+"\n"
                res=res+inconv
    ############### OUTCONV ###################
    elif entry=="OUTCONV":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            typerec=temp["types"][arg[1]]
            if typerec.has_key("OUTCONV") and (arg[0]=="OUT" or arg[0]=="INOUT"):
                outconv=typerec["OUTCONV"]
                outconv=outconv.replace("%C%", cpref+arg[2])
                outconv=outconv.replace("%I%", ipref+arg[2])
                if res != "":
                    res=res+"\n"
                res=res+outconv
    ############### CPAR ######################
    elif entry=="CPAR":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        for arg in func["ARG"]:
            if res != "":
                res=res+", "
            res=res+cpref+arg[2]
    ############### FCPAR #####################
    elif entry=="FCPAR":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        for arg in func["ARG"]:
            typerec=temp["types"][arg[1]]
            name=cpref+arg[2]
            if typerec.has_key("FCALL"):
                name=typerec["FCALL"].replace("%C%", name)
            if res != "":
                res=res+", "
            res=res+name

    ############### CINPAR ####################
    elif entry=="CINPAR":
        res=""
        cpref="c_"
        if temp.has_key("CPREF"):
            cpref=temp["CPREF"]
        for arg in func["ARG"]:
            if arg[0]=="IN" or arg[0]=="INOUT":
                if res != "":
                    res=res+", "
                res=res+cpref+arg[2]
    ############### IINPAR ####################
    elif entry=="IINPAR":
        res=""
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            if arg[0]=="IN" or arg[0]=="INOUT":
                if res != "":
                    res=res+", "
                res=res+ipref+arg[2]
    ############### DIINPAR ###################
    elif entry=="DIINPAR":
        res=""
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            if arg[0]=="IN" or arg[0]=="INOUT":
                if res != "":
                    res=res+", "
                res=res+ipref+arg[2]
                if func.has_key("DEFAULT"):
                    for defa in func["DEFAULT"]:
                        if defa[0]==temp["NAME"] or defa[0]==None:
                            if defa[1]==arg[2]:
                                res=res+"="+defa[2]
    ############### TIINPAR ###################
    elif entry=="TIINPAR":
        res=""
        ipref=""
        if temp.has_key("IPREF"):
            ipref=temp["IPREF"]
        for arg in func["ARG"]:
            if arg[0]=="IN" or arg[0]=="INOUT":
                typerec=temp["types"][arg[1]]
                if res != "":
                    res=res+", "
                res=res+typerec["ITYPE"]+" " + ipref+arg[2]
    ############### IRETTYPE ##################    
    elif entry=="IRETTYPE":
        res=None
        retvar=func["RETURN"][0]
        for arg in func["ARG"]:
            if arg[2]==retvar:
                typerec=temp["types"][arg[1]]
                res=typerec["ITYPE"]
    ###########################################
    else:
        print "Error, unknown entry in template:", entry
        sys.exit()
    
    return res

########################

def remove_comments(string):
    res=""
    i=0
    while i<len(string):
        if string[i]=="#":
            while i<len(string) and not string[i]=="\n":
                i=i+1
        elif string[i]=="\\" and i<len(string)-1:
            i=i+1
            res=res+string[i]
        elif string[i]=="\\":
            i=i+1
        else:
            res=res+string[i]
            i=i+1
        
    return res

########################

def is_comment(str):
    return re.match("^\s*$", str) or re.match("^\s*#.*$", str)

def skip_lines(lines):
    global lineno
    while len(lines) > 0 and is_comment(lines[0]):
        lines=lines[1:]
        lineno=lineno+1
    return lines

def parse_configfile(config_str, templates, types):
    lines=config_str.splitlines()
    lines=skip_lines(lines)
    while len(lines) != 0:
        lines=parse_template(lines, templates)
        lines=skip_lines(lines)

def parse_template(lines, templates):
    global lineno, verbose
    commands=dict()
    types=dict()
    # header line
    mo=re.match("TEMPLATE[ ]([\w\-.]+)\s*:\s*$", lines[0])
    if not mo:
        print("Template file parse error in line "+str(lineno))
        print(lines[0])
        sys.exit()
    template_name=mo.group(1)
    commands["NAME"] = template_name
    lines=lines[1:]
    lineno=lineno+1
    if verbose:
        print("  Parsing template " + template_name)

    # commands
    while len(lines) > 0 and re.match("^\t.*$", lines[0]):
        mo=re.match("^\t\s*([\w\-\.]+)\s*:\s*(.*)$", lines[0])
        command=mo.group(1)
        arg=mo.group(2)
        if (commands.has_key(command)):
            print("Warning: duplicate template command: " + command +
                  " in line " + str(lineno))
        lines=lines[1:]
        lineno=lineno+1
        # multi-line arguments
        while len(lines) > 0 and re.match("\t\t.*$", lines[0]):
            arg=arg.strip(" \t\n") + "\n" + lines[0][2:].rstrip(" \t\n")
            lines=lines[1:]
            lineno=lineno+1
        # done
        commands[command] = arg
        if command == "TYPES" and arg != "":
            if verbose:
                print("    Parsing types from " + arg)
            type_file = open(arg, "r")
            type_str = type_file.read()
            type_file.close()
            type_str=remove_comments(type_str)
            parse_typefile(type_str, types)

    # return the remainder
    if len(types) > 0: commands["types"] = types
    templates[template_name]=commands
    return lines

def parse_typefile(type_str, types):
    global lineno
    oldlineno=lineno
    lineno=0
    lines = type_str.splitlines()
    lines = skip_lines(lines)
    while len(lines) != 0:
        lines=parse_type(lines, types)
        lines=skip_lines(lines)
    lineno=oldlineno

def parse_type(lines, types):
    global lineno, verbose
    entries=dict()
    # header line
    mo=re.match("^TYPE[ ]([\w_]+)\s*:\s*$", lines[0])
    if not mo:
        print("Type file parse error in line "+str(lineno))
        print(lines[0])
        sys.exit()
    type_name=mo.group(1)
    if verbose:
        print("      Parsing type " + type_name)
    lines=lines[1:]
    lineno=lineno+1

    # entries
    while len(lines) > 0 and len(lines[0]) > 0 and lines[0][0]=="\t":
        mo=re.match("^\t\s*([\w\-\.]+)\s*:\s*(.*)$", lines[0])
        entry=mo.group(1)
        arg=mo.group(2)
        lines=lines[1:]
        lineno=lineno+1
        # multi-line arguments
        while len(lines) >0 and re.match("\t\t.*$", lines[0]):
            arg=arg.strip(" \t\n") + "\n" + lines[0][2:].rstrip(" \t\n")
            lines=lines[1:]
            lineno=lineno+1
        # done
        entries[entry] = arg

    types[type_name] = entries
    return lines

########################

def parse_functionfile(function_str, functions):
    lines=function_str.splitlines()
    lines=skip_lines(lines)
    while len(lines) != 0:
        lines=parse_function(lines, functions)
        lines=skip_lines(lines)

def parse_function(lines, functions):
    global lineno, verbose
    entries=dict()
    # header line
    mo=re.match("^FUNCTION[ ]([\w_]+)\s*:\s*$", lines[0])
    if not mo:
        print("Function file parse error in line "+str(lineno))
        print(lines[0])
        sys.exit()
    function_name=mo.group(1)
    entries["NAME"]=function_name
    lines=lines[1:]
    lineno=lineno+1
    if (verbose):
        print("  Parsing function " + function_name)

    # entries
    while len(lines) > 0 and len(lines[0]) > 0 and lines[0][0]=="\t":
        mo=re.match("^\t\s*([\w\-\.]+)\s*:\s*(.*)$", lines[0])
        entry=mo.group(1)
        arg=mo.group(2).strip(" \t\n")
        lines=lines[1:]
        lineno=lineno+1
        # multi-line arguments
        while len(lines) > 0 and re.match("\t\t.*$", lines[0]):
            arg=arg.strip(" \t\n") + "\n" + lines[0][2:].rstrip(" \t\n")
            lines=lines[1:]
            lineno=lineno+1
        # done
        arg=parse_function_entry(entry, arg)
        if entries.has_key(entry):
                entries[entry].append(arg)
        else:
            entries[entry] = [ arg ]

    # return the remainder
    functions[function_name] = entries
    return lines

def parse_function_entry(entry, arg):
    global lineno
    if entry=="ARG":
        words=arg.split()
        if len(words)==3:
            res=( words[0], words[1], words[2] )
        elif len(words)==2:
            res=("IN", words[0], words[1])
        else:
            print("Error in function entry, line" + lineno-1)
            sys.exit()
    elif entry=="RETURN":
        res=arg
    elif entry=="DEFAULT":
        mo=re.match("^([^\"]*)\s*\"(.*)\"$", arg)
        if not mo:
            print("Error in function entry, line" + lineno-1)
            sys.exit()
        value=mo.group(2)
        words=mo.group(1).split()
        if len(words)==1:
            res=( None, words[0], value)
        elif len(words)==2:
            res=( words[0], words[1], value)
        else:
            print("Error in function entry, line" + lineno-1)
            sys.exit()
    elif entry=="INAME":
        words=arg.split()
        if len(words)==2:
            res=(words[0], words[1])
        else:
            res=(None, words[0])
    else:
        print "Error: unknown function entry in line", lineno-1
        sys.exit()

    return res

def usage():
    print("Usage: " + sys.argv[0] + " [-hv] [--help] " +
          " -c config-file -f function-file")

########################

if __name__ == "__main__":
    main()

