#! /usr/bin/env python

import re
import seqdict
import sys
import getopt
import os

version="0.1"
date="Jul 29 2007"

def usage():
    print "Stimulus version", version, date
    print sys.argv[0], "-f <function-file> -t <type-file> -l language "
    print ' ' * len(sys.argv[0]), "-i <input-file> -o <output-file>"
    print ' ' * len(sys.argv[0]), "-h --help -v"

################################################################################
class PLexer:
    def __init__(self, stream):
        self.stream=stream
        self.ws_stack=[0]
        self.tokens=[]
        self.lineno=0

    def lineno(self):
        return self.lineno

    def token(self):
        keys=[]

        if (len(self.tokens)>0):
            return self.tokens.pop(0)

        # Read a line, skip empty lines and comments
        while True:
            line=self.stream.readline(); self.lineno = self.lineno+1
            if line=="":
                for k in keys:
                    self.tokens.append( ("key", k) )
                    keys=[]
                while len(self.ws_stack)>0:
                    self.tokens.append( ("dedent", "") )
                    self.ws_stack.pop()
                self.tokens.append( ("eof", "") )
                return self.tokens.pop(0)
            if re.match("^[ \t]*$", line): continue
            if re.match("^[ \t]*#", line): continue
            break

        if line[-1]=="\n": line=line[:(len(line)-1)]
        ws=re.match(r"^[ \t]*", line).span()[1]
        line=line.strip()
        if ws > self.ws_stack[-1]:
            self.tokens.append( ("indent", "") )
            self.ws_stack.append(ws)
        else:
            for k in keys:
                self.tokens.append( ("key", k) )
                keys=[]
            while ws < self.ws_stack[-1]:
                self.ws_stack.pop()
                self.tokens.append( ("dedent", "") )
            if ws != self.ws_stack[-1]:
                print "Bad indentation in line", self.lineno
                exit

        # Ok, we're done with the white space, now let's see
        # whether this line is continued
        while line[-1]=="\\":
            line=line[:(len(line)-1)]
            line=line+"\n  " + self.stream.readline().strip() ; self.lineno=self.lineno+1

        # We have the line now, check whether there is a ':' in it
        line=line.split(":", 1)
        if len(line)>1:
            line[0]=line[0].strip()
            line[1]=line[1].strip()
            if line[0]=="":
                print "Missing keyword in line", self.lineno
                exit
            keys=line[0].split(",")
            keys=[ k.strip() for k in keys ]
            if line[1] == "":
                self.tokens.append( ("key", keys.pop(0)) )
            else:
                for k in keys:
                    self.tokens.append( ("key", k))
                    self.tokens.append( ("indent", "") )
                    self.tokens.append( ("text", line[1]) )
                    self.tokens.append( ("dedent", "") )
        else:
            self.tokens.append( ("text", line[0].strip()) )
            for k in keys:
                self.tokens.append( ("dedent", "") )
                self.tokens.append( ("key", k) )
                self.tokens.append( ("indent", "") )
            keys=[]
                
        if self.tokens:
            return self.tokens.pop(0)

################################################################################
class PParser:
    def parse(self, stream):
        lex=PLexer(stream)
        val=seqdict.seqdict()
        val_stack=[val, None]
        nam_stack=[None, None]
        
        tok=lex.token()
        while not tok[0]=="eof":
            if tok[0]=="indent":
                val_stack.append(None)
                nam_stack.append(None)
            elif tok[0]=="dedent":
                v=val_stack.pop()
                n=nam_stack.pop()
                if n is None:                    
                    val_stack[-1]=v
                else:
                    val_stack[-1][n]=v
            elif tok[0]=="key":
                if not nam_stack[-1] is None:
                    val_stack[-2][nam_stack[-1]]=val_stack[-1]
                if tok[1][-5:]=="-list":
                    val_stack[-1]=seqdict.seqdict()
                    nam_stack[-1]=tok[1][:-5]
                else:
                    val_stack[-1]={}
                    nam_stack[-1]=tok[1]
            elif tok[0]=="text":
                val_stack[-1]=tok[1]
            tok=lex.token()

        return val

################################################################################

def main():
    # Command line arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 't:f:l:i:o:hv', ['help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    types=[]; functions=[]; inputs=[]; languages=[]; outputs=[]; verbose=False

    for o,a in optlist:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o == "-o":
            outputs.append(a)
        elif o == "-t":
            types.append(a)
        elif o == "-f":
            functions.append(a)        
        elif o == "-l":
            languages.append(a)
        elif o == "-i":
            inputs.append(a)
        elif o =="-v":
            verbose=True

    # Parameter checks
    # Note: the lists might be empty, but languages and outputs must
    # have the same length.
    if len(languages) != len(outputs):
        print "Error: number of languages and output files must match"
        sys.exit(4)
    for l in languages:
        if not l+"CodeGenerator" in globals():
            print "Error: unknown language:", l
            sys.exit(6)
    for f in types:
        if not os.access(f, os.R_OK):
            print "Error: cannot open type file:", f
            sys.exit(5)
    for f in functions:
        if not os.access(f, os.R_OK):
            print "Error: cannot open function file:", f
            sys.exit(5)
    for f in inputs:
        if not os.access(f, os.R_OK):
            print "Error: cannot open input file:", f
            sys.exit(5)
    # TODO: output files are not checked now

    # OK, do the trick:
    for l in range(len(languages)):
        cl=globals()[languages[l]+"CodeGenerator"]
        cg=cl(functions, types)
        cg.generate(inputs, outputs[l])

################################################################################
class CodeGenerator:
    def __init__(self, func, types):
        # Set name
        self.name=str(self.__class__).split(".")[-1]
        self.name=self.name[0:len(self.name)-len("CodeGenerator")]

        # Parse function and type files
        parser=PParser()
        self.func=seqdict.seqdict()
        for f in func:
            ff=open(f)
            newfunc=parser.parse(ff)
            self.func.extend(newfunc)
            ff.close()

        self.types=seqdict.seqdict()
        for t in types:
            ff=open(t)
            newtypes=parser.parse(ff)
            self.types.extend(newtypes)
            ff.close()

    def generate(self, inputs, output):
        out=open(output, "w")
        self.append_inputs(inputs, out)
        for f in self.func.keys():
            if 'FLAGS' in self.func[f]:
                flags=self.func[f]['FLAGS']
                flags=flags.split(",")
                flags=[ flag.strip() for flag in flags ]
            else:
                self.func[f]['FLAGS']=[]
            self.generate_function(f, out)
        out.close()

    def generate_function(self, f, out):
        print "Error: invalid code generator, this method should be overridden"
        sys.exit(1)

    def append_inputs(self, inputs, output):
        for i in inputs:
            ii=open(i)
            str=ii.read()
            while str != "":
                output.write(str)
                str=ii.read()
            ii.close()
        pass    

    def ignore(self, function):
        if 'IGNORE' in self.func[function]:
            ign=self.func[function]['IGNORE']
            ign=ign.split(",")
            ign=[i.strip() for i in ign]
            if self.name in ign: return True
        return False

################################################################################
# GNU R, see http://www.r-project.org
################################################################################

class RNamespaceCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate(self, inputs, output):
        """This is very simple, we include an 'export' line for every
        function which it not to be ignored by the RNamespace language.
        Function names are taken from NAME-R if present, otherwise
        underscores are converted to dots and the leading 'i' (from
        'igraph') is stripped to create the function name,
        ie. igraph_clusters is mapped to graph.clusters."""
        out=open(output, "w")
        self.append_inputs(inputs, out)
        for f in self.func.keys():
            if (self.ignore(f)):
                continue
            name=self.func[f].get("NAME-R", f[1:].replace("_", "."))
            out.write("export(" + name + ")\n")
        out.close()
        

class RRCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate_function(self, function, out):
        
        # Ignore?
        if self.ignore(function):
            return

        name=self.func[function].get("NAME-R", function[1:].replace("_", "."))
        if "PARAM" in self.func[function]:            
            params=self.func[function]["PARAM"]
        else:
            params=seqdict.seqdict()

        # Check types
        for p in params.keys():
            tname=params[p]['type']
            if not tname in self.types.keys():
                print "Error: Unknown type encountered:", tname
                sys.exit(7)
            params[p].setdefault('mode', 'IN')

        ## Header
        ## do_par handles the translation of a single argument in the
        ## header. Pretty simple, the only difficulty is that we
        ## might need to add default values. Default values are taken
        ## from a language specific dictionary, this is compiled from
        ## the type file(s).
        
        ## So we take all arguments with mode 'IN' or 'INOUT' and
        ## check whether they have a default value. If yes then we 
        ## check if the default value is given in the type file. If
        ## yes then we use the value given there, otherwise the
        ## default value is ignored silently. (Not very nice.)
        
        ## We also check whether the function has a 'PROGRESS' flag.
        ## If yes then there is an extra 'verbose'
        ## argument. Alternatively we could just always add this
        ## argument independently of the 'PROGRESS' flag.
        out.write(name)
        out.write(" <- function(")
        def do_par(pname):
            tname=params[pname]['type']
            t=self.types[tname]
            default=""
            if 'default' in params[pname]:
                if 'DEFAULT' in t and params[pname]['default'] in t['DEFAULT']:
                    default="=" + t['DEFAULT'][ params[pname]['default'] ]
            return pname.replace("_", ".") + default
            
        head=[ do_par(n) for n,p in params.items()
               if p['mode'] in ['IN','INOUT'] ]
        if 'PROGRESS' in self.func[function]['FLAGS']:
            head.append('verbose=igraph.par("verbose")')
        out.write(", ".join(head))
        out.write(") {\n")

        ## Argument checks, INCONV
        ## We take 'IN' and 'INOUT' mode arguments and if they have an
        ## INCONV field then we use that. This is typically for
        ## argument checks, like we check here that the argument
        ## supplied for a graph is indeed an igraph graph object. We
        ## also covert numeric vectors to 'double' here.

        ## The INCONV fields are simply concatenated by newline
        ## characters.
        out.write("  # Argument checks\n")
        def do_par(pname):
            t=self.types[params[pname]['type']]
            m=params[pname]['mode']
            if m in ['IN', 'INOUT'] and 'INCONV' in t:
                res="  " + t['INCONV']
            else:
                res=""
            return res.replace("%I%", pname.replace("_", "."))

        inconv=[ do_par(n) for n in params.keys() ]
        inconv=[ i for i in inconv if i != "" ]
        out.write("\n".join(inconv)+"\n\n")

        ## Function call
        ## This is a bit more difficult than INCONV. Here we supply
        ## each argument to the .Call function, if the argument has a
        ## 'CALL' field then it is used, otherwise we simply use its
        ## name. We also need to take care about the PROGRESS bar
        ## argument. Note that arguments with empty CALL fields are
        ## completely ignored, so giving an empty CALL field is
        ## different than not giving it at all.

        ## The tail of the function is also written and we're ready.
        def do_par(pname):
            t=self.types[params[pname]['type']]
            call=pname.replace("_", ".")
            if 'CALL' in t:
                call=t['CALL']
                if call:
                    call=call.replace('%I%', pname.replace("_", "."))
                else:
                    call=""
            return call
              
        out.write("  # Function call\n")
        out.write("  .Call(\"R_" + function + "\", ")
        call=[ do_par(n) for n,p in params.items() if p['mode'] in ['IN', 'INOUT'] ]
        call=[ c for c in call if c != "" ]
        if 'PROGRESS' in self.func[function]['FLAGS']:
            call.append('as.logical(verbose)')
        out.write(", ".join(call))
        out.write(",\n        PACKAGE=\"igraph\")\n")

        out.write("}\n\n")

class RCCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate_function(self, function, out):

        # Ignore?
        if self.ignore(function):
            return

        if "PARAM" in self.func[function]:            
            params=self.func[function]["PARAM"]
        else:
            params=seqdict.seqdict()

        # Check types
        for p in params.keys():
            tname=params[p]['type']
            if not tname in self.types.keys():
                print "Error: Unknown type encountered:", tname
                sys.exit(7)
            params[p].setdefault('mode', 'IN')

        ## Compile the output
        ## This code generator is quite difficult, so we use different
        ## functions to generate the approprite chunks and then
        ## compile them together using a simple template.
        ## See the documentation of each chunk below.
        res={}
        res['func']=function
        res['header']=self.chunk_header(function, params)
        res['decl']=self.chunk_declaration(function, params)
        res['before']=self.chunk_before(function, params)
        res['inconv']=self.chunk_inconv(function, params)
        res['call']=self.chunk_call(function, params)
        res['outconv']=self.chunk_outconv(function, params)
        res['after']=self.chunk_after(function, params)

        # Replace into the template
        text="""
/*-------------------------------------------/
/ %(func)-42s /
/-------------------------------------------*/
%(header)s {
                                        /* Declarations */
%(decl)s

%(before)s
                                        /* Convert input */
%(inconv)s
                                        /* Call igraph */
%(call)s
                                        /* Convert output */
%(outconv)s

%(after)s

  UNPROTECT(1);
  return(result);
}\n""" % res

        out.write(text)
        
    def chunk_header(self, function, params):
        """The header. All functions return with a 'SEXP', so this is
        easy. We just take the 'IN' and 'INOUT' arguments, all will
        have type SEXP, and concatenate them by commas. If the
        function has a 'PROGRESS' flag that is used. The function name
        is created by prefixing the original name with 'R_'."""
        inout=[ "SEXP "+n for n,p in params.items() if p['mode'] in ['IN','INOUT'] ]
        if 'PROGRESS' in self.func[function]['FLAGS']:
            inout.append("SEXP verbose")
        return "SEXP R_" + function + "(" + ", ".join(inout) + ")"
        
    def chunk_declaration(self, function, params):
        """There are a couple of things to declare. First a C type is
        needed for every argument, these will be supplied in the C
        igraph call. Then, all 'OUT' arguments need a SEXP variable as
        well, the result will be stored here. The return type
        of the C function also needs to be declared, that comes
        next. The result and names SEXP variables will contain the
        final result, these are last. ('names' is not always used, but
        it is easier to always declare it.)
        """
        def do_par(pname):
            cname="c_"+pname
            t=self.types[params[pname]['type']]
            if 'DECL' in t:
                decl="  " + t['DECL']
            elif 'CTYPE' in t:
                decl="  " + t['CTYPE'] + " " + cname + ";"
            else:
                decl=""            
            return decl.replace("%C%", cname).replace("%I%", pname)

        inout=[ do_par(n) for n in params.keys() ]
        out=[ "  SEXP "+n+";" for n,p in params.items()
              if p['mode']=='OUT' ]

        rt=self.types[self.func[function]['RETURN']]    
        if 'DECL' in rt:
            retdecl="  " + rt['DECL']
        elif 'CTYPE' in rt:
            retdecl="  " + rt['CTYPE'] + " c_result;"
        else:
            retdecl=""

        return "\n".join(inout + out + [retdecl] + ["  SEXP result, names;"])

    def chunk_before(self, function, params):
        """Quite confusingly, there is a separate igraph_before()
        function for functions with progress bar, so we call that if
        the function has a 'PROGRESS' flag. This should be eliminated
        in the future."""
        if 'PROGRESS' in self.func[function]['FLAGS']:
            return '  R_igraph_before2(verbose, "");'
        else:
            return '  R_igraph_before();'        

    def chunk_inconv(self, function, params):
        """Input conversions. Not only for types with mode 'IN' and
        'INOUT', eg. for 'OUT' vector types we need to allocate the
        required memory here, do all the initializations, etc. Types
        without INCONV fields are ignored. The usual %C%, %I% is
        performed at the end.
        """
        def do_par(pname):
            cname="c_"+pname
            t=self.types[params[pname]['type']]
            mode=params[pname]['mode']
            if 'INCONV' in t and mode in t['INCONV']:
                inconv="  " + t['INCONV'][mode]
            else:
                inconv=""

            return inconv.replace("%C%", cname).replace("%I%", pname)

        inconv=[ do_par(n) for n in params.keys() ]
        inconv=[ i for i in inconv if i != "" ]

        return "\n".join(inconv)

    def chunk_call(self, function, params):
        """Every single argument is included, independently of their
        mode. If a type has a 'CALL' field then that is used after the
        usual %C% and %I% substitutions, otherwise the standard 'c_'
        prefixed C argument name is used.
        """
        types=[ self.types[params[n]['type']] for n in params.keys() ]
        call=map( lambda t, n: t.get('CALL', "c_"+n), types, params.keys() )
        call=map( lambda c, n: c.replace("%C%", "c_"+n).replace("%I%", n),
                  call, params.keys() )
        return "  c_result=" + function + "(" + ", ".join(call) + ");\n"

    def chunk_outconv(self, function, params):
        """The output conversions, this is quite difficult. A function
        may report its results in two ways: by returning it directly
        or by setting a variable to which a pointer was passed. igraph
        usually uses the latter and returns error codes, except for
        some simple functions like 'igraph_vcount()' which cannot
        fail.

        First we add the output conversion for all types. This is
        easy. Note that even 'IN' arguments may have output
        conversion, eg. this is the place to free memory allocated to
        them in the 'INCONV' part.

        Then we check how many 'OUT' or 'INOUT' arguments we
        have. There are three cases. If there is a single such
        argument then that is already converted and we need to return
        that. If there is no such argument then the output of the
        function was returned, so we perform the output conversion for
        the returned type and this will be the result. If there are
        more than one 'OUT' and 'INOUT' arguments then they are
        collected in a named list. The names come from the argument
        names.
        """
        def do_par(pname):
            cname="c_"+pname
            t=self.types[params[pname]['type']]
            mode=params[pname]['mode']
            if 'OUTCONV' in t and mode in t['OUTCONV']:
                outconv="  " + t['OUTCONV'][mode]
            else:
                outconv=""
            return outconv.replace("%C%", cname).replace("%I%", pname)

        outconv=[ do_par(n) for n in params.keys() ]
        outconv=[ o for o in outconv if o != "" ]

        retpars=[ n for n,p in params.items() if p['mode'] in ['OUT', 'INOUT'] ]
        if len(retpars)==0:
            # return the return value of the function
            rt=self.types[self.func[function]['RETURN']]
            if 'OUTCONV' in rt:
                retconv="  " + rt['OUTCONV']['OUT']
            else:
                retconv=""
            retconv=retconv.replace("%C%", "c_result").replace("%I%", "result")
            ret="\n".join(outconv) + "\n" + retconv
        elif len(retpars)==1:
            # return the single output value
            retconv="  result=" + retpars[0] + ";"
            ret="\n".join(outconv) + "\n" + retconv
        else:
            # create a list of output values
            sets=map ( lambda c, n: "  SET_VECTOR_ELT(result, "+str(c)+", "+n+");",
                       range(len(retpars)), retpars )
            names=map ( lambda c, n: "  SET_STRING_ELT(names, "+str(c)+
                        ", CREATE_STRING_VECTOR(\""+n+"\"));",
                        range(len(retpars)), retpars )
            ret="\n".join(["  PROTECT(result=NEW_LIST(" + str(len(retpars)) + "));",
                           "  PROTECT(names=NEW_CHARACTER(" + str(len(retpars)) + "));"]+
                          outconv + sets + names +
                          ["  SET_NAMES(result, names);" ] + 
                          ["  UNPROTECT("+str(len(sets)+1)+");" ])
            
        return ret

    def chunk_after(self, function, params):
        """The pair of igraph_before(), different for functions with
        progreess reporting."""
        if 'PROGRESS' in self.func[function]['FLAGS']:
            return '  R_igraph_after2(verbose);'
        else:
            return '  R_igraph_after();'

################################################################################
if __name__ == "__main__":
    main()

