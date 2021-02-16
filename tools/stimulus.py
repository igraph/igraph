#! /usr/bin/env python3

import re
import sys
from collections import OrderedDict
import getopt
import os
import abc

version = "0.2"
date = "Jan 13 2021"


def usage():
    print("Stimulus version", version, date)
    print(sys.argv[0], "-f <function-file> -t <type-file> -l language ")
    print(" " * len(sys.argv[0]), "-i <input-file> -o <output-file>")
    print(" " * len(sys.argv[0]), "-h --help -v")


################################################################################
class StimulusError(Exception):
    def __init__(self, message):
        self.msg = message

    def __str__(self):
        return str(self.msg)


################################################################################
class PLexer:
    def __init__(self, stream):
        self.stream = stream
        self.ws_stack = [0]
        self.tokens = []
        self.lineno = 0

    def lineno(self):
        return self.lineno

    def token(self):
        keys = []

        if len(self.tokens) > 0:
            return self.tokens.pop(0)

        # Read a line, skip empty lines and comments
        while True:
            line = self.stream.readline()
            self.lineno = self.lineno + 1
            if line == "":
                for k in keys:
                    self.tokens.append(("key", k))
                    keys = []
                while len(self.ws_stack) > 0:
                    self.tokens.append(("dedent", ""))
                    self.ws_stack.pop()
                self.tokens.append(("eof", ""))
                return self.tokens.pop(0)
            if re.match("^[ \t]*$", line):
                continue
            if re.match("^[ \t]*#", line):
                continue
            break

        if line[-1] == "\n":
            line = line[: (len(line) - 1)]
        ws = re.match(r"^[ \t]*", line).span()[1]
        line = line.strip()
        if ws > self.ws_stack[-1]:
            self.tokens.append(("indent", ""))
            self.ws_stack.append(ws)
        else:
            for k in keys:
                self.tokens.append(("key", k))
                keys = []
            while ws < self.ws_stack[-1]:
                self.ws_stack.pop()
                self.tokens.append(("dedent", ""))
            if ws != self.ws_stack[-1]:
                print("Bad indentation in line", self.lineno)
                exit

        # Ok, we're done with the white space, now let's see
        # whether this line is continued
        while line[-1] == "\\":
            line = line[: (len(line) - 1)]
            line = line + "\n  " + self.stream.readline().strip()
            self.lineno = self.lineno + 1

        # We have the line now, check whether there is a ':' in it
        line = line.split(":", 1)
        if len(line) > 1:
            line[0] = line[0].strip()
            line[1] = line[1].strip()
            if line[0] == "":
                print("Missing keyword in line", self.lineno)
                exit
            keys = line[0].split(",")
            keys = [k.strip() for k in keys]
            if line[1] == "":
                self.tokens.append(("key", keys.pop(0)))
            else:
                for k in keys:
                    self.tokens.append(("key", k))
                    self.tokens.append(("indent", ""))
                    self.tokens.append(("text", line[1]))
                    self.tokens.append(("dedent", ""))
        else:
            self.tokens.append(("text", line[0].strip()))
            for k in keys:
                self.tokens.append(("dedent", ""))
                self.tokens.append(("key", k))
                self.tokens.append(("indent", ""))
            keys = []

        if self.tokens:
            return self.tokens.pop(0)


################################################################################
class PParser:
    def parse(self, stream):
        lex = PLexer(stream)
        val = OrderedDict()
        val_stack = [val, None]
        nam_stack = [None, None]

        tok = lex.token()
        while not tok[0] == "eof":
            if tok[0] == "indent":
                val_stack.append(None)
                nam_stack.append(None)
            elif tok[0] == "dedent":
                v = val_stack.pop()
                n = nam_stack.pop()
                if n is None:
                    val_stack[-1] = v
                else:
                    val_stack[-1][n] = v
            elif tok[0] == "key":
                if not nam_stack[-1] is None:
                    val_stack[-2][nam_stack[-1]] = val_stack[-1]
                if tok[1][-5:] == "-list":
                    val_stack[-1] = OrderedDict()
                    nam_stack[-1] = tok[1][:-5]
                else:
                    val_stack[-1] = {}
                    nam_stack[-1] = tok[1]
            elif tok[0] == "text":
                val_stack[-1] = tok[1]
            tok = lex.token()

        return val


################################################################################


def main():
    # Command line arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "t:f:l:i:o:hv", ["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    types = []
    functions = []
    inputs = []
    languages = []
    outputs = []
    verbose = False

    for o, a in optlist:
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
        elif o == "-v":
            verbose = True

    # Parameter checks
    # Note: the lists might be empty, but languages and outputs must
    # have the same length.
    if len(languages) != len(outputs):
        print("Error: number of languages and output files must match")
        sys.exit(4)
    for l in languages:
        if not l + "CodeGenerator" in globals():
            print("Error: unknown language:", l)
            sys.exit(6)
    for f in types:
        if not os.access(f, os.R_OK):
            print("Error: cannot open type file:", f)
            sys.exit(5)
    for f in functions:
        if not os.access(f, os.R_OK):
            print("Error: cannot open function file:", f)
            sys.exit(5)
    for f in inputs:
        if not os.access(f, os.R_OK):
            print("Error: cannot open input file:", f)
            sys.exit(5)
    # TODO: output files are not checked now

    # OK, do the trick:
    for language, output in zip(languages, outputs):
        cl = globals()[language + "CodeGenerator"]
        cg = cl(functions, types)
        cg.generate(inputs, output)


################################################################################
class CodeGenerator(metaclass=abc.ABCMeta):
    def __init__(self, func, types):
        # Set name, note this only works correctly if derived classes always
        # extend it as by prepending the language to the CodeGenerator class
        # name
        self.name = type(self).__name__
        self.name = self.name[0 : len(self.name) - len("CodeGenerator")]

        # Parse function and type files
        parser = PParser()
        self.func = OrderedDict()
        for f in func:
            ff=open(f, "rU")
            newfunc=parser.parse(ff)
            self.func.extend(newfunc)
            ff.close()

        self.types = OrderedDict()
        for t in types:
            ff=open(t, "rU")
            newtypes=parser.parse(ff)
            self.types.extend(newtypes)
            ff.close()

        # The default return type is 'ERROR'
        for f in self.func.keys():
            if "RETURN" not in self.func[f]:
                self.func[f]["RETURN"] = "ERROR"

    def generate(self, inputs, output):
        out = open(output, "w")
        self.append_inputs(inputs, out)
        for f in self.func.keys():
            if "FLAGS" in self.func[f]:
                flags = self.func[f]["FLAGS"]
                flags = flags.split(",")
                flags = [flag.strip() for flag in flags]
            else:
                self.func[f]["FLAGS"] = []
            self.generate_function(f, out)
        out.close()

    @abc.abstractmethod
    def generate_function(self, f, out):
        raise NotImplementedError(
            "Error: invalid code generator, this method should be overridden"
        )

    def parse_params(self, function):
        if "PARAMS" not in self.func[function]:
            return OrderedDict()

        params = self.func[function]["PARAMS"]
        params = params.split(",")
        params = [p.strip() for p in params]
        params = [p.split(" ", 1) for p in params]
        for p in range(len(params)):
            if params[p][0] in ["OUT", "IN", "INOUT"]:
                params[p] = [params[p][0]] + params[p][1].split(" ", 1)
            else:
                params[p] = ["IN", params[p][0]] + params[p][1].split(" ", 1)
            if "=" in params[p][2]:
                params[p] = params[p][:2] + params[p][2].split("=", 1)
        params = [[p.strip() for p in pp] for pp in params]
        res = OrderedDict()
        for p in params:
            if len(p) == 3:
                res[p[2]] = {"mode": p[0], "type": p[1]}
            else:
                res[p[2]] = {"mode": p[0], "type": p[1], "default": p[3]}
        return res

    def parse_deps(self, function):
        if "DEPS" not in self.func[function]:
            return OrderedDict()

        deps = self.func[function]["DEPS"]
        deps = deps.split(",")
        deps = [d.strip() for d in deps]
        deps = [d.split("ON", 1) for d in deps]
        deps = [[dd.strip() for dd in d] for d in deps]
        deps = [[d[0]] + d[1].split(" ", 1) for d in deps]
        deps = [[dd.strip() for dd in d] for d in deps]
        res = OrderedDict()
        for d in deps:
            res[d[0]] = d[1:]
        return res

    def append_inputs(self, inputs, output):
        for i in inputs:
            ii=open(i, "rU")
            str=ii.read()
            while str != "":
                output.write(str)
                str = ii.read()
            ii.close()
        pass

    def ignore(self, function):
        if "IGNORE" in self.func[function]:
            ign = self.func[function]["IGNORE"]
            ign = ign.split(",")
            ign = [i.strip() for i in ign]
            if self.name in ign:
                return True
        return False


################################################################################
# GNU R, see http://www.r-project.org
# TODO: free memory when CTRL+C pressed, even on windows
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
        out = open(output, "w")
        self.append_inputs(inputs, out)
        for f in self.func.keys():
            if self.ignore(f):
                continue
            name = self.func[f].get("NAME-R", f[1:].replace("_", "."))
            out.write("export(" + name + ")\n")
        out.close()


class RRCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate_function(self, function, out):

        # Ignore?
        if self.ignore(function):
            return

        name = self.func[function].get("NAME-R", function[1:].replace("_", "."))
        params = self.parse_params(function)
        self.deps = self.parse_deps(function)

        # Check types
        for p in params.keys():
            tname = params[p]["type"]
            if not tname in self.types.keys():
                print "Error: Unknown type encountered:", tname
                sys.exit(7)

            params[p].setdefault('mode', 'IN')

        ## Roxygen to export the function
        internal = self.func[function].get("INTERNAL")
        if internal is None or internal == "False":
            out.write("#' @export\n")

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

        out.write(name)
        out.write(" <- function(")

        def do_par(pname):
            tname = params[pname]["type"]
            t = self.types[tname]
            default = ""
            header = pname.replace("_", ".")
            if "HEADER" in t:
                header = t["HEADER"]
            if header:
                header = header.replace("%I%", pname.replace("_", "."))
            else:
                header = ""
            if "default" in params[pname]:
                if "DEFAULT" in t and params[pname]["default"] in t["DEFAULT"]:
                    default = "=" + t["DEFAULT"][params[pname]["default"]]
                else:
                    default = "=" + params[pname]["default"]
            header = header + default
            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    header = header.replace("%I" + str(i + 1) + "%", dep)
            if re.search("%I[0-9]*%", header):
                print(
                    "Error: Missing HEADER dependency for "
                    + tname
                    + " "
                    + pname
                    + " in function "
                    + name
                )
            return header

        head = [do_par(n) for n, p in params.items() if p["mode"] in ["IN", "INOUT"]]
        head = [h for h in head if h != ""]
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
            tname = params[pname]["type"]
            t = self.types[tname]
            m = params[pname]["mode"]
            if m in ["IN", "INOUT"] and "INCONV" in t:
                if m in t["INCONV"]:
                    res = "  " + t["INCONV"][m]
                else:
                    res = "  " + t["INCONV"]
            else:
                res = ""
            res = res.replace("%I%", pname.replace("_", "."))

            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    res = res.replace("%I" + str(i + 1) + "%", dep)
            if re.search("%I[0-9]*%", res):
                print(
                    "Error: Missing IN dependency for "
                    + tname
                    + " "
                    + pname
                    + " in function "
                    + name
                )
            return res

        inconv = [do_par(n) for n in params.keys()]
        inconv = [i for i in inconv if i != ""]
        out.write("\n".join(inconv) + "\n\n")

        ## Function call
        ## This is a bit more difficult than INCONV. Here we supply
        ## each argument to the .Call function, if the argument has a
        ## 'CALL' field then it is used, otherwise we simply use its
        ## name.
        ## argument. Note that arguments with empty CALL fields are
        ## completely ignored, so giving an empty CALL field is
        ## different than not giving it at all.

        ## Function call
        def do_par(pname):
            t = self.types[params[pname]["type"]]
            call = pname.replace("_", ".")
            if "CALL" in t:
                call = t["CALL"]
                if call:
                    call = call.replace("%I%", pname.replace("_", "."))
                else:
                    call = ""
            return call

        out.write("  on.exit( .Call(C_R_igraph_finalizer) )\n")
        out.write("  # Function call\n")
        out.write("  res <- .Call(C_R_" + function + ", ")
        call = [do_par(n) for n, p in params.items() if p["mode"] in ["IN", "INOUT"]]
        call = [c for c in call if c != ""]
        out.write(", ".join(call))
        out.write(")\n")

        ## Output conversions
        def do_opar(pname, realname=None, iprefix=""):
            if realname is None:
                realname = pname
            tname = params[pname]["type"]
            t = self.types[tname]
            mode = params[pname]["mode"]
            if "OUTCONV" in t and mode in t["OUTCONV"]:
                outconv = "  " + t["OUTCONV"][mode]
            else:
                outconv = ""
            outconv = outconv.replace("%I%", iprefix + realname)

            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    outconv = outconv.replace("%I" + str(i + 1) + "%", dep)
            if re.search("%I[0-9]*%", outconv):
                print(outconv)
                print(self.deps)
                print(
                    "Error: Missing OUT dependency for "
                    + tname
                    + " "
                    + pname
                    + " in function "
                    + name
                )
            return re.sub("%I[0-9]+%", "", outconv)

        retpars = [n for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]

        if len(retpars) <= 1:
            outconv = [do_opar(n, "res") for n in params.keys()]
        else:
            outconv = [do_opar(n, iprefix="res$") for n in params.keys()]

        outconv = [o for o in outconv if o != ""]

        if len(retpars) == 0:
            # returning the return value of the function
            rt = self.types[self.func[function]["RETURN"]]
            if "OUTCONV" in rt:
                retconv = "  " + rt["OUTCONV"]["OUT"]
            else:
                retconv = ""
            retconv = retconv.replace("%I%", "res")
            # TODO: %I1% etc, is not handled here!
            ret = "\n".join(outconv) + "\n" + retconv + "\n"
        elif len(retpars) == 1:
            # returning a single output value
            ret = "\n".join(outconv) + "\n"
        else:
            # returning a list of output values
            None
            ret = "\n".join(outconv) + "\n"
        out.write(ret)

        ## Some graph attributes to add
        if "GATTR-R" in self.func[function].keys():
            gattrs = self.func[function]["GATTR-R"].split(",")
            gattrs = [ga.split(" IS ", 1) for ga in gattrs]
            sstr = "  res <- set.graph.attribute(res, '{name}', '{val}')\n"
            for ga in gattrs:
                aname = ga[0].strip()
                aval = ga[1].strip().replace("'", "\\'")
                out.write(sstr.format(name=aname, val=aval))

        ## Add some parameters as graph attributes
        if "GATTR-PARAM-R" in self.func[function].keys():
            pars = self.func[function]["GATTR-PARAM-R"].split(",")
            pars = [p.strip().replace("_", ".") for p in pars]
            sstr = "  res <- set.graph.attribute(res, '{par}', {par})\n"
            for p in pars:
                out.write(sstr.format(par=p))

        ## Set the class if requested
        if "CLASS-R" in self.func[function].keys():
            myclass = self.func[function]["CLASS-R"]
            out.write('  class(res) <- "' + myclass + '"\n')

        ## See if there is a postprocessor
        if "PP-R" in self.func[function].keys():
            pp = self.func[function]["PP-R"]
            out.write("  res <- " + pp + "(res)\n")

        out.write("  res\n}\n\n")


class RCCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate_function(self, function, out):

        # Ignore?
        if self.ignore(function):
            return

        params = self.parse_params(function)
        self.deps = self.parse_deps(function)

        # Check types
        for p in params.keys():
            tname = params[p]["type"]
            if not tname in self.types.keys():
                print("Error: Unknown type " + tname + " in " + function)
                return
            params[p].setdefault("mode", "IN")

        ## Compile the output
        ## This code generator is quite difficult, so we use different
        ## functions to generate the approprite chunks and then
        ## compile them together using a simple template.
        ## See the documentation of each chunk below.
        res = {}
        res["func"] = function
        res["header"] = self.chunk_header(function, params)
        res["decl"] = self.chunk_declaration(function, params)
        res["inconv"] = self.chunk_inconv(function, params)
        res["call"] = self.chunk_call(function, params)
        res["outconv"] = self.chunk_outconv(function, params)

        # Replace into the template
        text = (
            """
/*-------------------------------------------/
/ %(func)-42s /
/-------------------------------------------*/
%(header)s {
                                        /* Declarations */
%(decl)s
                                        /* Convert input */
%(inconv)s
                                        /* Call igraph */
%(call)s
                                        /* Convert output */
%(outconv)s

  UNPROTECT(1);
  return(result);
}\n"""
            % res
        )

        out.write(text)

    def chunk_header(self, function, params):
        """The header. All functions return with a 'SEXP', so this is
        easy. We just take the 'IN' and 'INOUT' arguments, all will
        have type SEXP, and concatenate them by commas. The function name
        is created by prefixing the original name with 'R_'."""

        def do_par(pname):
            t = self.types[params[pname]["type"]]
            if "HEADER" in t:
                if t["HEADER"]:
                    return t["HEADER"].replace("%I%", pname)
                else:
                    return ""
            else:
                return pname

        inout = [do_par(n) for n, p in params.items() if p["mode"] in ["IN", "INOUT"]]
        inout = ["SEXP " + n for n in inout if n != ""]
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
            cname = "c_" + pname
            t = self.types[params[pname]["type"]]
            if "DECL" in t:
                decl = "  " + t["DECL"]
            elif "CTYPE" in t:
                ctype = t["CTYPE"]
                if type(ctype) == dict:
                    mode = params[pname]["mode"]
                    decl = "  " + ctype[mode] + " " + cname + ";"
                else:
                    decl = "  " + ctype + " " + cname + ";"
            else:
                decl = ""
            return decl.replace("%C%", cname).replace("%I%", pname)

        inout = [do_par(n) for n in params.keys()]
        out = ["  SEXP " + n + ";" for n, p in params.items() if p["mode"] == "OUT"]

        retpars = [n for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]

        rt = self.types[self.func[function]["RETURN"]]
        if "DECL" in rt:
            retdecl = "  " + rt["DECL"]
        elif "CTYPE" in rt and len(retpars) == 0:
            ctype = rt["CTYPE"]
            if type(ctype) == dict:
                mode = params[pname]["mode"]
                retdecl = "  " + ctype[mode] + " " + "c_result;"
            else:
                retdecl = "  " + rt["CTYPE"] + " c_result;"
        else:
            retdecl = ""

        if len(retpars) <= 1:
            res = "\n".join(inout + out + [retdecl] + ["  SEXP result;"])
        else:
            res = "\n".join(inout + out + [retdecl] + ["  SEXP result, names;"])
        return res

    def chunk_inconv(self, function, params):
        """Input conversions. Not only for types with mode 'IN' and
        'INOUT', eg. for 'OUT' vector types we need to allocate the
        required memory here, do all the initializations, etc. Types
        without INCONV fields are ignored. The usual %C%, %I% is
        performed at the end.
        """

        def do_par(pname):
            cname = "c_" + pname
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "INCONV" in t and mode in t["INCONV"]:
                inconv = "  " + t["INCONV"][mode]
            else:
                inconv = ""

            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    inconv = inconv.replace("%C" + str(i + 1) + "%", "c_" + dep)

            return inconv.replace("%C%", cname).replace("%I%", pname)

        inconv = [do_par(n) for n in params.keys()]
        inconv = [i for i in inconv if i != ""]

        return "\n".join(inconv)

    def chunk_call(self, function, params):
        """Every single argument is included, independently of their
        mode. If a type has a 'CALL' field then that is used after the
        usual %C% and %I% substitutions, otherwise the standard 'c_'
        prefixed C argument name is used.
        """

        def docall(t, n):
            if type(t) == dict:
                mode = params[n]["mode"]
                if mode in t:
                    return t[mode]
                else:
                    return ""
            else:
                return t

        types = [self.types[params[n]["type"]] for n in params.keys()]
        call = list(map(
            lambda t, n: docall(t.get("CALL", "c_" + n), n), types, params.keys()
        ))
        call = list(map(
            lambda c, n: c.replace("%C%", "c_" + n).replace("%I%", n),
            call,
            params.keys(),
        ))
        retpars = [n for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]
        call = [c for c in call if c != ""]
        res = "  " + function + "(" + ", ".join(call) + ");\n"
        if len(retpars) == 0:
            res = "  c_result=" + res
        return res

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
            cname = "c_" + pname
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "OUTCONV" in t and mode in t["OUTCONV"]:
                outconv = "  " + t["OUTCONV"][mode]
            else:
                outconv = ""

            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    outconv = outconv.replace("%C" + str(i + 1) + "%", "c_" + dep)
            return outconv.replace("%C%", cname).replace("%I%", pname)

        outconv = [do_par(n) for n in params.keys()]
        outconv = [o for o in outconv if o != ""]

        retpars = [n for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]
        if len(retpars) == 0:
            # return the return value of the function
            rt = self.types[self.func[function]["RETURN"]]
            if "OUTCONV" in rt:
                retconv = "  " + rt["OUTCONV"]["OUT"]
            else:
                retconv = ""
            retconv = retconv.replace("%C%", "c_result").replace("%I%", "result")
            ret = "\n".join(outconv) + "\n" + retconv
        elif len(retpars) == 1:
            # return the single output value
            retconv = "  result=" + retpars[0] + ";"
            ret = "\n".join(outconv) + "\n" + retconv
        else:
            # create a list of output values
            sets = list(
                map(
                    lambda c, n: "  SET_VECTOR_ELT(result, " + str(c) + ", " + n + ");",
                    range(len(retpars)),
                    retpars,
                )
            )
            names = list(
                map(
                    lambda c, n: "  SET_STRING_ELT(names, "
                    + str(c)
                    + ', CREATE_STRING_VECTOR("'
                    + n
                    + '"));',
                    range(len(retpars)),
                    retpars,
                )
            )
            ret = "\n".join(
                [
                    "  PROTECT(result=NEW_LIST(" + str(len(retpars)) + "));",
                    "  PROTECT(names=NEW_CHARACTER(" + str(len(retpars)) + "));",
                ]
                + outconv
                + sets
                + names
                + ["  SET_NAMES(result, names);"]
                + ["  UNPROTECT(" + str(len(sets) + 1) + ");"]
            )

        return ret


################################################################################
# Java interface, experimental version using JNI (Java Native Interface)
# TODO: - everything :) This is just a PoC implementation.
################################################################################


class JavaCodeGenerator(CodeGenerator):
    """Class containing the common parts of JavaJavaCodeGenerator and
    JavaCCodeGenerator"""

    package = "net.sf.igraph"

    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def camelcase(s):
        """Returns a camelCase version of the given string (as used in Java
        libraries"""
        parts = s.split("_")
        result = [parts.pop(0)]
        for part in parts:
            result.append(part.capitalize())
        return "".join(result)

    camelcase = staticmethod(camelcase)

    def get_function_metadata(self, f, type_param="JAVATYPE"):
        """Returns metadata for the given function based on the parameters.
        f is the name of the function. The result is a dict with the following
        keys:

        - java_modifiers: Java modifiers to be used in the .java file
        - return_type: return type of the function
        - name: name of the function
        - argument_types: list of argument types
        - self_name: name of the "self" argument
        - is_static: whether the function is static
        - is_constructor: whether the function is a constructor
        """
        params = self.parse_params(f)
        is_static, is_constructor = False, False

        # We will collect data related to the current function in a dict
        data = {}
        data["name"] = self.func[f].get("NAME-JAVA", JavaCodeGenerator.camelcase(f[7:]))
        data["java_modifiers"] = ["public"]

        # Check parameter types to determine Java calling semantics
        types = {"IN": [], "OUT": [], "INOUT": []}
        for p in params.keys():
            types[params[p]["mode"]].append(params[p])

        if len(types["OUT"]) + len(types["INOUT"]) == 1:
            # If a single one is OUT or INOUT and all others are
            # INs, then this is our lucky day - the method fits the Java
            # semantics
            if len(types["OUT"]) > 0:
                return_type_name = types["OUT"][0]["type"]
            else:
                return_type_name = types["INOUT"][0]["type"]
        elif len(types["OUT"]) + len(types["INOUT"]) == 0 and "RETURN" in self.func[f]:
            # There are only input parameters and the return type is specified,
            # this also fits the Java semantics
            return_type_name = self.func[f]["RETURN"]
        else:
            raise StimulusError(
                "{}: calling convention unsupported yet".format(data["name"])
            )

        # Loop through the input parameters
        method_arguments = []
        found_self = False
        for p in params.keys():
            if params[p]["mode"] != "IN":
                continue
            type_name = params[p]["type"]
            if not found_self and type_name == "GRAPH":
                # this will be the 'self' argument
                found_self = True
                data["self_name"] = p
                continue
            tdesc = self.types.get(type_name, {})
            if type_param not in tdesc:
                raise StimulusError(
                    "{}: unknown input type {} (needs {}), skipping".format(
                        data["name"], type_name, type_param
                    )
                )
            method_arguments.append(" ".join([tdesc[type_param], p]))
        data["argument_types"] = method_arguments

        if not found_self:
            # Loop through INOUT arguments if we found no "self" yet
            for p in params.keys():
                if params[p]["mode"] == "INOUT" and params[p]["type"] == "GRAPH":
                    found_self = True
                    data["self_name"] = p
                    break

        tdesc = self.types.get(return_type_name, {})
        if type_param not in tdesc:
            raise StimulusError(
                "{}: unknown return type {}, skipping".format(
                    data["name"], return_type_name
                )
            )
        data["return_type"] = tdesc[type_param]

        if not found_self:
            data["java_modifiers"].append("static")
            data["name"] = data["name"][0].upper() + data["name"][1:]

        data["java_modifiers"] = " ".join(data["java_modifiers"])
        data["is_static"] = not found_self
        data["is_constructor"] = is_constructor

        return data


class JavaJavaCodeGenerator(JavaCodeGenerator):
    def __init__(self, func, types):
        JavaCodeGenerator.__init__(self, func, types)

    def generate(self, inputs, output):
        out = open(output, "w")

        if len(inputs) > 1:
            raise StimulusError("Java code generator supports only a single input")

        input = open(inputs[0], "rU")
        for line in input:
            if "%STIMULUS%" not in line:
                out.write(line)
                continue

            for f in self.func.keys():
                if self.ignore(f):
                    continue
                try:
                    func_metadata = self.get_function_metadata(f)
                    func_metadata["arguments"] = ", ".join(
                        func_metadata["argument_types"]
                    )
                    out.write(
                        "    %(java_modifiers)s native %(return_type)s %(name)s(%(arguments)s);\n"
                        % func_metadata
                    )
                except StimulusError as e:
                    out.write("    // %s\n" % str(e))

        out.close()


class JavaCCodeGenerator(JavaCodeGenerator):
    def __init__(self, func, types):
        JavaCodeGenerator.__init__(self, func, types)

    def generate_function(self, function, out):
        # Ignore?
        if self.ignore(function):
            return

        try:
            self.metadata = self.get_function_metadata(function, "CTYPE")
        except StimulusError as e:
            out.write("/* %s */\n" % str(e))
            return

        params = self.parse_params(function)
        self.deps = self.parse_deps(function)

        # Check types
        for p in params.keys():
            tname = params[p]["type"]
            if not tname in self.types.keys():
                print("Error: Unknown type " + tname + " in " + function)
                return
            params[p].setdefault("mode", "IN")

        ## Compile the output
        ## This code generator is quite difficult, so we use different
        ## functions to generate the approprite chunks and then
        ## compile them together using a simple template.
        ## See the documentation of each chunk below.
        try:
            res = {}
            res["func"] = function
            res["header"] = self.chunk_header(function, params)
            res["decl"] = self.chunk_declaration(function, params)
            res["before"] = self.chunk_before(function, params)
            res["inconv"] = self.chunk_inconv(function, params)
            res["call"] = self.chunk_call(function, params)
            res["outconv"] = self.chunk_outconv(function, params)
            res["after"] = self.chunk_after(function, params)
        except StimulusError as e:
            out.write("/* %s */\n" % str(e))
            return

        # Replace into the template
        text = (
            """
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

  return result;
}\n"""
            % res
        )

        out.write(text)

    def chunk_header(self, function, params):
        """The header.

        The name of the function is the igraph function name minus the
        igraph_ prefix, camelcased and prefixed with the underscored
        Java classname: net_sf_igraph_Graph_. The arguments
        are mapped from the JAVATYPE key of the type dict. Static
        methods also need a 'jclass cls' argument, ordinary methods
        need 'jobject jobj'. Besides that, the Java environment pointer
        is also passed.
        """
        data = self.get_function_metadata(function, "JAVATYPE")
        types = []

        data["funcname"] = "Java_%s_Graph_%s" % (
            self.package.replace(".", "_"),
            data["name"],
        )

        if data["is_static"]:
            data["argument_types"].insert(0, "jclass cls")
        else:
            data["argument_types"].insert(0, "jobject " + data["self_name"])
        data["argument_types"].insert(0, "JNIEnv *env")

        data["types"] = ", ".join(data["argument_types"])

        res = "JNIEXPORT %(return_type)s JNICALL %(funcname)s(%(types)s)" % data
        return res

    def chunk_declaration(self, function, params):
        """The declaration part of the function body

        There are a couple of things to declare. First a C type is
        needed for every argument, these will be supplied in the C
        igraph call. Then, all 'OUT' arguments need an appropriate variable as
        well, the result will be stored here. The return type
        of the C function also needs to be declared, that comes
        next. The result variable will contain the final result. Finally,
        if the method is not static but we are returning a new Graph object
        (e.g. in the case of igraph_linegraph), we need a jclass variable
        to store the Java class object."""

        def do_cpar(pname):
            cname = "c_" + pname
            t = self.types[params[pname]["type"]]
            if "CDECL" in t:
                decl = "  " + t["CDECL"]
            elif "CTYPE" in t:
                decl = "  " + t["CTYPE"] + " " + cname + ";"
            else:
                decl = ""
            return decl.replace("%C%", cname).replace("%I%", pname)

        def do_jpar(pname):
            jname = "j_" + pname
            t = self.types[params[pname]["type"]]
            if "JAVADECL" in t:
                decl = "  " + t["JAVADECL"]
            elif "JAVATYPE" in t:
                decl = "  " + t["JAVATYPE"] + " " + jname + ";"
            else:
                decl = ""
            return decl.replace("%J%", jname).replace("%I%", pname)

        inout = [do_cpar(n) for n in params.keys()]
        out = [do_jpar(n) for n, p in params.items() if p["mode"] == "OUT"]

        rt = self.types[self.func[function]["RETURN"]]
        if "CDECL" in rt:
            retdecl = "  " + rt["CDECL"]
        elif "CTYPE" in rt:
            retdecl = "  " + rt["CTYPE"] + " c__result;"
        else:
            retdecl = ""

        rnames = [n for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]
        jretdecl = ""
        if len(rnames) > 0:
            n = rnames[0]
            rtname = params[n]["type"]
        else:
            rtname = self.func[function]["RETURN"]
        rt = self.types[rtname]
        if "JAVADECL" in rt:
            jretdecl = "  " + rt["JAVADECL"]
        elif "JAVATYPE" in rt:
            jretdecl = "  " + rt["JAVATYPE"] + " result;"

        decls = inout + out + [retdecl, jretdecl]
        if not self.metadata["is_static"] and rtname == "GRAPH":
            self.metadata["need_class_decl"] = True
            decls.append(
                "  jclass cls = (*env)->GetObjectClass(env, %s);"
                % self.metadata["self_name"]
            )
        else:
            self.metadata["need_class_decl"] = False
        return "\n".join([i for i in decls if i != ""])

    def chunk_before(self, function, params):
        """We simply call Java_igraph_before"""
        return "  Java_igraph_before();"

    def chunk_inconv(self, function, params):
        """Input conversions. Not only for types with mode 'IN' and
        'INOUT', eg. for 'OUT' vector types we need to allocate the
        required memory here, do all the initializations, etc. Types
        without INCONV fields are ignored. The usual %C%, %I% is
        performed at the end.
        """

        def do_par(pname):
            cname = "c_" + pname
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "INCONV" in t and mode in t["INCONV"]:
                inconv = "  " + t["INCONV"][mode]
            else:
                inconv = ""

            if pname in self.deps.keys():
                deps = self.deps[pname]
                for i, dep in enumerate(deps):
                    inconv = inconv.replace("%C" + str(i + 1) + "%", "c_" + dep)

            return inconv.replace("%C%", cname).replace("%I%", pname)

        inconv = [do_par(n) for n in params.keys()]
        inconv = [i for i in inconv if i != ""]

        return "\n".join(inconv)

    def chunk_call(self, function, params):
        """Every single argument is included, independently of their
        mode. If a type has a 'CALL' field then that is used after the
        usual %C% and %I% substitutions, otherwise the standard 'c_'
        prefixed C argument name is used.
        """
        types = [self.types[params[n]["type"]] for n in params.keys()]
        call = list(map(lambda t, n: t.get("CALL", "c_" + n), types, params.keys()))
        call = list(map(
            lambda c, n: c.replace("%C%", "c_" + n).replace("%I%", n),
            call,
            params.keys(),
        ))
        lines = [
            "  if ((*env)->ExceptionCheck(env)) {",
            "    c__result = IGRAPH_EINVAL;",
            "  } else {",
            "    c__result = " + function + "(" + ", ".join(call) + ");",
            "  }",
        ]
        return "\n".join(lines)

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
        the returned type and this will be the result. The case of
        more than one 'OUT' and 'INOUT' arguments is not yet supported by
        the Java interface.
        """

        def do_par(pname):
            cname = "c_" + pname
            jname = "j_" + pname
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "OUTCONV" in t and mode in t["OUTCONV"]:
                outconv = "  " + t["OUTCONV"][mode]
            else:
                outconv = ""
            return outconv.replace("%C%", cname).replace("%I%", jname)

        outconv = [do_par(n) for n in params.keys()]
        outconv = [o for o in outconv if o != ""]

        retpars = [(n, p) for n, p in params.items() if p["mode"] in ["OUT", "INOUT"]]
        if len(retpars) == 0:
            # return the return value of the function
            rt = self.types[self.func[function]["RETURN"]]
            if "OUTCONV" in rt:
                retconv = "  " + rt["OUTCONV"]["OUT"]
            else:
                retconv = ""
            retconv = retconv.replace("%C%", "c__result").replace("%I%", "result")
            if len(retconv) > 0:
                outconv.append(retconv)
            ret = "\n".join(outconv)
        elif len(retpars) == 1:
            # return the single output value
            if retpars[0][1]["mode"] == "OUT":
                # OUT parameter
                retconv = "  result = j_" + retpars[0][0] + ";"
            else:
                # INOUT parameter
                retconv = "  result = " + retpars[0][0] + ";"
            outconv.append(retconv)

            outconv.insert(0, "if (c__result == 0) {")
            outconv.extend(["} else {", "  result = 0;", "}"])
            outconv = ["  %s" % line for line in outconv]
            ret = "\n".join(outconv)
        else:
            raise StimulusError(
                "{}: the case of multiple outputs not supported yet".format(function)
            )

        return ret

    def chunk_after(self, function, params):
        """We simply call Java_igraph_after"""
        return "  Java_igraph_after();"


################################################################################
# Shell interface, igraph functions directly from the command line
# TODO: - read/write default input/output from/to stdin/stdout
#       - short options
#       - prefixed output (?)
#       - default values depending on other parameters
#       - other input/output graph formats, to be controlled by
#         environment variables (?): IGRAPH_INGRAPH, IGRAPH_OUTGRAPH
################################################################################


class ShellLnCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate(self, inputs, output):
        out = open(output, "w")
        self.append_inputs(inputs, out)
        for f in self.func.keys():
            if self.ignore(f):
                continue
            out.write(f + "\n")
        out.close()


class ShellCodeGenerator(CodeGenerator):
    def __init__(self, func, types):
        CodeGenerator.__init__(self, func, types)

    def generate(self, inputs, output):
        out = open(output, "w")
        self.append_inputs(inputs, out)
        out.write("\n/* Function prototypes first */\n\n")

        for f in self.func.keys():
            if self.ignore(f):
                continue
            if "FLAGS" in self.func[f]:
                flags = self.func[f]["FLAGS"]
                flags = flags.split(",")
                flags = [flag.strip() for flag in flags]
            else:
                self.func[f]["FLAGS"] = []
            self.generate_prototype(f, out)

        out.write("\n/* The main function */\n\n")
        out.write("int main(int argc, char **argv) {\n\n")
        out.write("  const char *base=basename(argv[0]);\n\n  ")
        for f in self.func.keys():
            if self.ignore(f):
                continue
            out.write(
                'if (!strcasecmp(base, "'
                + f
                + '")) {\n    return shell_'
                + f
                + "(argc, argv);\n  } else "
            )
        out.write('{\n    printf("Unknown function, exiting\\n");\n')
        out.write("  }\n\n  shell_igraph_usage(argc, argv);\n  return 0;\n\n}\n")

        out.write("\n/* The functions themselves at last */\n")
        for f in self.func.keys():
            if self.ignore(f):
                continue
            self.generate_function(f, out)

        out.close()

    def generate_prototype(self, function, out):
        out.write("int shell_" + function + "(int argc, char **argv);\n")

    def generate_function(self, function, out):
        params = self.parse_params(function)

        # Check types, also enumerate them
        args = OrderedDict()
        for p in params.keys():
            tname = params[p]["type"]
            if not tname in self.types.keys():
                print "W: Unknown type encountered:", tname
                return

            params[p].setdefault('mode', 'IN')
            t=self.types[tname]
            mode=params[p]['mode']
            if 'INCONV' in t or 'OUTCONV' in t:
                args[p]=params[p].copy()
                args[p]['shell_no']=len(args)-1
                if mode=="INOUT":
                    args[p]['mode']='IN'
                    args[p+'-out']=params[p].copy()
                    args[p+'-out']['mode']='OUT'
                    args[p+'-out']['shell_no']=len(args)-1
                    if 'INCONV' not in t or 'IN' not in t['INCONV']:
                        print "Warning: no INCONV for type", tname, ", mode IN"
                    if 'OUTCONV' not in t or 'OUT' not in t['OUTCONV']:
                        print "Warning: no OUTCONV for type", tname, ", mode OUT"
            if mode =='IN' and ('INCONV' not in t or mode not in t['INCONV']):
                print "Warning: no INCONV for type", tname, ", mode", mode
            if mode == 'OUT' and ('OUTCONV' not in t or mode not in t['OUTCONV']):
                print "Warning: no OUTCONV for type", tname, ", mode", mode

        res={'nargs': len(args)}
        res['func']=function
        res['args']=self.chunk_args(function, args)
        res['decl']=self.chunk_decl(function, params)
        res['inconv']=self.chunk_inconv(function, args)
        res['call']=self.chunk_call(function, params)
        res['outconv']=self.chunk_outconv(function, args)
        res['default']=self.chunk_default(function, args)
        res['usage']=self.chunk_usage(function, args)
        text="""
/*-------------------------------------------/
/ %(func)-42s /
/-------------------------------------------*/
void shell_%(func)s_usage(char **argv) {
%(usage)s
  exit(1);
}

int shell_%(func)s(int argc, char **argv) {

%(decl)s

  int shell_seen[%(nargs)s];
  int shell_index=-1;
  struct option shell_options[]= { %(args)s
                                   { "help",no_argument,0,%(nargs)s },
                                   { 0,0,0,0 }
                                 };

  /* 0 - not seen, 1 - seen as argument, 2 - seen as default */
  memset(shell_seen, 0, %(nargs)s*sizeof(int));
%(default)s

  /* Parse arguments and read input */
  while (getopt_long(argc, argv, "", shell_options, &shell_index) != -1) {

    if (shell_index==-1) {
      exit(1);
    }

    if (shell_seen[shell_index]==1) {
      fprintf(stderr, "Error, `--%%s' argument given twice.\\n",
              shell_options[shell_index].name);
      exit(1);
    }
    shell_seen[shell_index]=1;
%(inconv)s
    shell_index=-1;
  }

  /* Check that we have all arguments */
  for (shell_index=0; shell_index<%(nargs)s; shell_index++) {
    if (!shell_seen[shell_index]) {
      fprintf(stderr, "Error, argument missing: `--%%s'.\\n",
              shell_options[shell_index].name);
      exit(1);
    }
  }

  /* Do the operation */
%(call)s

  /* Write the result */
%(outconv)s

  return 0;
}\n"""
            % res
        )
        out.write(text)

    def chunk_args(self, function, params):
        res = [
            ['"' + n + '"', "required_argument", "0", str(p["shell_no"])]
            for n, p in params.items()
        ]
        res = ["{ " + ",".join(e) + " }," for e in res]
        return "\n                                   ".join(res)

    def chunk_decl(self, function, params):
        def do_par(pname):
            t = self.types[params[pname]["type"]]
            if "DECL" in t:
                decl = "  " + t["DECL"].replace("%C%", pname)
            elif "CTYPE" in t:
                decl = "  " + t["CTYPE"] + " " + pname
            else:
                decl = ""
            if "default" in params[pname]:
                if "DEFAULT" in t and params[pname]["default"] in t["DEFAULT"]:
                    default = "=" + t["DEFAULT"][params[pname]["default"]]
                else:
                    default = "=" + params[pname]["default"]
            else:
                default = ""
            if decl:
                return decl + default + ";"
            else:
                return ""

        decl = [do_par(n) for n in params.keys()]
        inout = [
            "  char* shell_arg_" + n + "=0;"
            for n, p in params.items()
            if p["mode"] in ["INOUT", "OUT"]
        ]
        rt = self.types[self.func[function]["RETURN"]]
        if "DECL" in rt:
            retdecl = "  " + rt["DECL"]
        elif "CTYPE" in rt:
            retdecl = "  " + rt["CTYPE"] + " shell_result;"
        else:
            retdecl = ""

        if self.func[function]["RETURN"] != "ERROR":
            retchar = '  char *shell_arg_shell_result="-";'
        else:
            retchar = ""
        return "\n".join(decl + inout + [retdecl, retchar])

    def chunk_default(self, function, params):
        def do_par(pname):
            t = self.types[params[pname]["type"]]
            if "default" in params[pname]:
                res = "  shell_seen[" + str(params[pname]["shell_no"]) + "]=2;"
            else:
                res = ""
            return res

        res = [do_par(n) for n in params.keys()]
        res = [n for n in res if n != ""]
        return "\n".join(res)

    def chunk_inconv(self, function, params):
        def do_par(pname):
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "INCONV" in t and mode in t["INCONV"]:
                inconv = "" + t["INCONV"][mode]
            else:
                inconv = ""
            if pname.endswith("-out"):
                pname = pname[0:-4]
            return inconv.replace("%C%", pname)

        inconv = [
            "    case " + str(p["shell_no"]) + ": /* " + n + " */\n      " + do_par(n)
            for n, p in params.items()
        ]
        inconv = [n + "\n      break;" for n in inconv]
        inconv = ["".join(n) for n in inconv]
        text = (
            "\n    switch (shell_index) {\n"
            + "\n".join(inconv)
            + "\n    case "
            + str(len(inconv))
            + ":\n      shell_"
            + function
            + "_usage(argv);\n      break;"
            + "\n    default:\n      break;\n    }\n"
        )
        return text

    def chunk_call(self, function, params):
        types = [self.types[params[n]["type"]] for n in params.keys()]
        call = list(map(lambda t, n: t.get("CALL", n), types, params.keys()))
        call = list(map(lambda c, n: c.replace("%C%", n), call, params.keys()))
        return "  shell_result=" + function + "(" + ", ".join(call) + ");"

    def chunk_outconv(self, function, params):
        def do_par(pname):
            t = self.types[params[pname]["type"]]
            mode = params[pname]["mode"]
            if "OUTCONV" in t and mode in t["OUTCONV"]:
                outconv = "  " + t["OUTCONV"][mode]
            else:
                outconv = ""
            if pname.endswith("-out"):
                pname = pname[0:-4]
            return outconv.replace("%C%", pname)

        outconv = [do_par(n) for n in params.keys()]
        rt = self.types[self.func[function]["RETURN"]]
        if "OUTCONV" in rt and "OUT" in rt["OUTCONV"]:
            rtout = "  " + rt["OUTCONV"]["OUT"]
        else:
            rtout = ""
        outconv.append(rtout.replace("%C%", "shell_result"))
        outconv = [o for o in outconv if o != ""]
        return "\n".join(outconv)

    def chunk_usage(self, function, params):
        res = ["--" + n + "=<" + n + ">" for n in params.keys()]
        return '  printf("%s ' + " ".join(res) + '\\n", basename(argv[0]));'


################################################################################
if __name__ == "__main__":
    main()
