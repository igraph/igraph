/* glpapi19.c (stand-alone LP/MIP solver) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "glpapi.h"
#include "glpgmp.h"

struct csa
{     /* common storage area */
      glp_prob *prob;
      /* LP/MIP problem object */
      glp_bfcp bfcp;
      /* basis factorization control parameters */
      glp_smcp smcp;
      /* simplex method control parameters */
      glp_iptcp iptcp;
      /* interior-point method control parameters */
      glp_iocp iocp;
      /* integer optimizer control parameters */
      glp_tran *tran;
      /* model translator workspace */
      glp_graph *graph;
      /* network problem object */
      int format;
      /* problem file format: */
#define FMT_MPS_DECK    1  /* fixed MPS */
#define FMT_MPS_FILE    2  /* free MPS */
#define FMT_LP          3  /* CPLEX LP */
#define FMT_GLP         4  /* GLPK LP/MIP */
#define FMT_MATHPROG    5  /* MathProg */
#define FMT_MIN_COST    6  /* DIMACS min-cost flow */
#define FMT_MAX_FLOW    7  /* DIMACS maximum flow */
      const char *in_file;
      /* name of input problem file */
#define DATA_MAX 10
      /* maximal number of input data files */
      int ndf;
      /* number of input data files specified */
      const char *in_data[1+DATA_MAX];
      /* name(s) of input data file(s) */
      const char *out_dpy;
      /* name of output file to send display output; NULL means the
         display output is sent to the terminal */
      int seed;
      /* seed value to be passed to the MathProg translator; initially
         set to 1; 0x80000000 means the value is omitted */
      int solution;
      /* solution type flag: */
#define SOL_BASIC       1  /* basic */
#define SOL_INTERIOR    2  /* interior-point */
#define SOL_INTEGER     3  /* mixed integer */
      const char *in_res;
      /* name of input solution file in raw format */
      int dir;
      /* optimization direction flag:
         0       - not specified
         GLP_MIN - minimization
         GLP_MAX - maximization */
      int scale;
      /* automatic problem scaling flag */
      const char *out_sol;
      /* name of output solution file in printable format */
      const char *out_res;
      /* name of output solution file in raw format */
      const char *out_ranges;
      /* name of output file to write sensitivity analysis report */
      int check;
      /* input data checking flag; no solution is performed */
      const char *new_name;
      /* new name to be assigned to the problem */
      const char *out_mps;
      /* name of output problem file in fixed MPS format */
      const char *out_freemps;
      /* name of output problem file in free MPS format */
      const char *out_cpxlp;
      /* name of output problem file in CPLEX LP format */
      const char *out_glp;
      /* name of output problem file in GLPK format */
      const char *out_pb;
      /* name of output problem file in OPB format */
      const char *out_npb;
      /* name of output problem file in normalized OPB format */
      const char *log_file;
      /* name of output file to hardcopy terminal output */
      int crash;
      /* initial basis option: */
#define USE_STD_BASIS   1  /* use standard basis */
#define USE_ADV_BASIS   2  /* use advanced basis */
#define USE_CPX_BASIS   3  /* use Bixby's basis */
#define USE_INI_BASIS   4  /* use initial basis from ini_file */
      const char *ini_file;
      /* name of input file containing initial basis */
      int exact;
      /* flag to use glp_exact rather than glp_simplex */
      int xcheck;
      /* flag to check final basis with glp_exact */
      int nomip;
      /* flag to consider MIP as pure LP */
};

static void print_help(const char *my_name)
{     /* print help information */
      xprintf("Usage: %s [options...] filename\n", my_name);
      xprintf("\n");
      xprintf("General options:\n");
      xprintf("   --mps             read LP/MIP problem in fixed MPS fo"
         "rmat\n");
      xprintf("   --freemps         read LP/MIP problem in free MPS for"
         "mat (default)\n");
      xprintf("   --lp              read LP/MIP problem in CPLEX LP for"
         "mat\n");
      xprintf("   --glp             read LP/MIP problem in GLPK format "
         "\n");
      xprintf("   --math            read LP/MIP model written in GNU Ma"
         "thProg modeling\n");
      xprintf("                     language\n");
      xprintf("   -m filename, --model filename\n");
      xprintf("                     read model section and optional dat"
         "a section from\n");
      xprintf("                     filename (same as --math)\n");
      xprintf("   -d filename, --data filename\n");
      xprintf("                     read data section from filename (fo"
         "r --math only);\n");
      xprintf("                     if model file also has data section"
         ", it is ignored\n");
      xprintf("   -y filename, --display filename\n");
      xprintf("                     send display output to filename (fo"
         "r --math only);\n");
      xprintf("                     by default the output is sent to te"
         "rminal\n");
      xprintf("   --seed value      initialize pseudo-random number gen"
         "erator used in\n");
      xprintf("                     MathProg model with specified seed "
         "(any integer);\n");
      xprintf("                     if seed value is ?, some random see"
         "d will be used\n");
      xprintf("   --mincost         read min-cost flow problem in DIMAC"
         "S format\n");
      xprintf("   --maxflow         read maximum flow problem in DIMACS"
         " format\n");
      xprintf("   --simplex         use simplex method (default)\n");
      xprintf("   --interior        use interior point method (LP only)"
         "\n");
      xprintf("   -r filename, --read filename\n");
      xprintf("                     read solution from filename rather "
         "to find it with\n");
      xprintf("                     the solver\n");
      xprintf("   --min             minimization\n");
      xprintf("   --max             maximization\n");
      xprintf("   --scale           scale problem (default)\n");
      xprintf("   --noscale         do not scale problem\n");
      xprintf("   -o filename, --output filename\n");
      xprintf("                     write solution to filename in print"
         "able format\n");
      xprintf("   -w filename, --write filename\n");
      xprintf("                     write solution to filename in plain"
         " text format\n");
      xprintf("   --ranges filename\n");
      xprintf("                     write sensitivity analysis report t"
         "o filename in\n");
      xprintf("                     printable format (simplex only)\n");
      xprintf("   --tmlim nnn       limit solution time to nnn seconds "
         "\n");
      xprintf("   --memlim nnn      limit available memory to nnn megab"
         "ytes\n");
      xprintf("   --check           do not solve problem, check input d"
         "ata only\n");
      xprintf("   --name probname   change problem name to probname\n");
      xprintf("   --wmps filename   write problem to filename in fixed "
         "MPS format\n");
      xprintf("   --wfreemps filename\n");
      xprintf("                     write problem to filename in free M"
         "PS format\n");
      xprintf("   --wlp filename    write problem to filename in CPLEX "
         "LP format\n");
      xprintf("   --wglp filename   write problem to filename in GLPK f"
         "ormat\n");
#if 0
      xprintf("   --wpb filename    write problem to filename in OPB fo"
         "rmat\n");
      xprintf("   --wnpb filename   write problem to filename in normal"
         "ized OPB format\n");
#endif
      xprintf("   --log filename    write copy of terminal output to fi"
         "lename\n");
      xprintf("   -h, --help        display this help information and e"
         "xit\n");
      xprintf("   -v, --version     display program version and exit\n")
         ;
      xprintf("\n");
      xprintf("LP basis factorization options:\n");
      xprintf("   --luf             LU + Forrest-Tomlin update\n");
      xprintf("                     (faster, less stable; default)\n");
      xprintf("   --cbg             LU + Schur complement + Bartels-Gol"
         "ub update\n");
      xprintf("                     (slower, more stable)\n");
      xprintf("   --cgr             LU + Schur complement + Givens rota"
         "tion update\n");
      xprintf("                     (slower, more stable)\n");
      xprintf("\n");
      xprintf("Options specific to simplex solver:\n");
      xprintf("   --primal          use primal simplex (default)\n");
      xprintf("   --dual            use dual simplex\n");
      xprintf("   --std             use standard initial basis of all s"
         "lacks\n");
      xprintf("   --adv             use advanced initial basis (default"
         ")\n");
      xprintf("   --bib             use Bixby's initial basis\n");
      xprintf("   --ini filename    use as initial basis previously sav"
         "ed with -w\n");
      xprintf("                     (disables LP presolver)\n");
      xprintf("   --steep           use steepest edge technique (defaul"
         "t)\n");
      xprintf("   --nosteep         use standard \"textbook\" pricing\n"
         );
      xprintf("   --relax           use Harris' two-pass ratio test (de"
         "fault)\n");
      xprintf("   --norelax         use standard \"textbook\" ratio tes"
         "t\n");
      xprintf("   --presol          use presolver (default; assumes --s"
         "cale and --adv)\n");
      xprintf("   --nopresol        do not use presolver\n");
      xprintf("   --exact           use simplex method based on exact a"
         "rithmetic\n");
      xprintf("   --xcheck          check final basis using exact arith"
         "metic\n");
      xprintf("\n");
      xprintf("Options specific to interior-point solver:\n");
      xprintf("   --nord            use natural (original) ordering\n");
      xprintf("   --qmd             use quotient minimum degree orderin"
         "g\n");
      xprintf("   --amd             use approximate minimum degree orde"
         "ring (default)\n");
      xprintf("   --symamd          use approximate minimum degree orde"
         "ring\n");
      xprintf("\n");
      xprintf("Options specific to MIP solver:\n");
      xprintf("   --nomip           consider all integer variables as c"
         "ontinuous\n");
      xprintf("                     (allows solving MIP as pure LP)\n");
      xprintf("   --first           branch on first integer variable\n")
         ;
      xprintf("   --last            branch on last integer variable\n");
      xprintf("   --mostf           branch on most fractional variable "
         "\n");
      xprintf("   --drtom           branch using heuristic by Driebeck "
         "and Tomlin\n");
      xprintf("                     (default)\n");
      xprintf("   --pcost           branch using hybrid pseudocost heur"
         "istic (may be\n");
      xprintf("                     useful for hard instances)\n");
      xprintf("   --dfs             backtrack using depth first search "
         "\n");
      xprintf("   --bfs             backtrack using breadth first searc"
         "h\n");
      xprintf("   --bestp           backtrack using the best projection"
         " heuristic\n");
      xprintf("   --bestb           backtrack using node with best loca"
         "l bound\n");
      xprintf("                     (default)\n");
      xprintf("   --intopt          use MIP presolver (default)\n");
      xprintf("   --nointopt        do not use MIP presolver\n");
      xprintf("   --binarize        replace general integer variables b"
         "y binary ones\n");
      xprintf("                     (assumes --intopt)\n");
      xprintf("   --fpump           apply feasibility pump heuristic\n")
         ;
      xprintf("   --gomory          generate Gomory's mixed integer cut"
         "s\n");
      xprintf("   --mir             generate MIR (mixed integer roundin"
         "g) cuts\n");
      xprintf("   --cover           generate mixed cover cuts\n");
      xprintf("   --clique          generate clique cuts\n");
      xprintf("   --cuts            generate all cuts above\n");
      xprintf("   --mipgap tol      set relative mip gap tolerance to t"
         "ol\n");
      xprintf("\n");
      xprintf("For description of the MPS and CPLEX LP formats see Refe"
         "rence Manual.\n");
      xprintf("For description of the modeling language see \"GLPK: Mod"
         "eling Language\n");
      xprintf("GNU MathProg\". Both documents are included in the GLPK "
         "distribution.\n");
      xprintf("\n");
      xprintf("See GLPK web page at <http://www.gnu.org/software/glpk/g"
         "lpk.html>.\n");
      xprintf("\n");
      xprintf("Please report bugs to <bug-glpk@gnu.org>.\n");
      return;
}

static void print_version(int briefly)
{     /* print version information */
      xprintf("GLPSOL: GLPK LP/MIP Solver, v%s\n", glp_version());
      if (briefly) goto done;
      xprintf("\n");
      xprintf("Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, "
         "2007, 2008,\n");
      xprintf("2009, 2010 Andrew Makhorin, Department for Applied Infor"
         "matics, Moscow\n");
      xprintf("Aviation Institute, Moscow, Russia. All rights reserved."
         "\n");
      xprintf("\n");
      xprintf("This program has ABSOLUTELY NO WARRANTY.\n");
      xprintf("\n");
      xprintf("This program is free software; you may re-distribute it "
         "under the terms\n");
      xprintf("of the GNU General Public License version 3 or later.\n")
         ;
done: return;
}

static int parse_cmdline(struct csa *csa, int argc, const char *argv[])
{     /* parse command-line parameters */
      int k;
#define p(str) (strcmp(argv[k], str) == 0)
      for (k = 1; k < argc; k++)
      {  if (p("--mps"))
            csa->format = FMT_MPS_DECK;
         else if (p("--freemps"))
            csa->format = FMT_MPS_FILE;
         else if (p("--lp") || p("--cpxlp"))
            csa->format = FMT_LP;
         else if (p("--glp"))
            csa->format = FMT_GLP;
         else if (p("--math") || p("-m") || p("--model"))
            csa->format = FMT_MATHPROG;
         else if (p("-d") || p("--data"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No input data file specified\n");
               return 1;
            }
            if (csa->ndf == DATA_MAX)
            {  xprintf("Too many input data files\n");
               return 1;
            }
            csa->in_data[++(csa->ndf)] = argv[k];
         }
         else if (p("-y") || p("--display"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No display output file specified\n");
               return 1;
            }
            if (csa->out_dpy != NULL)
            {  xprintf("Only one display output file allowed\n");
               return 1;
            }
            csa->out_dpy = argv[k];
         }
         else if (p("--seed"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' ||
               argv[k][0] == '-' && !isdigit((unsigned char)argv[k][1]))
            {  xprintf("No seed value specified\n");
               return 1;
            }
            if (strcmp(argv[k], "?") == 0)
               csa->seed = 0x80000000;
            else if (str2int(argv[k], &csa->seed))
            {  xprintf("Invalid seed value `%s'\n", argv[k]);
               return 1;
            }
         }
         else if (p("--mincost"))
            csa->format = FMT_MIN_COST;
         else if (p("--maxflow"))
            csa->format = FMT_MAX_FLOW;
         else if (p("--simplex"))
            csa->solution = SOL_BASIC;
         else if (p("--interior"))
            csa->solution = SOL_INTERIOR;
#if 1 /* 28/V-2010 */
         else if (p("--alien"))
            csa->iocp.alien = GLP_ON;
#endif
         else if (p("-r") || p("--read"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No input solution file specified\n");
               return 1;
            }
            if (csa->in_res != NULL)
            {  xprintf("Only one input solution file allowed\n");
               return 1;
            }
            csa->in_res = argv[k];
         }
         else if (p("--min"))
            csa->dir = GLP_MIN;
         else if (p("--max"))
            csa->dir = GLP_MAX;
         else if (p("--scale"))
            csa->scale = 1;
         else if (p("--noscale"))
            csa->scale = 0;
         else if (p("-o") || p("--output"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No output solution file specified\n");
               return 1;
            }
            if (csa->out_sol != NULL)
            {  xprintf("Only one output solution file allowed\n");
               return 1;
            }
            csa->out_sol = argv[k];
         }
         else if (p("-w") || p("--write"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No output solution file specified\n");
               return 1;
            }
            if (csa->out_res != NULL)
            {  xprintf("Only one output solution file allowed\n");
               return 1;
            }
            csa->out_res = argv[k];
         }
         else if (p("--ranges") || p("--bounds"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No output file specified to write sensitivity a"
                  "nalysis report\n");
               return 1;
            }
            if (csa->out_ranges != NULL)
            {  xprintf("Only one output file allowed to write sensitivi"
                  "ty analysis report\n");
               return 1;
            }
            csa->out_ranges = argv[k];
         }
         else if (p("--tmlim"))
         {  int tm_lim;
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No time limit specified\n");
               return 1;
            }
            if (str2int(argv[k], &tm_lim) || tm_lim < 0)
            {  xprintf("Invalid time limit `%s'\n", argv[k]);
               return 1;
            }
            if (tm_lim <= INT_MAX / 1000)
               csa->smcp.tm_lim = csa->iocp.tm_lim = 1000 * tm_lim;
            else
               csa->smcp.tm_lim = csa->iocp.tm_lim = INT_MAX;
         }
         else if (p("--memlim"))
         {  int mem_lim;
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No memory limit specified\n");
               return 1;
            }
            if (str2int(argv[k], &mem_lim) || mem_lim < 1)
            {  xprintf("Invalid memory limit `%s'\n", argv[k]);
               return 1;
            }
            glp_mem_limit(mem_lim);
         }
         else if (p("--check"))
            csa->check = 1;
         else if (p("--name"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No problem name specified\n");
               return 1;
            }
            if (csa->new_name != NULL)
            {  xprintf("Only one problem name allowed\n");
               return 1;
            }
            csa->new_name = argv[k];
         }
         else if (p("--wmps"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No fixed MPS output file specified\n");
               return 1;
            }
            if (csa->out_mps != NULL)
            {  xprintf("Only one fixed MPS output file allowed\n");
               return 1;
            }
            csa->out_mps = argv[k];
         }
         else if (p("--wfreemps"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No free MPS output file specified\n");
               return 1;
            }
            if (csa->out_freemps != NULL)
            {  xprintf("Only one free MPS output file allowed\n");
               return 1;
            }
            csa->out_freemps = argv[k];
         }
         else if (p("--wlp") || p("--wcpxlp") || p("--wlpt"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No CPLEX LP output file specified\n");
               return 1;
            }
            if (csa->out_cpxlp != NULL)
            {  xprintf("Only one CPLEX LP output file allowed\n");
               return 1;
            }
            csa->out_cpxlp = argv[k];
         }
         else if (p("--wglp"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No GLPK LP/MIP output file specified\n");
               return 1;
            }
            if (csa->out_glp != NULL)
            {  xprintf("Only one GLPK LP/MIP output file allowed\n");
               return 1;
            }
            csa->out_glp = argv[k];
         }
         else if (p("--wpb"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No problem output file specified\n");
               return 1;
            }
            if (csa->out_pb != NULL)
            {  xprintf("Only one OPB output file allowed\n");
               return 1;
            }
            csa->out_pb = argv[k];
         }
         else if (p("--wnpb"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No problem output file specified\n");
               return 1;
            }
            if (csa->out_npb != NULL)
            {  xprintf("Only one normalized OPB output file allowed\n");
               return 1;
            }
            csa->out_npb = argv[k];
         }
         else if (p("--log"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No log file specified\n");
               return 1;
            }
            if (csa->log_file != NULL)
            {  xprintf("Only one log file allowed\n");
               return 1;
            }
            csa->log_file = argv[k];
         }
         else if (p("-h") || p("--help"))
         {  print_help(argv[0]);
            return -1;
         }
         else if (p("-v") || p("--version"))
         {  print_version(0);
            return -1;
         }
         else if (p("--luf"))
            csa->bfcp.type = GLP_BF_FT;
         else if (p("--cbg"))
            csa->bfcp.type = GLP_BF_BG;
         else if (p("--cgr"))
            csa->bfcp.type = GLP_BF_GR;
         else if (p("--primal"))
            csa->smcp.meth = GLP_PRIMAL;
         else if (p("--dual"))
            csa->smcp.meth = GLP_DUAL;
         else if (p("--std"))
            csa->crash = USE_STD_BASIS;
         else if (p("--adv"))
            csa->crash = USE_ADV_BASIS;
         else if (p("--bib"))
            csa->crash = USE_CPX_BASIS;
         else if (p("--ini"))
         {  csa->crash = USE_INI_BASIS;
            csa->smcp.presolve = GLP_OFF;
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No initial basis file specified\n");
               return 1;
            }
            if (csa->ini_file != NULL)
            {  xprintf("Only one initial basis file allowed\n");
               return 1;
            }
            csa->ini_file = argv[k];
         }
         else if (p("--steep"))
            csa->smcp.pricing = GLP_PT_PSE;
         else if (p("--nosteep"))
            csa->smcp.pricing = GLP_PT_STD;
         else if (p("--relax"))
            csa->smcp.r_test = GLP_RT_HAR;
         else if (p("--norelax"))
            csa->smcp.r_test = GLP_RT_STD;
         else if (p("--presol"))
            csa->smcp.presolve = GLP_ON;
         else if (p("--nopresol"))
            csa->smcp.presolve = GLP_OFF;
         else if (p("--exact"))
            csa->exact = 1;
         else if (p("--xcheck"))
            csa->xcheck = 1;
         else if (p("--nord"))
            csa->iptcp.ord_alg = GLP_ORD_NONE;
         else if (p("--qmd"))
            csa->iptcp.ord_alg = GLP_ORD_QMD;
         else if (p("--amd"))
            csa->iptcp.ord_alg = GLP_ORD_AMD;
         else if (p("--symamd"))
            csa->iptcp.ord_alg = GLP_ORD_SYMAMD;
         else if (p("--nomip"))
            csa->nomip = 1;
         else if (p("--first"))
            csa->iocp.br_tech = GLP_BR_FFV;
         else if (p("--last"))
            csa->iocp.br_tech = GLP_BR_LFV;
         else if (p("--drtom"))
            csa->iocp.br_tech = GLP_BR_DTH;
         else if (p("--mostf"))
            csa->iocp.br_tech = GLP_BR_MFV;
         else if (p("--pcost"))
            csa->iocp.br_tech = GLP_BR_PCH;
         else if (p("--dfs"))
            csa->iocp.bt_tech = GLP_BT_DFS;
         else if (p("--bfs"))
            csa->iocp.bt_tech = GLP_BT_BFS;
         else if (p("--bestp"))
            csa->iocp.bt_tech = GLP_BT_BPH;
         else if (p("--bestb"))
            csa->iocp.bt_tech = GLP_BT_BLB;
         else if (p("--intopt"))
            csa->iocp.presolve = GLP_ON;
         else if (p("--nointopt"))
            csa->iocp.presolve = GLP_OFF;
         else if (p("--binarize"))
            csa->iocp.presolve = csa->iocp.binarize = GLP_ON;
         else if (p("--fpump"))
            csa->iocp.fp_heur = GLP_ON;
         else if (p("--gomory"))
            csa->iocp.gmi_cuts = GLP_ON;
         else if (p("--mir"))
            csa->iocp.mir_cuts = GLP_ON;
         else if (p("--cover"))
            csa->iocp.cov_cuts = GLP_ON;
         else if (p("--clique"))
            csa->iocp.clq_cuts = GLP_ON;
         else if (p("--cuts"))
            csa->iocp.gmi_cuts = csa->iocp.mir_cuts =
            csa->iocp.cov_cuts = csa->iocp.clq_cuts = GLP_ON;
         else if (p("--mipgap"))
         {  double mip_gap;
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  xprintf("No relative gap tolerance specified\n");
               return 1;
            }
            if (str2num(argv[k], &mip_gap) || mip_gap < 0.0)
            {  xprintf("Invalid relative mip gap tolerance `%s'\n",
                  argv[k]);
               return 1;
            }
            csa->iocp.mip_gap = mip_gap;
         }
         else if (argv[k][0] == '-' ||
                 (argv[k][0] == '-' && argv[k][1] == '-'))
         {  xprintf("Invalid option `%s'; try %s --help\n",
               argv[k], argv[0]);
            return 1;
         }
         else
         {  if (csa->in_file != NULL)
            {  xprintf("Only one input problem file allowed\n");
               return 1;
            }
            csa->in_file = argv[k];
         }
      }
#undef p
      return 0;
}

typedef struct { double rhs, pi; } v_data;
typedef struct { double low, cap, cost, x; } a_data;

int glp_main(int argc, const char *argv[])
{     /* stand-alone LP/MIP solver */
      struct csa _csa, *csa = &_csa;
      int ret;
      glp_long start;
      /* perform initialization */
      csa->prob = glp_create_prob();
      glp_get_bfcp(csa->prob, &csa->bfcp);
      glp_init_smcp(&csa->smcp);
      csa->smcp.presolve = GLP_ON;
      glp_init_iptcp(&csa->iptcp);
      glp_init_iocp(&csa->iocp);
      csa->iocp.presolve = GLP_ON;
      csa->tran = NULL;
      csa->graph = NULL;
      csa->format = FMT_MPS_FILE;
      csa->in_file = NULL;
      csa->ndf = 0;
      csa->out_dpy = NULL;
      csa->seed = 1;
      csa->solution = SOL_BASIC;
      csa->in_res = NULL;
      csa->dir = 0;
      csa->scale = 1;
      csa->out_sol = NULL;
      csa->out_res = NULL;
      csa->out_ranges = NULL;
      csa->check = 0;
      csa->new_name = NULL;
      csa->out_mps = NULL;
      csa->out_freemps = NULL;
      csa->out_cpxlp = NULL;
      csa->out_glp = NULL;
      csa->out_pb = NULL;
      csa->out_npb = NULL;
      csa->log_file = NULL;
      csa->crash = USE_ADV_BASIS;
      csa->ini_file = NULL;
      csa->exact = 0;
      csa->xcheck = 0;
      csa->nomip = 0;
      /* parse command-line parameters */
      ret = parse_cmdline(csa, argc, argv);
      if (ret < 0)
      {  ret = EXIT_SUCCESS;
         goto done;
      }
      if (ret > 0)
      {  ret = EXIT_FAILURE;
         goto done;
      }
      /*--------------------------------------------------------------*/
      /* remove all output files specified in the command line */
      if (csa->out_dpy != NULL) remove(csa->out_dpy);
      if (csa->out_sol != NULL) remove(csa->out_sol);
      if (csa->out_res != NULL) remove(csa->out_res);
      if (csa->out_ranges != NULL) remove(csa->out_ranges);
      if (csa->out_mps != NULL) remove(csa->out_mps);
      if (csa->out_freemps != NULL) remove(csa->out_freemps);
      if (csa->out_cpxlp != NULL) remove(csa->out_cpxlp);
      if (csa->out_glp != NULL) remove(csa->out_glp);
      if (csa->out_pb != NULL) remove(csa->out_pb);
      if (csa->out_npb != NULL) remove(csa->out_npb);
      if (csa->log_file != NULL) remove(csa->log_file);
      /*--------------------------------------------------------------*/
      /* open log file, if required */
      if (csa->log_file != NULL)
      {  if (glp_open_tee(csa->log_file))
         {  xprintf("Unable to create log file\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /*--------------------------------------------------------------*/
      /* print version information */
      print_version(1);
      /*--------------------------------------------------------------*/
      /* print parameters specified in the command line */
      if (argc > 1)
      {  int k, len = INT_MAX;
         xprintf("Parameter(s) specified in the command line:");
         for (k = 1; k < argc; k++)
         {  if (len > 72)
               xprintf("\n"), len = 0;
            xprintf(" %s", argv[k]);
            len += 1 + strlen(argv[k]);
         }
         xprintf("\n");
      }
      /*--------------------------------------------------------------*/
      /* read problem data from the input file */
      if (csa->in_file == NULL)
      {  xprintf("No input problem file specified; try %s --help\n",
            argv[0]);
         ret = EXIT_FAILURE;
         goto done;
      }
      if (csa->format == FMT_MPS_DECK)
      {  ret = glp_read_mps(csa->prob, GLP_MPS_DECK, NULL,
            csa->in_file);
         if (ret != 0)
err1:    {  xprintf("MPS file processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      else if (csa->format == FMT_MPS_FILE)
      {  ret = glp_read_mps(csa->prob, GLP_MPS_FILE, NULL,
            csa->in_file);
         if (ret != 0) goto err1;
      }
      else if (csa->format == FMT_LP)
      {  ret = glp_read_lp(csa->prob, NULL, csa->in_file);
         if (ret != 0)
         {  xprintf("CPLEX LP file processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      else if (csa->format == FMT_GLP)
      {  ret = glp_read_prob(csa->prob, 0, csa->in_file);
         if (ret != 0)
         {  xprintf("GLPK LP/MIP file processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      else if (csa->format == FMT_MATHPROG)
      {  int k;
         /* allocate the translator workspace */
         csa->tran = glp_mpl_alloc_wksp();
         /* set seed value */
         if (csa->seed == 0x80000000)
         {  csa->seed = glp_time().lo;
            xprintf("Seed value %d will be used\n", csa->seed);
         }
         _glp_mpl_init_rand(csa->tran, csa->seed);
         /* read model section and optional data section */
         if (glp_mpl_read_model(csa->tran, csa->in_file, csa->ndf > 0))
err2:    {  xprintf("MathProg model processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
         /* read optional data section(s), if necessary */
         for (k = 1; k <= csa->ndf; k++)
         {  if (glp_mpl_read_data(csa->tran, csa->in_data[k]))
               goto err2;
         }
         /* generate the model */
         if (glp_mpl_generate(csa->tran, csa->out_dpy)) goto err2;
         /* build the problem instance from the model */
         glp_mpl_build_prob(csa->tran, csa->prob);
      }
      else if (csa->format == FMT_MIN_COST)
      {  csa->graph = glp_create_graph(sizeof(v_data), sizeof(a_data));
         ret = glp_read_mincost(csa->graph, offsetof(v_data, rhs),
            offsetof(a_data, low), offsetof(a_data, cap),
            offsetof(a_data, cost), csa->in_file);
         if (ret != 0)
         {  xprintf("DIMACS file processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
         glp_mincost_lp(csa->prob, csa->graph, GLP_ON,
            offsetof(v_data, rhs), offsetof(a_data, low),
            offsetof(a_data, cap), offsetof(a_data, cost));
         glp_set_prob_name(csa->prob, csa->in_file);
      }
      else if (csa->format == FMT_MAX_FLOW)
      {  int s, t;
         csa->graph = glp_create_graph(sizeof(v_data), sizeof(a_data));
         ret = glp_read_maxflow(csa->graph, &s, &t,
            offsetof(a_data, cap), csa->in_file);
         if (ret != 0)
         {  xprintf("DIMACS file processing error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
         glp_maxflow_lp(csa->prob, csa->graph, GLP_ON, s, t,
            offsetof(a_data, cap));
         glp_set_prob_name(csa->prob, csa->in_file);
      }
      else
         xassert(csa != csa);
      /*--------------------------------------------------------------*/
      /* change problem name, if required */
      if (csa->new_name != NULL)
         glp_set_prob_name(csa->prob, csa->new_name);
      /* change optimization direction, if required */
      if (csa->dir != 0)
         glp_set_obj_dir(csa->prob, csa->dir);
      /* sort elements of the constraint matrix */
      glp_sort_matrix(csa->prob);
      /*--------------------------------------------------------------*/
      /* write problem data in fixed MPS format, if required */
      if (csa->out_mps != NULL)
      {  ret = glp_write_mps(csa->prob, GLP_MPS_DECK, NULL,
            csa->out_mps);
         if (ret != 0)
         {  xprintf("Unable to write problem in fixed MPS format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem data in free MPS format, if required */
      if (csa->out_freemps != NULL)
      {  ret = glp_write_mps(csa->prob, GLP_MPS_FILE, NULL,
            csa->out_freemps);
         if (ret != 0)
         {  xprintf("Unable to write problem in free MPS format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem data in CPLEX LP format, if required */
      if (csa->out_cpxlp != NULL)
      {  ret = glp_write_lp(csa->prob, NULL, csa->out_cpxlp);
         if (ret != 0)
         {  xprintf("Unable to write problem in CPLEX LP format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem data in GLPK format, if required */
      if (csa->out_glp != NULL)
      {  ret = glp_write_prob(csa->prob, 0, csa->out_glp);
         if (ret != 0)
         {  xprintf("Unable to write problem in GLPK format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem data in OPB format, if required */
      if (csa->out_pb != NULL)
      {  ret = lpx_write_pb(csa->prob, csa->out_pb, 0, 0);
         if (ret != 0)
         {  xprintf("Unable to write problem in OPB format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem data in normalized OPB format, if required */
      if (csa->out_npb != NULL)
      {  ret = lpx_write_pb(csa->prob, csa->out_npb, 1, 1);
         if (ret != 0)
         {  xprintf(
               "Unable to write problem in normalized OPB format\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /*--------------------------------------------------------------*/
      /* if only problem data check is required, skip computations */
      if (csa->check)
      {  ret = EXIT_SUCCESS;
         goto done;
      }
      /*--------------------------------------------------------------*/
      /* determine the solution type */
      if (!csa->nomip &&
          glp_get_num_int(csa->prob) + glp_get_num_bin(csa->prob) > 0)
      {  if (csa->solution == SOL_INTERIOR)
         {  xprintf("Interior-point method is not able to solve MIP pro"
               "blem; use --simplex\n");
            ret = EXIT_FAILURE;
            goto done;
         }
         csa->solution = SOL_INTEGER;
      }
      /*--------------------------------------------------------------*/
      /* if solution is provided, read it and skip computations */
      if (csa->in_res != NULL)
      {  if (csa->solution == SOL_BASIC)
            ret = glp_read_sol(csa->prob, csa->in_res);
         else if (csa->solution == SOL_INTERIOR)
            ret = glp_read_ipt(csa->prob, csa->in_res);
         else if (csa->solution == SOL_INTEGER)
            ret = glp_read_mip(csa->prob, csa->in_res);
         else
            xassert(csa != csa);
         if (ret != 0)
         {  xprintf("Unable to read problem solution\n");
            ret = EXIT_FAILURE;
            goto done;
         }
         goto skip;
      }
      /*--------------------------------------------------------------*/
      /* scale the problem data, if required */
      if (csa->scale)
      {  if (csa->solution == SOL_BASIC && !csa->smcp.presolve ||
             csa->solution == SOL_INTERIOR ||
             csa->solution == SOL_INTEGER && !csa->iocp.presolve)
            glp_scale_prob(csa->prob, GLP_SF_AUTO);
      }
      /*--------------------------------------------------------------*/
      /* construct starting LP basis */
      if (csa->solution == SOL_BASIC && !csa->smcp.presolve ||
          csa->solution == SOL_INTEGER && !csa->iocp.presolve)
      {  if (csa->crash == USE_STD_BASIS)
            glp_std_basis(csa->prob);
         else if (csa->crash == USE_ADV_BASIS)
            glp_adv_basis(csa->prob, 0);
         else if (csa->crash == USE_CPX_BASIS)
            glp_cpx_basis(csa->prob);
         else if (csa->crash == USE_INI_BASIS)
         {  ret = glp_read_sol(csa->prob, csa->ini_file);
            if (ret != 0)
            {  xprintf("Unable to read initial basis\n");
               ret = EXIT_FAILURE;
               goto done;
            }
         }
         else
            xassert(csa != csa);
      }
      /*--------------------------------------------------------------*/
      /* solve the problem */
      start = xtime();
      if (csa->solution == SOL_BASIC)
      {  if (!csa->exact)
         {  glp_set_bfcp(csa->prob, &csa->bfcp);
            glp_simplex(csa->prob, &csa->smcp);
            if (csa->xcheck)
            {  if (csa->smcp.presolve &&
                   glp_get_status(csa->prob) != GLP_OPT)
                  xprintf("If you need to check final basis for non-opt"
                     "imal solution, use --nopresol\n");
               else
                  glp_exact(csa->prob, &csa->smcp);
            }
            if (csa->out_sol != NULL || csa->out_res != NULL)
            {  if (csa->smcp.presolve &&
                   glp_get_status(csa->prob) != GLP_OPT)
               xprintf("If you need actual output for non-optimal solut"
                  "ion, use --nopresol\n");
            }
         }
         else
            glp_exact(csa->prob, &csa->smcp);
      }
      else if (csa->solution == SOL_INTERIOR)
         glp_interior(csa->prob, &csa->iptcp);
      else if (csa->solution == SOL_INTEGER)
      {  if (!csa->iocp.presolve)
         {  glp_set_bfcp(csa->prob, &csa->bfcp);
            glp_simplex(csa->prob, &csa->smcp);
         }
#if 0
         csa->iocp.msg_lev = GLP_MSG_DBG;
         csa->iocp.pp_tech = GLP_PP_NONE;
#endif
         glp_intopt(csa->prob, &csa->iocp);
      }
      else
         xassert(csa != csa);
      /*--------------------------------------------------------------*/
      /* display statistics */
      xprintf("Time used:   %.1f secs\n", xdifftime(xtime(), start));
      {  glp_long tpeak;
         char buf[50];
         glp_mem_usage(NULL, NULL, NULL, &tpeak);
         xprintf("Memory used: %.1f Mb (%s bytes)\n",
            xltod(tpeak) / 1048576.0, xltoa(tpeak, buf));
      }
      /*--------------------------------------------------------------*/
skip: /* postsolve the model, if necessary */
      if (csa->tran != NULL)
      {  if (csa->solution == SOL_BASIC)
            ret = glp_mpl_postsolve(csa->tran, csa->prob, GLP_SOL);
         else if (csa->solution == SOL_INTERIOR)
            ret = glp_mpl_postsolve(csa->tran, csa->prob, GLP_IPT);
         else if (csa->solution == SOL_INTEGER)
            ret = glp_mpl_postsolve(csa->tran, csa->prob, GLP_MIP);
         else
            xassert(csa != csa);
         if (ret != 0)
         {  xprintf("Model postsolving error\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /*--------------------------------------------------------------*/
      /* write problem solution in printable format, if required */
      if (csa->out_sol != NULL)
      {  if (csa->solution == SOL_BASIC)
            ret = lpx_print_sol(csa->prob, csa->out_sol);
         else if (csa->solution == SOL_INTERIOR)
            ret = lpx_print_ips(csa->prob, csa->out_sol);
         else if (csa->solution == SOL_INTEGER)
            ret = lpx_print_mip(csa->prob, csa->out_sol);
         else
            xassert(csa != csa);
         if (ret != 0)
         {  xprintf("Unable to write problem solution\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write problem solution in printable format, if required */
      if (csa->out_res != NULL)
      {  if (csa->solution == SOL_BASIC)
            ret = glp_write_sol(csa->prob, csa->out_res);
         else if (csa->solution == SOL_INTERIOR)
            ret = glp_write_ipt(csa->prob, csa->out_res);
         else if (csa->solution == SOL_INTEGER)
            ret = glp_write_mip(csa->prob, csa->out_res);
         else
            xassert(csa != csa);
         if (ret != 0)
         {  xprintf("Unable to write problem solution\n");
            ret = EXIT_FAILURE;
            goto done;
         }
      }
      /* write sensitivity analysis report, if required */
      if (csa->out_ranges != NULL)
      {  if (csa->solution == SOL_BASIC)
         {  if (glp_get_status(csa->prob) == GLP_OPT)
            {  if (glp_bf_exists(csa->prob))
ranges:        {  ret = glp_print_ranges(csa->prob, 0, NULL, 0,
                     csa->out_ranges);
                  if (ret != 0)
                  {  xprintf("Unable to write sensitivity analysis repo"
                        "rt\n");
                     ret = EXIT_FAILURE;
                     goto done;
                  }
               }
               else
               {  ret = glp_factorize(csa->prob);
                  if (ret == 0) goto ranges;
                  xprintf("Cannot produce sensitivity analysis report d"
                     "ue to error in basis factorization (glp_factorize"
                     " returned %d); try --nopresol\n", ret);
               }
            }
            else
               xprintf("Cannot produce sensitivity analysis report for "
                  "non-optimal basic solution\n");
         }
         else
            xprintf("Cannot produce sensitivity analysis report for int"
               "erior-point or MIP solution\n");
      }
      /*--------------------------------------------------------------*/
      /* all seems to be ok */
      ret = EXIT_SUCCESS;
      /*--------------------------------------------------------------*/
done: /* delete the LP/MIP problem object */
      if (csa->prob != NULL)
         glp_delete_prob(csa->prob);
      /* free the translator workspace, if necessary */
      if (csa->tran != NULL)
         glp_mpl_free_wksp(csa->tran);
      /* delete the network problem object, if necessary */
      if (csa->graph != NULL)
         glp_delete_graph(csa->graph);
      xassert(gmp_pool_count() == 0);
      gmp_free_mem();
      /* close log file, if necessary */
      if (csa->log_file != NULL) glp_close_tee();
      /* check that no memory blocks are still allocated */
      {  int count;
         glp_long total;
         glp_mem_usage(&count, NULL, &total, NULL);
         if (count != 0)
            xerror("Error: %d memory block(s) were lost\n", count);
         xassert(count == 0);
         xassert(total.lo == 0 && total.hi == 0);
      }
      /* free the GLPK environment */
      glp_free_env();
      /* return to the control program */
      return ret;
}

/* eof */
