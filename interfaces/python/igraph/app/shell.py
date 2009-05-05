"""Command-line user interface of igraph

The command-line interface launches a Python shell with the igraph
module automatically imported into the main namespace. This is mostly a
convenience module and it is used only by the C{igraph} command line
script which executes a suitable Python shell and automatically imports
C{igraph}'s classes and functions in the top-level namespace.

Supported Python shells are:

  - IDLE shell (class L{IDLEShell})
  - IPython shell (class L{IPythonShell})
  - Classic Python shell (class L{ClassicPythonShell})

The shells are tried in the above mentioned preference order one by
one, unless the C{global.shells} configuration key is set which
overrides the default order. IDLE shell is only tried in Windows
unless explicitly stated by C{global.shells}, since Linux and
Mac OS X users are likely to invoke igraph from the command line.
"""
from igraph import *
from igraph import __version__
import sys, re

class TerminalController:
    """
    A class that can be used to portably generate formatted output to
    a terminal.  
    
    `TerminalController` defines a set of instance variables whose
    values are initialized to the control sequence necessary to
    perform a given action.  These can be simply included in normal
    output to the terminal:

        >>> term = TerminalController()
        >>> print 'This is '+term.GREEN+'green'+term.NORMAL

    Alternatively, the `render()` method can used, which replaces
    '${action}' with the string required to perform 'action':

        >>> term = TerminalController()
        >>> print term.render('This is ${GREEN}green${NORMAL}')

    If the terminal doesn't support a given action, then the value of
    the corresponding instance variable will be set to ''.  As a
    result, the above code will still work on terminals that do not
    support color, except that their output will not be colored.
    Also, this means that you can test whether the terminal supports a
    given action by simply testing the truth value of the
    corresponding instance variable:

        >>> term = TerminalController()
        >>> if term.CLEAR_SCREEN:
        ...     print 'This terminal supports clearning the screen.'

    Finally, if the width and height of the terminal are known, then
    they will be stored in the `COLS` and `LINES` attributes.
    
    @author: Edward Loper
    """
    # Cursor movement:
    BOL = ''             #: Move the cursor to the beginning of the line
    UP = ''              #: Move the cursor up one line
    DOWN = ''            #: Move the cursor down one line
    LEFT = ''            #: Move the cursor left one char
    RIGHT = ''           #: Move the cursor right one char

    # Deletion:
    CLEAR_SCREEN = ''    #: Clear the screen and move to home position
    CLEAR_EOL = ''       #: Clear to the end of the line.
    CLEAR_BOL = ''       #: Clear to the beginning of the line.
    CLEAR_EOS = ''       #: Clear to the end of the screen

    # Output modes:
    BOLD = ''            #: Turn on bold mode
    BLINK = ''           #: Turn on blink mode
    DIM = ''             #: Turn on half-bright mode
    REVERSE = ''         #: Turn on reverse-video mode
    NORMAL = ''          #: Turn off all modes

    # Cursor display:
    HIDE_CURSOR = ''     #: Make the cursor invisible
    SHOW_CURSOR = ''     #: Make the cursor visible

    # Terminal size:
    COLS = None          #: Width of the terminal (None for unknown)
    LINES = None         #: Height of the terminal (None for unknown)

    # Foreground colors:
    BLACK = BLUE = GREEN = CYAN = RED = MAGENTA = YELLOW = WHITE = ''
    
    # Background colors:
    BG_BLACK = BG_BLUE = BG_GREEN = BG_CYAN = ''
    BG_RED = BG_MAGENTA = BG_YELLOW = BG_WHITE = ''
    
    _STRING_CAPABILITIES = """
    BOL=cr UP=cuu1 DOWN=cud1 LEFT=cub1 RIGHT=cuf1
    CLEAR_SCREEN=clear CLEAR_EOL=el CLEAR_BOL=el1 CLEAR_EOS=ed BOLD=bold
    BLINK=blink DIM=dim REVERSE=rev UNDERLINE=smul NORMAL=sgr0
    HIDE_CURSOR=cinvis SHOW_CURSOR=cnorm""".split()
    _COLORS = """BLACK BLUE GREEN CYAN RED MAGENTA YELLOW WHITE""".split()
    _ANSICOLORS = "BLACK RED GREEN YELLOW BLUE MAGENTA CYAN WHITE".split()

    def __init__(self, term_stream=sys.stdout):
        """
        Create a `TerminalController` and initialize its attributes
        with appropriate values for the current terminal.
        `term_stream` is the stream that will be used for terminal
        output; if this stream is not a tty, then the terminal is
        assumed to be a dumb terminal (i.e., have no capabilities).
        """
        # Curses isn't available on all platforms
        try: import curses
        except: return

        # If the stream isn't a tty, then assume it has no capabilities.
        if not term_stream.isatty(): return

        # Check the terminal type.  If we fail, then assume that the
        # terminal has no capabilities.
        try: curses.setupterm()
        except: return

        # Look up numeric capabilities.
        self.COLS = curses.tigetnum('cols')
        self.LINES = curses.tigetnum('lines')
        
        # Look up string capabilities.
        for capability in self._STRING_CAPABILITIES:
            (attrib, cap_name) = capability.split('=')
            setattr(self, attrib, self._tigetstr(cap_name) or '')

        # Colors
        set_fg = self._tigetstr('setf')
        if set_fg:
            for i,color in zip(range(len(self._COLORS)), self._COLORS):
                setattr(self, color, curses.tparm(set_fg, i) or '')
        set_fg_ansi = self._tigetstr('setaf')
        if set_fg_ansi:
            for i,color in zip(range(len(self._ANSICOLORS)), self._ANSICOLORS):
                setattr(self, color, curses.tparm(set_fg_ansi, i) or '')
        set_bg = self._tigetstr('setb')
        if set_bg:
            for i,color in zip(range(len(self._COLORS)), self._COLORS):
                setattr(self, 'BG_'+color, curses.tparm(set_bg, i) or '')
        set_bg_ansi = self._tigetstr('setab')
        if set_bg_ansi:
            for i,color in zip(range(len(self._ANSICOLORS)), self._ANSICOLORS):
                setattr(self, 'BG_'+color, curses.tparm(set_bg_ansi, i) or '')

    def _tigetstr(self, cap_name):
        # String capabilities can include "delays" of the form "$<2>".
        # For any modern terminal, we should be able to just ignore
        # these, so strip them out.
        import curses
        cap = curses.tigetstr(cap_name) or ''
        return re.sub(r'\$<\d+>[/*]?', '', cap)

    def render(self, template):
        """
        Replace each $-substitutions in the given template string with
        the corresponding terminal control string (if it's defined) or
        '' (if it's not).
        """
        return re.sub(r'\$\$|\${\w+}', self._render_sub, template)

    def _render_sub(self, match):
        s = match.group()
        if s == '$$': return s
        else: return getattr(self, s[2:-1])


class ProgressBar:
    """
    A 2-line progress bar, which looks like::
    
                                Header
        20% [===========----------------------------------]

    The progress bar is colored, if the terminal supports color
    output; and adjusts to the width of the terminal.
    """
    BAR = '%3d%% ${GREEN}[${BOLD}%s%s${NORMAL}${GREEN}]${NORMAL}'
    HEADER = '${BOLD}${CYAN}%s${NORMAL}\n'
        
    def __init__(self, term):
        self.term = term
        if not (self.term.CLEAR_EOL and self.term.UP and self.term.BOL):
            raise ValueError("Terminal isn't capable enough -- you "
                             "should use a simpler progress display.")
        self.width = self.term.COLS or 75
        self.bar = term.render(self.BAR)
        self.header = self.term.render(self.HEADER % "".center(self.width))
        self.cleared = 1 #: true if we haven't drawn the bar yet.

    def update(self, percent, message):
        if self.cleared:
            sys.stdout.write(self.header)
            self.cleared = 0
        n = int((self.width-10)*(percent/100.0))
        sys.stdout.write(
            self.term.BOL + self.term.UP + self.term.CLEAR_EOL +
            self.term.render(self.HEADER % message.center(self.width)) +
            (self.bar % (percent, '='*n, '-'*(self.width-10-n)))
            )

    def clear(self):
        if not self.cleared:
            sys.stdout.write(self.term.BOL + self.term.CLEAR_EOL +
                             self.term.UP + self.term.CLEAR_EOL)
            self.cleared = 1

class Shell(object):
    """Superclass of the embeddable shells supported by igraph"""

    def __init__(self):
        raise ValueError, "abstract class"
    def __call__(self, namespace=None):
        raise ValueError, "abstract class"
    def supports_progress_bar(self):
        return hasattr(self, "_progress_handler")
    def get_progress_handler(self):
        if self.supports_progress_bar(): return self._progress_handler
        return None

class IDLEShell(Shell):
    """IDLE embedded shell interface.

    This class allows igraph to be embedded in IDLE (the Tk Python IDE).
    
    @todo: no progress bar support yet. Shell/Restart Shell command should
      re-import igraph again."""
    
    def __init__(self):
        """Constructor.

        Imports IDLE's embedded shell. The implementation of this method is
        ripped from idlelib.PyShell.main() after removing the unnecessary
        parts."""
        import idlelib.PyShell
        import sys
        
        idlelib.PyShell.use_subprocess = True
        
        try:
            sys.ps1
        except AttributeError:
            sys.ps1 = '>>> '

        root = idlelib.PyShell.Tk(className="Idle")
        idlelib.PyShell.fixwordbreaks(root)
        root.withdraw()
        flist = idlelib.PyShell.PyShellFileList(root)
        if not flist.open_shell(): raise NotImplementedError
        self._shell = flist.pyshell
        self._root = root

    def __call__(self, namespace=None):
        """Starts the shell"""
        self._shell.interp.execsource("from igraph import *")
        self._root.mainloop()
        self._root.destroy()


class IPythonShell(Shell):
    """IPython embedded shell interface.

    This class allows igraph to be embedded in IPython's interactive shell."""

    def __init__(self):
        """Constructor.

        Imports IPython's embedded shell with separator lines removed."""
        from IPython.Shell import IPShellEmbed
        self._shell = IPShellEmbed(['-nosep'])
        try:
            self.__class__.progress_bar = ProgressBar(TerminalController())
        except ValueError:
            # Terminal is not capable enough, disable progress handler
            del self.__class__._progress_handler

    def __call__(self, namespace=None):
        """Starts the embedded shell.
        
        @param namespace: global namespace to use"""
        print "igraph %s running inside" % __version__,
        print self._shell.IP.BANNER,
        self._shell(local_ns=namespace)

    def _progress_handler(message, percentage):
        """Progress bar handler, called when C{igraph} reports the progress
        of an operation

        @param message: message provided by C{igraph}
        @param percentage: percentage provided by C{igraph}
        """
        if percentage >= 100:
            IPythonShell.progress_bar.clear()
        else:
            IPythonShell.progress_bar.update(percentage, message)
    _progress_handler = staticmethod(_progress_handler)

class ClassicPythonShell(Shell):
    """Classic Python shell interface.

    This class allows igraph to be embedded in Python's shell."""
    
    def __init__(self):
        """Constructor.

        Imports Python's classic shell"""
        from code import InteractiveConsole
        try:
            self.__class__.progress_bar = ProgressBar(TerminalController())
        except ValueError:
            # Terminal is not capable enough, disable progress handler
            del self.__class__._progress_handler

    def __call__(self, namespace=None):
        """Starts the embedded shell.
        
        @param namespace: global namespace to use"""
        from code import InteractiveConsole
        self._shell = InteractiveConsole(locals=namespace)
        print >>sys.stderr, "igraph %s running inside " % __version__,
        self._shell.runsource("from igraph import *")
        self._shell.interact()

    def _progress_handler(message, percentage):
        """Progress bar handler, called when C{igraph} reports the progress
        of an operation

        @param message: message provided by C{igraph}
        @param percentage: percentage provided by C{igraph}
        """
        if percentage >= 100:
            ClassicPythonShell.progress_bar.clear()
        else:
            ClassicPythonShell.progress_bar.update(percentage, message)
    _progress_handler = staticmethod(_progress_handler)

def main():
    if config.filename:
        print >>sys.stderr, "Using configuration from %s" % config.filename
    else:
        print >>sys.stderr, "No configuration file, using defaults"

    if config.has_key("shells"):
        parts = [part.strip() for part in config["shells"].split(",")]
        shell_classes = []
        available_classes = dict([(k, v) for k, v in globals().iteritems() \
            if isinstance(v, type) and issubclass(v, Shell)])
        for part in parts:
            klass = available_classes.get(part, None)
            if klass is None:
                print >>sys.stderr, "Warning: unknown shell class `%s'" % part
            else:
                shell_classes.append(klass)
    else:
        shell_classes = [IPythonShell, ClassicPythonShell]
        import platform
        if platform.system() == "Windows":
            shell_classes.insert(0, IDLEShell)

    shell = None
    for shell_class in shell_classes:
        try:
            shell = shell_class()
            break
        except:
            # Try the next one
            pass

    if isinstance(shell, Shell):
        if config["verbose"] and shell.supports_progress_bar():
            set_progress_handler(shell.get_progress_handler())
        shell()
    else:
        print >>sys.stderr, "No suitable Python shell was found."
        print >>sys.stderr, "Check configuration variable `general.shells'."

if __name__ == '__main__': main()

