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

# pylint: disable-msg=W0401
# W0401: wildcard import. That's exactly what we need for the shell.
from igraph import __version__, set_progress_handler, set_status_handler
from igraph.configuration import Configuration
import sys, re

# pylint: disable-msg=C0103,R0903
# C0103: invalid name. Disabled because this is a third-party class.
# R0903: too few public methods.
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
        This is green

    Alternatively, the `render()` method can used, which replaces
    '${action}' with the string required to perform 'action':

        >>> term = TerminalController()
        >>> print term.render('This is ${GREEN}green${NORMAL}')
        This is green

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
        ...

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
        try:
            import curses
        except ImportError:
            return

        # If the stream isn't a tty, then assume it has no capabilities.
        if not term_stream.isatty():
            return

        # Check the terminal type.  If we fail, then assume that the
        # terminal has no capabilities.
        try:
            curses.setupterm()
        except StandardError:
            return

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
            for i, color in zip(range(len(self._COLORS)), self._COLORS):
                setattr(self, color, curses.tparm(set_fg, i) or '')
        set_fg_ansi = self._tigetstr('setaf')
        if set_fg_ansi:
            for i, color in zip(range(len(self._ANSICOLORS)), self._ANSICOLORS):
                setattr(self, color, curses.tparm(set_fg_ansi, i) or '')
        set_bg = self._tigetstr('setb')
        if set_bg:
            for i, color in zip(range(len(self._COLORS)), self._COLORS):
                setattr(self, 'BG_'+color, curses.tparm(set_bg, i) or '')
        set_bg_ansi = self._tigetstr('setab')
        if set_bg_ansi:
            for i, color in zip(range(len(self._ANSICOLORS)), self._ANSICOLORS):
                setattr(self, 'BG_'+color, curses.tparm(set_bg_ansi, i) or '')

    @staticmethod
    def _tigetstr(cap_name):
        """Rewrites string capabilities to remove "delays" which are not
        required for modern terminals"""
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
        """Helper function for L{render}"""
        s = match.group()
        if s == '$$':
            return s
        else:
            return getattr(self, s[2:-1])


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
        self.progress_bar = term.render(self.BAR)
        self.header = self.term.render(self.HEADER % "".center(self.width))
        self.cleared = True #: true if we haven't drawn the bar yet.

        self.last_percent = 0
        self.last_message = ""

    def update(self, percent=None, message=None):
        """Updates the progress bar.

        @param percent: the percentage to be shown. If C{None}, the previous
          value will be used.
        @param message: the message to be shown above the progress bar. If
          C{None}, the previous message will be used.
        """
        if self.cleared:
            sys.stdout.write("\n"+self.header)
            self.cleared = False

        if message is None:
            message = self.last_message
        else:
            self.last_message = message

        if percent is None:
            percent = self.last_percent
        else:
            self.last_percent = percent

        n = int((self.width-10)*(percent/100.0))
        sys.stdout.write(
            self.term.BOL + self.term.UP + self.term.UP + self.term.CLEAR_EOL +
            self.term.render(self.HEADER % message.center(self.width)) +
            (self.progress_bar % (percent, '='*n, '-'*(self.width-10-n))) + "\n"
            )

    def update_message(self, message):
        """Updates the message of the progress bar.

        @param message: the message to be shown above the progress bar
        """
        return self.update(message=message.strip())

    def clear(self):
        """Clears the progress bar (i.e. removes it from the screen)"""
        if not self.cleared:
            sys.stdout.write(self.term.BOL + self.term.CLEAR_EOL +
                             self.term.UP + self.term.CLEAR_EOL +
                             self.term.UP + self.term.CLEAR_EOL)
            self.cleared = True
            self.last_percent = 0
            self.last_message = ""


class Shell(object):
    """Superclass of the embeddable shells supported by igraph"""

    def __init__(self):
        pass

    def __call__(self):
        raise NotImplementedError("abstract class")

    def supports_progress_bar(self):
        """Checks whether the shell supports progress bars.

        This is done by checking for the existence of an attribute
        called C{_progress_handler}."""
        return hasattr(self, "_progress_handler")

    def supports_status_messages(self):
        """Checks whether the shell supports status messages.

        This is done by checking for the existence of an attribute
        called C{_status_handler}."""
        return hasattr(self, "_status_handler")

    # pylint: disable-msg=E1101
    def get_progress_handler(self):
        """Returns the progress handler (if exists) or None (if not)."""
        if self.supports_progress_bar():
            return self._progress_handler
        return None

    # pylint: disable-msg=E1101
    def get_status_handler(self):
        """Returns the status handler (if exists) or None (if not)."""
        if self.supports_status_messages():
            return self._status_handler
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
        Shell.__init__(self)

        import idlelib.PyShell
        
        idlelib.PyShell.use_subprocess = True
        
        try:
            sys.ps1
        except AttributeError:
            sys.ps1 = '>>> '

        root = idlelib.PyShell.Tk(className="Idle")
        idlelib.PyShell.fixwordbreaks(root)
        root.withdraw()
        flist = idlelib.PyShell.PyShellFileList(root)
        if not flist.open_shell():
            raise NotImplementedError
        self._shell = flist.pyshell
        self._root = root

    def __call__(self):
        """Starts the shell"""
        self._shell.interp.execsource("from igraph import *")
        self._root.mainloop()
        self._root.destroy()


class ConsoleProgressBarMixin(object):
    """Mixin class for console shells that support a progress bar."""

    def __init__(self):
        try:
            self.__class__.progress_bar = ProgressBar(TerminalController())
        except ValueError:
            # Terminal is not capable enough, disable progress handler
            del self.__class__._progress_handler
            del self.__class__._status_handler

    @classmethod
    def _progress_handler(cls, message, percentage):
        """Progress bar handler, called when C{igraph} reports the progress
        of an operation

        @param message: message provided by C{igraph}
        @param percentage: percentage provided by C{igraph}
        """
        if percentage >= 100:
            cls.progress_bar.clear()
        else:
            cls.progress_bar.update(percentage, message)

    @classmethod
    def _status_handler(cls, message):
        """Status message handler, called when C{igraph} sends a status
        message to be displayed.

        @param message: message provided by C{igraph}
        """
        cls.progress_bar.update_message(message)


class IPythonShell(Shell, ConsoleProgressBarMixin):
    """IPython embedded shell interface.

    This class allows igraph to be embedded in IPython's interactive shell."""

    def __init__(self):
        """Constructor.

        Imports IPython's embedded shell with separator lines removed."""
        Shell.__init__(self)
        ConsoleProgressBarMixin.__init__(self)

        # We cannot use IPShellEmbed here because generator expressions do not
        # work there (e.g., set(g.degree(x) for x in [1,2,3])) where g comes
        # from an external context
        import sys

        from IPython import __version__ as ipython_version
        self.ipython_version = ipython_version

        try:
            # IPython >= 0.11 supports this
            from IPython.frontend.terminal.ipapp import TerminalIPythonApp
            self._shell = TerminalIPythonApp.instance()
            sys.argv.append("--nosep")
        except ImportError:
            # IPython 0.10 and earlier
            import IPython.Shell
            self._shell = IPython.Shell.start()
            self._shell.IP.runsource("from igraph import *")
            sys.argv.append("-nosep")

    def __call__(self):
        """Starts the embedded shell."""
        print "igraph %s running inside " % __version__,
        if self._shell.__class__.__name__ == "TerminalIPythonApp":
            self._shell.initialize()
            self._shell.shell.ex("from igraph import *")
            self._shell.start()
        else:
            self._shell.mainloop()


class ClassicPythonShell(Shell, ConsoleProgressBarMixin):
    """Classic Python shell interface.

    This class allows igraph to be embedded in Python's shell."""
    
    def __init__(self):
        """Constructor.

        Imports Python's classic shell"""
        Shell.__init__(self)
        self._shell = None

        try:
            self.__class__.progress_bar = ProgressBar(TerminalController())
        except ValueError:
            # Terminal is not capable enough, disable progress handler
            del self.__class__._progress_handler

    def __call__(self):
        """Starts the embedded shell."""
        if self._shell is None:
            from code import InteractiveConsole
            self._shell = InteractiveConsole()
            print >> sys.stderr, "igraph %s running inside " % __version__,
            self._shell.runsource("from igraph import *")

        self._shell.interact()


def main():
    """The main entry point for igraph when invoked from the command
    line shell"""
    config = Configuration.instance()

    if config.filename:
        print >> sys.stderr, "Using configuration from %s" % config.filename
    else:
        print >> sys.stderr, "No configuration file, using defaults"

    if config.has_key("shells"):
        parts = [part.strip() for part in config["shells"].split(",")]
        shell_classes = []
        available_classes = dict([(k, v) for k, v in globals().iteritems() \
            if isinstance(v, type) and issubclass(v, Shell)])
        for part in parts:
            klass = available_classes.get(part, None)
            if klass is None:
                print >> sys.stderr, "Warning: unknown shell class `%s'" % part
                continue
            shell_classes.append(klass)
    else:
        shell_classes = [IPythonShell, ClassicPythonShell]
        import platform
        if platform.system() == "Windows":
            shell_classes.insert(0, IDLEShell)

    shell = None
    for shell_class in shell_classes:
        # pylint: disable-msg=W0703
        # W0703: catch "Exception"
        try:
            shell = shell_class()
            break
        except StandardError:
            # Try the next one
            pass

    if isinstance(shell, Shell):
        if config["verbose"]:
            if shell.supports_progress_bar():
                set_progress_handler(shell.get_progress_handler())
            if shell.supports_status_messages():
                set_status_handler(shell.get_status_handler())
        shell()
    else:
        print >> sys.stderr, "No suitable Python shell was found."
        print >> sys.stderr, "Check configuration variable `general.shells'."

if __name__ == '__main__':
    sys.exit(main())

