"""
Configuration framework for igraph.

igraph has some parameters which usually affect the behaviour of many functions.
This module provides the framework for altering and querying igraph parameters
as well as saving them to and retrieving them from disk.
"""

__license__ = """
Copyright (C) 2006-2007  Gabor Csardi <csardi@rmki.kfki.hu>,
Tamas Nepusz <ntamas@rmki.kfki.hu>

MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
02110-1301 USA
"""
from ConfigParser import SafeConfigParser
import platform
import os.path

def get_platform_image_viewer():
    """Returns the path of an image viewer on the given platform"""
    plat = platform.system()
    if plat == "Darwin":
        # Most likely Mac OS X
        return "open"
    elif plat == "Linux":
        # Linux has a whole lot of choices, try to find one
        choices = ["gthumb", "gqview", "kuickshow", "xnview", "display"]
        paths = ["/usr/bin", "/bin"]
        for path in paths:
            for choice in choices:
                full_path = os.path.join(path, choice)
                if os.path.isfile(full_path): return full_path
        return ""
    elif plat == "Windows" or plat == "Microsoft":    # Thanks to Dale Hunscher
        # Use the built-in Windows image viewer, if available
        return "start"
    else:
        # Unknown system
        return ""

class Configuration(object):
    """Class representing igraph configuration details.
    
    General things
    ==============

    The configuration of igraph is stored in the form of name-value pairs.
    This object provides an interface to the configuration data using the
    syntax known from dict:
    
      >>> c=Configuration()
      >>> c["general.verbose"] = True
      >>> print c["general.verbose"]
      True

    Configuration keys are organized into sections, and the name to be used
    for a given key is always in the form C{section.keyname}, like
    C{general.verbose} in the example above. In that case, C{general} is the
    name of the configuration section, and C{verbose} is the name of the key.
    If the name of the section is omitted, it defaults to C{general}, so
    C{general.verbose} can be referred to as C{verbose}:

      >>> c=Configuration()
      >>> c["verbose"] = True
      >>> print c["general.verbose"]
      True

    User-level configuration is stored in C{~/.igraphrc} per default on Linux
    and Mac OS X systems, or in C{C:\Documents and Settings\username\.igraphrc}
    on Windows systems. However, this configuration is read only when C{igraph}
    is launched through its shell interface defined in L{igraph.app.shell}.
    This behaviour might change before version 1.0.

    Known configuration keys
    ========================

    The known configuration keys are presented below, sorted by section. When
    referring to them in program code, don't forget to add the section name,
    expect in the case of section C{general}.

    General settings
    ----------------

    These settings are all stored in section C{general}.

        - B{shells}: the list of preferred Python shells to be used with the
          command-line C{igraph} script. The shells in the list are tried one
          by one until any of them is found on the system. C{igraph} functions
          are then imported into the main namespace of the shell and the shell
          is launched. Known shells and their respective class names to be
          used can be found in L{igraph.app.shell}. Example:
          C{IPythonShell, ClassicPythonShell}. This is the default, by the way.
        - B{verbose}: whether L{igraph} should talk more than really necessary.
          For instance, if set to C{True}, some functions display progress bars.

    Application settings
    --------------------

    These settings specify the external applications that are possibly
    used by C{igraph}. They are all stored in section C{apps}.

        - B{image_viewer}: image viewer application. If set to an empty string,
          it will be determined automatically from the platform C{igraph} runs
          on. On Mac OS X, it defaults to the Preview application. On Linux,
          it chooses a viewer from several well-known Linux viewers like
          C{gthumb}, C{kuickview} and so on (see the source code for the full
          list). On Windows, it defaults to the system's built-in image viewer.

    Plotting settings
    -----------------

    These settings specify the default values used by plotting functions.
    They are all stored in section C{plotting}.

        - B{layout}: default graph layout algorithm to be used.
        - B{palette}: default palette to be used for converting integer
          numbers to colors. See L{colors.Palette} for more information.
          Valid palette names are stored in C{colors.palettes}.
    """
    class Types:
        """Static class for the implementation of custom getter/setter functions
        for configuration keys"""
        def setboolean(object, section, key, value):
            if str(value).lower() in ["0", "false", "no", "off"]:
                value="false"
            elif str(value).lower() in ["1", "true", "yes", "on"]:
                value="true"
            else:
                raise ValueError, "value cannot be coerced to boolean type"
            object.set(section, key, value)
        setboolean=staticmethod(setboolean)

    _types = {
        "boolean": {
            "getter": SafeConfigParser.getboolean,
            "setter": Types.setboolean
        }
    }

    _sections = ("general", "apps", "plotting", )
    _definitions = {
        "general.shells": { "default": "IPythonShell,ClassicPythonShell" },
        "general.verbose": { "default": True, "type": "boolean" },

        "apps.image_viewer": { "default": get_platform_image_viewer() },

        "plotting.layout": { "default": "random" },
        "plotting.palette": { "default": "gray" }
    }

    def __init__(self, filename=None):
        """Creates a new configuration instance.

        @param filename: file or file pointer to be read. Can be omitted.
        """
        self._config = SafeConfigParser()
        self._filename = None

        # Create default sections
        for sec in self._sections: self._config.add_section(sec)
        # Create default values
        for name, definition in self._definitions.iteritems():
            if definition.has_key("default"):
                self[name]=definition["default"]

        if filename is not None: self.load(filename)

    def _get_filename(self):
        """Returns the filename associated to the object.

        It is usually the name of the configuration file that was used when
        creating the object. L{Configuration.load} always overwrites it with
        the filename given to it. If C{None}, the configuration was either
        created from scratch or it was updated from a stream without name
        information."""
        return self._filename
    filename=property(_get_filename, doc=_get_filename.__doc__)

    def _item_to_section_key(item):
        """Converts an item description to a section-key pair.
        
        @param item: the item to be converted
        @return: if C{item} contains a period (C{.}), it is splitted into two parts
          at the first period, then the two parts are returned, so the part before
          the period is the section. If C{item} does not contain a period, the
          section is assumed to be C{general}, and the second part of the returned
          pair contains C{item} unchanged"""
        if "." in item:
            section, key = item.split(".", 1)
        else:
            section, key = "general", item
        return section, key
    _item_to_section_key=staticmethod(_item_to_section_key)

    def __getitem__(self, item):
        """Returns the given configuration item.

        @param item: the configuration key to retrieve.
        @return: the configuration value"""
        section, key = self._item_to_section_key(item)
        definition = self._definitions.get("%s.%s" % (section,key), {})
        getter = None
        if definition.has_key("type"):
            getter = self._types[definition["type"]].get("getter", None)
        if getter is None: getter = self._config.__class__.get
        return getter(self._config, section, key)

    def __setitem__(self, item, value):
        """Sets the given configuration item.

        @param item: the configuration key to set
        @param value: the new value of the configuration key
        """
        section, key = self._item_to_section_key(item)
        definition = self._definitions.get("%s.%s" % (section,key), {})
        setter = None
        if definition.has_key("type"):
            setter = self._types[definition["type"]].get("setter", None)
        if setter is None: setter = self._config.__class__.set
        return setter(self._config, section, key, value)

    def __delitem__(self, item):
        """Deletes the given item from the configuration.

        If the item has a default value, the default value is written back instead
        of the current value. Without a default value, the item is really deleted.
        """
        section, key = self._item_to_section_key(item)
        definition = self._definitions.get("%s.%s" % (section, key), {})
        if definition.has_key("default"):
            self[item]=definition["default"]
        else:
            self._config.remove_option(section, key)

    def has_key(self, item):
        """Checks if the configuration has a given key.

        @param item: the key being sought"""
        if "." in item:
            section, key = item.split(".", 1)
        else:
            section, key = "general", item
        return self._config.has_option(section, key)
            
    def load(self, fp=None):
        """Loads the configuration from the given file.

        @param fp: name of a file or a file object. The configuration will be loaded
          from here. Can be omitted, in this case, the user-level configuration is
          loaded.
        """
        fp = fp or get_user_config_file()
        if not isinstance(fp, file):
            fp=open(fp, "r")
            file_was_open=True
        self._config.readfp(fp)
        self._filename = getattr(fp, "name", None)
        if file_was_open: fp.close()

    def save(self, fp=None):
        """Saves the configuration.

        @param fp: name of a file or a file object. The configuration will be saved
          there. Can be omitted, in this case, the user-level configuration file will
          be overwritten.
        """
        fp = fp or get_user_config_file()
        if not isinstance(fp, file):
            fp=open(fp, "w")
            file_was_open=True
        self._config.write(fp)
        if file_was_open: fp.close()

def get_user_config_file():
    """Returns the path where the user-level configuration file is stored"""
    import os.path
    return os.path.expanduser("~/.igraphrc")


def init():
    """Default mechanism to initiate igraph configuration

    This method loads the user-specific configuration file from the
    user's home directory, or if it does not exist, creates a default
    configuration
    
    @return: the L{Configuration} object loaded or created"""
    cfile = get_user_config_file()
    try:
        config = Configuration(cfile)
        return config
    except IOError:
        # No config file yet, whatever
        config = Configuration()
        return config


