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

    _sections = ("general", )
    _definitions = {
        "general.shells": { "default": "IPythonShell,ClassicPythonShell" },
        "general.verbose": { "default": True, "type": "boolean" }
    }

    def __init__(self, filename=None):
        """Creates a new configuration instance.

        @param fp: file or file pointer to be read. Can be omitted.
        """
        self._config = SafeConfigParser()

        # Create default sections
        for sec in self._sections: self._config.add_section(sec)
        # Create default values
        for name, definition in self._definitions.iteritems():
            if definition.has_key("default"):
                self[name]=definition["default"]

        if filename is not None: self.load(filename)

    def __getitem__(self, item):
        """Returns the given configuration item.

        @param item: the configuration key to retrieve.
        @return: the configuration value"""
        if "." in item:
            section, key = item.split(".", 1)
        else:
            section, key = "general", item
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
        if "." in item:
            section, key = item.split(".", 1)
        else:
            section, key = "general", item
        definition = self._definitions.get("%s.%s" % (section,key), {})
        setter = None
        if definition.has_key("type"):
            setter = self._types[definition["type"]].get("setter", None)
        if setter is None: setter = self._config.__class__.set
        return setter(self._config, section, key, value)

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

