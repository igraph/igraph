"""
This package contains third party libraries that igraph depends on and that
are small enough to be distributed with igraph itself.

The primary entry point of this module is ``vendor_import``, a function that
first tries to import a particular library using the standard Python mechanism
and falls back to the version of the library provided within ``igraph.vendor``
if the standard Python import fails.

The libraries contained within igraph are as follows:

	- `texttable`, a library to print ASCII tables, by Gerome Fournier.
	  See <http://foutaise.org/code/>.
"""

__license__ = "GPL"

__all__ = ["vendor_import"]
__docformat__ = "restructuredtext en"

def vendor_import(module_name):
	"""Tries to import a module name ``module_name`` using the standard Python
	`import` statement and return the imported module. If the import fails,
	tries to import a module of the same name from within ``igraph.vendor``
	and return that module instead.
	"""

	parts = module_name.split(".")

	try:
		result = __import__(module_name, level=0)
	except ImportError:
		result = __import__("igraph.vendor.%s" % module_name, level=0)
		parts[0:0] = ["igraph", "vendor"]

	parts.pop(0)
	while parts:
		result = getattr(result, parts.pop(0))

	return result
