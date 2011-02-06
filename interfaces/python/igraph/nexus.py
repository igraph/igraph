"""Interface to the Nexus online graph repository.

The classes in this file facilitate access to the Nexus online graph
repository at U{http://nexus.igraph.org}."""

__all__ = ["Nexus"]

from textwrap import TextWrapper
from urllib import urlencode
from urllib2 import urlopen
from urlparse import urlparse, urlunparse

from igraph.compat import property
from igraph.configuration import Configuration
from igraph.utils import multidict

class NexusConnection(object):
    """Connection to a remote Nexus server."""

    def __init__(self, nexus_url=None):
        """Constructs a connection to a remote Nexus server.

        @param nexus_url: the root URL of the remote server. Leave it at its
          default value (C{None}) unless you have set up your own Nexus server
          and you want to connect to that. C{None} fetches the URL from
          igraph's configuration file or uses the default URL if no URL
          is specified in the configuration file.
        """
        self.url = Configuration.instance()["remote.nexus.url"]

    def get(self, id):
        """Retrieves the dataset with the given ID from Nexus.

        @param id: the ID of the dataset to retrieve.
        @return: an instance of L{Graph}.
        """
        from igraph import load

        params = dict(format="Python-igraph", id=id)
        response = self._get_response("/api/dataset", params)
        return load(response, format="pickle")

    def info(self, id):
        """Retrieves informations about the dataset with the given numeric
        or string ID from Nexus.

        @param id: the numeric or string ID of the dataset to retrieve.
        @return: an instance of L{NexusDatasetInfo}.
        """
        params = dict(format="text", id=id)
        response = self._get_response("/api/dataset_info", params)
        return NexusDatasetInfo.FromMultiDict(self._parse_text_response(response))

    def list(self, tags=None, operator="or", order="date"):
        """Retrieves a list of datasets matching a set of tags from Nexus.

        @param tags: the tags the returned datasets should have. C{None}
          retrieves all the datasets, a single string retrieves datasets
          having that given tag. Multiple tags may also be specified as
          a list, tuple or any other iterable.
        @param operator: when multiple tags are given, this argument
          specifies whether the retrieved datasets should match all
          the tags (C{"and"}) or any of them (C{"or"}).
        @param order: the order of entries; it must be one of C{"date"},
          C{"name"} or C{"popularity"}.
        @return: a L{NexusDatasetInfoList} object, which basically acts like a
          list and yields L{NexusDatasetInfo} objects. The list is populated
          lazily; i.e. the requests will be fired only when needed.
        """
        params = dict(format="text", order=order)
        if tags is not None:
            if not hasattr(tags, "__iter__") or isinstance(tags, basestring):
                params["tag"] = str(tags)
            else:
                params["tag"] = "|".join(str(tag) for tag in tags)
                params["operator"] = operator

        return NexusDatasetInfoList(self, "/api/dataset_info", params)

    def search(self, query, order="date"):
        """Retrieves a list of datasets matching a query string from Nexus.

        @param query: the query string. Searches are case insensitive and
          Nexus searches for complete words only. The special word OR
          can be used to find datasets that contain any of the given words
          (instead of all of them). Exact phrases must be enclosed in
          quotes in the search string. See the Nexus webpage for more
          information at L{http://nexus.igraph.org/web/docs#searching}.
        @param order: the order of entries; it must be one of C{"date"},
          C{"name"} or C{"popularity"}.
        @return: a L{NexusDatasetInfoList} object, which basically acts like a
          list and yields L{NexusDatasetInfo} objects. The list is populated
          lazily; i.e. the requests will be fired only when needed.
        """
        params = dict(q=query, order=order, format="text")
        return NexusDatasetInfoList(self, "/api/search", params)

    def _get_response(self, path, params={}):
        """Sends a request to Nexus at the given path with the given parameters
        and returns a file-like object for the response."""
        url = "%s%s?%s" % (self.url, path, urlencode(params))
        return urlopen(url)

    @staticmethod
    def _parse_text_response(response):
        """Parses a plain text formatted response from Nexus.

        Plain text formatted responses consist of key-value pairs, separated
        by C{":"}. Values may span multiple lines; in this case, the key is
        omitted after the first line and the extra lines start with
        whitespace.

        Examples:
        
            >>> d = Nexus._parse_text_response("Id: 17\\nName: foo")
            >>> sorted(d.items())
            [('Id', '17'), ('Name', 'foo')]
            >>> d = Nexus._parse_text_response("Id: 42\\nName: foo\\n  .\\n bar")
            >>> sorted(d.items())
            [('Id', '42'), ('Name', 'foo\\n\\nbar')]
        """
        if isinstance(response, basestring):
            response = response.split("\n")

        result = multidict()
        key, value = None, []
        for line in response:
            line = line.rstrip()
            if key is not None and line[0] in ' \t':
                # Line continuation
                line = line.lstrip()
                if line == '.':
                    line = ''
                value.append(line)
            else:
                # Key-value pair
                if key is not None:
                    result.add(key, "\n".join(value))
                key, value = line.split(":", 1)
                value = [value.strip()]

        if key is not None:
            result.add(key, "\n".join(value))

        return result

    @property
    def url(self):
        """Returns the root URL of the Nexus repository the connection is
        communicating with."""
        return self._url

    @url.setter
    def url(self, value):
        """Sets the root URL of the Nexus repository the connection is
        communicating with."""
        value = str(value)
        parts = urlparse(value, "http", False)
        self._url = urlunparse(parts)
        if self._url and self._url[-1] == "/":
            self._url = self._url[:-1]


class NexusDatasetInfo(object):
    """Information about a dataset in the Nexus repository."""

    def __init__(self, id=None, sid=None, name=None, vertices=None,
            edges=None, tags=None, attributes=None, rest=None):
        self._conn = None
        self.id = id
        self.sid = sid
        self.name = name
        self.vertices = vertices
        self.edges = edges
        self.tags = tags
        self.attributes = attributes
        if rest:
            self.rest = multidict(rest)
        else:
            self.rest = None

    def __repr__(self):
        params = "(id=%(id)r, sid=%(sid)r, name=%(name)r, vertices=%(vertices)r, "\
                 "edges=%(edges)r, tags=%(tags)r, attributes=%(attributes)r, "\
                 "rest=%(rest)r)" % self.__dict__
        return "%s%s" % (self.__class__.__name__, params)

    def __str__(self):
        lines = ["[Nexus dataset #%s, %s/%s, %s]" % (self.id, self.vertices, self.edges,
            self.name)]
        if self.tags:
            lines.append("Tags: %s" % "; ".join(self.tags))

        if self.rest:
            wrapper = TextWrapper(width=76, subsequent_indent='  ')

            keys = sorted(self.rest.iterkeys())
            if "Attribute" in self.rest:
                keys.remove("Attribute")
                keys.append("Attribute")

            for key in keys:
                for value in self.rest.getlist(key):
                    paragraphs = str(value).splitlines()
                    wrapper.initial_indent = "%s: " % key
                    for paragraph in paragraphs:
                        ls = wrapper.wrap(paragraph)
                        if ls:
                            lines.extend(wrapper.wrap(paragraph))
                        else:
                            lines.append("  .")
                        wrapper.initial_indent = "  "

        return "\n".join(lines)

    def _update_from_multidict(self, params):
        """Updates the dataset object from a multidict representation of
        key-value pairs, similar to the ones provided by the Nexus API in
        plain text response."""
        self.id = params.get("Id", None)
        self.sid = params.get("Sid", None)
        self.name = params.get("Name", None)
        self.vertices = params.get("Vertices", None)
        self.edges = params.get("Edges", None)
        self.tags = params.get("Tags", None)

        if self.rest is None:
            self.rest = multidict()
        for k in set(params.iterkeys()) - set("Id Name Vertices Edges Tags".split()):
            for v in params.getlist(k):
                self.rest.add(k, v)

        if self.id:
            self.id = int(self.id)
        if self.vertices:
            self.vertices = int(self.vertices)
        if self.edges:
            self.edges = int(self.edges)
        if self.tags is not None:
            self.tags = self.tags.split(";")

    @classmethod
    def FromMultiDict(cls, dict):
        """Constructs a Nexus dataset object from a multidict representation
        of key-value pairs, similar to the ones provided by the Nexus API in
        plain text response."""
        result = cls()
        result._update_from_multidict(dict)
        return result

    def download(self):
        """Retrieves the actual dataset from Nexus.

        Returns a L{Graph} instance."""
        if self.id is None:
            raise ValueError("dataset ID is empty")
        conn = self._conn or Nexus
        return conn.get(self.id)

    get = download


class NexusDatasetInfoList(object):
    """A read-only list-like object that can be used to retrieve the items
    from a Nexus search result.
    """

    def __init__(self, connection, method, params):
        """Constructs a Nexus dataset list that will use the given connection
        and the given parameters to retrieve the search results.

        @param connection: a Nexus connection object
        @param method: the URL of the Nexus API method to call
        @param params: the parameters to pass in the GET requests, in the
          form of a Python dictionary.
        """
        self._conn = connection
        self._method = str(method)
        self._params = params
        self._length = None
        self._datasets = []
        self._blocksize = 10

    def _fetch_results(self, index):
        """Fetches the results from Nexus such that the result item with the
        given index will be available (unless the result list is shorter than
        the given index of course)."""
        # Calculate the start offset
        page = index // self._blocksize
        offset = page * self._blocksize
        self._params["offset"] = offset
        self._params["limit"] = self._blocksize

        # Ensure that self._datasets has the necessary length
        diff = (page+1) * self._blocksize - len(self._datasets)
        if diff > 0:
            self._datasets.extend([None] * diff)

        response = self._conn._get_response(self._method, self._params)
        current_dataset = None
        for line in response:
            key, value = line.strip().split(": ", 1)

            if key == "Totalsize":
                # Total number of items in the search result
                self._length = int(value)
            elif key == "Id":
                # Starting a new dataset
                if current_dataset:
                    self._datasets[offset] = current_dataset
                    offset += 1
                current_dataset = NexusDatasetInfo(id=int(value))
                current_dataset._conn = self._conn
            elif key == "Sid":
                current_dataset.sid = value
            elif key == "Name":
                current_dataset.name = value
            elif key == "Vertices":
                current_dataset.vertices = int(value)
            elif key == "Edges":
                current_dataset.edges = int(value)
            elif key == "Tags":
                current_dataset.tags = value.split(";")

        if current_dataset:
            self._datasets[offset] = current_dataset


    def __getitem__(self, index):
        if len(self._datasets) <= index:
            self._fetch_results(index)
        elif self._datasets[index] is None:
            self._fetch_results(index)
        return self._datasets[index]

    def __iter__(self):
        for i in xrange(len(self)):
            yield self[i]

    def __len__(self):
        """Returns the number of result items."""
        if self._length is None:
            self._fetch_results(0)
        return self._length

Nexus = NexusConnection()