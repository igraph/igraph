# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""Interface to the Nexus online graph repository.

The classes in this file facilitate access to the Nexus online graph
repository at U{http://nexus.igraph.org}.

The main entry point of this package is the C{Nexus} variable, which is
an instance of L{NexusConnection}. Use L{NexusConnection.get} to get a particular
network from Nexus, L{NexusConnection.list} to list networks having a given set of
tags, L{NexusConnection.search} to search in the dataset descriptions, or
L{NexusConnection.info} to show the info sheet of a dataset."""

from cStringIO import StringIO
from gzip import GzipFile
from itertools import izip
from textwrap import TextWrapper
from urllib import urlencode
from urlparse import urlparse, urlunparse
from textwrap import TextWrapper

from igraph.compat import property
from igraph.configuration import Configuration
from igraph.utils import multidict

import re
import urllib2

__all__ = ["Nexus", "NexusConnection"]

__license__ = u"""\
Copyright (C) 2006-2012  Tamás Nepusz <ntamas@gmail.com>
Pázmány Péter sétány 1/a, 1117 Budapest, Hungary

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

class NexusConnection(object):
    """Connection to a remote Nexus server.

    In most cases, you will not have to instantiate this object, just use
    the global L{Nexus} variable which is an instance of L{NexusConnection}
    and connects to the Nexus repository at U{http://nexus.igraph.org}.

    Example:

      >>> print Nexus.info("karate")            #doctest:+ELLIPSIS
      Nexus dataset 'karate' (#1)
      vertices/edges: 34/78
      name: Zachary's karate club
      tags: social network; undirected; weighted
      ...
      >>> karate = Nexus.get("karate")
      >>> from igraph import summary
      >>> summary(karate)
      IGRAPH UNW- 34 78 -- Zachary's karate club network
      + attr: Author (g), Citation (g), name (g), Faction (v), id (v), name (v), weight (e)

    @undocumented: _get_response, _parse_dataset_id, _parse_text_response,
      _ensure_uncompressed"""

    def __init__(self, nexus_url=None):
        """Constructs a connection to a remote Nexus server.

        @param nexus_url: the root URL of the remote server. Leave it at its
          default value (C{None}) unless you have set up your own Nexus server
          and you want to connect to that. C{None} fetches the URL from
          igraph's configuration file or uses the default URL if no URL
          is specified in the configuration file.
        """
        self.debug = False
        self.url = nexus_url
        self._opener = urllib2.build_opener()

    def get(self, id):
        """Retrieves the dataset with the given ID from Nexus.

        Dataset IDs are formatted as follows: the name of a dataset on its own
        means that a single network should be returned if the dataset contains
        a single network, or multiple networks should be returned if the dataset
        contains multiple networks. When the name is followed by a dot and a
        network ID, only a single network will be returned: the one that has the
        given network ID. When the name is followed by a dot and a star, a
        dictionary mapping network IDs to networks will be returned even if the
        original dataset contains a single network only.

        E.g., getting C{"karate"} would return a single network since the
        Zachary karate club dataset contains one network only. Getting
        C{"karate.*"} on the other hand would return a dictionary with one
        entry that contains the Zachary karate club network.

        @param id: the ID of the dataset to retrieve.
        @return: an instance of L{Graph} (if a single graph has to be returned)
          or a dictionary mapping network IDs to instances of L{Graph}.
        """
        from igraph import load

        dataset_id, network_id = self._parse_dataset_id(id)

        params = dict(format="Python-igraph", id=dataset_id)
        response = self._get_response("/api/dataset", params, compressed=True)
        response = self._ensure_uncompressed(response)
        result = load(response, format="pickle")

        if network_id is None:
            # If result contains a single network only, return that network.
            # Otherwise return the whole dictionary
            if not isinstance(result, dict):
                return result
            if len(result) == 1:
                return result[result.keys()[0]]
            return result

        if network_id == "*":
            # Return a dict no matter what
            if not isinstance(result, dict):
                result = dict(dataset_id=result)
            return result

        return result[network_id]

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
          information at U{http://nexus.igraph.org/web/docs#searching}.
        @param order: the order of entries; it must be one of C{"date"},
          C{"name"} or C{"popularity"}.
        @return: a L{NexusDatasetInfoList} object, which basically acts like a
          list and yields L{NexusDatasetInfo} objects. The list is populated
          lazily; i.e. the requests will be fired only when needed.
        """
        params = dict(q=query, order=order, format="text")
        return NexusDatasetInfoList(self, "/api/search", params)

    @staticmethod
    def _ensure_uncompressed(response):
        """Expects an HTTP response object, checks its Content-Encoding header,
        decompresses the data and returns an in-memory buffer holding the
        uncompressed data."""
        compressed = response.headers.get("Content-Encoding") == "gzip"
        if not compressed:
            content_disp = response.headers.get("Content-Disposition", "")
            compressed = bool(re.match(r'attachment; *filename=.*\.gz\"?$',
                    content_disp))
        if compressed:
            return GzipFile(fileobj=StringIO(response.read()), mode="rb")
        print response.headers
        return response

    def _get_response(self, path, params={}, compressed=False):
        """Sends a request to Nexus at the given path with the given parameters
        and returns a file-like object for the response. `compressed` denotes
        whether we accept compressed responses."""
        if self.url is None:
            url = Configuration.instance()["remote.nexus.url"]
        else:
            url = self.url
        url = "%s%s?%s" % (url, path, urlencode(params))
        request = urllib2.Request(url)
        if compressed:
            request.add_header("Accept-Encoding", "gzip")
        if self.debug:
            print "[debug] Sending request: %s" % url
        return self._opener.open(request)

    @staticmethod
    def _parse_dataset_id(id):
        """Parses a dataset ID used in the `get` request.

        Returns the dataset ID and the network ID (the latter being C{None}
        if the original ID did not contain a network ID ).
        """
        dataset_id, _, network_id = str(id).partition(".")
        if not network_id:
            network_id = None
        return dataset_id, network_id

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
            if not line:
                continue
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
        if value is None:
            self._url = None
        else:
            value = str(value)
            parts = urlparse(value, "http", False)
            self._url = urlunparse(parts)
            if self._url and self._url[-1] == "/":
                self._url = self._url[:-1]


class NexusDatasetInfo(object):
    """Information about a dataset in the Nexus repository.
    
    @undocumented: _update_from_multidict, vertices_edges"""

    def __init__(self, id=None, sid=None, name=None, networks=None,
            vertices=None, edges=None, tags=None, attributes=None, rest=None):
        self._conn = None
        self.id = id
        self.sid = sid
        self.name = name
        self.vertices = vertices
        self.edges = edges
        self.tags = tags
        self.attributes = attributes
        if networks is None:
            self.networks = []
        elif not isinstance(networks, (str, unicode)):
            self.networks = list(networks)
        else:
            self.networks = [networks]
        if rest:
            self.rest = multidict(rest)
        else:
            self.rest = None

    @property
    def vertices_edges(self):
        if self.vertices is None or self.edges is None:
            return ""
        elif isinstance(self.vertices, (list, tuple)) and isinstance(self.edges, (list, tuple)):
            return " ".join("%s/%s" % (v,e) for v, e in izip(self.vertices, self.edges))
        else:
            return "%s/%s" % (self.vertices, self.edges)

    @vertices_edges.setter
    def vertices_edges(self, value):
        if value is None:
            self.vertices, self.edges = None, None
            return

        value = value.strip().split(" ")
        if len(value) == 0:
            self.vertices, self.edges = None, None
        elif len(value) == 1:
            self.vertices, self.edges = map(int, value[0].split("/"))
        else:
            self.vertices = []
            self.edges = []
            for ve in value:
                v, e = ve.split("/", 1)
                self.vertices.append(int(v))
                self.edges.append(int(e))

    def __repr__(self):
        params = "(id=%(id)r, sid=%(sid)r, name=%(name)r, networks=%(networks)r, "\
                "vertices=%(vertices)r, edges=%(edges)r, tags=%(tags)r, "\
                "attributes=%(attributes)r, rest=%(rest)r)" % self.__dict__
        return "%s%s" % (self.__class__.__name__, params)

    def __str__(self):
        if self.networks and len(self.networks) > 1:
            lines = ["Nexus dataset '%s' (#%s) with %d networks" % \
                    (self.sid, self.id, len(self.networks))]
        else:
            lines = ["Nexus dataset '%(sid)s' (#%(id)s)" % self.__dict__]

        lines.append("vertices/edges: %s" % self.vertices_edges)

        if self.name:
            lines.append("name: %s" % self.name)
        if self.tags:
            lines.append("tags: %s" % "; ".join(self.tags))

        if self.rest:
            wrapper = TextWrapper(width=76, subsequent_indent='  ')

            keys = sorted(self.rest.iterkeys())
            if "attribute" in self.rest:
                keys.remove("attribute")
                keys.append("attribute")

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
        self.id = params.get("id")
        self.sid = params.get("sid")
        self.name = params.get("name")
        self.vertices = params.get("vertices")
        self.edges = params.get("edges")
        self.tags = params.get("tags")

        networks = params.get("networks")
        if networks:
            self.networks = networks.split()

        keys_to_ignore = set("id sid name vertices edges tags networks".split())

        if self.vertices is None and self.edges is None:
            # Try "vertices/edges"
            self.vertices_edges = params.get("vertices/edges")
            keys_to_ignore.add("vertices/edges")

        if self.rest is None:
            self.rest = multidict()
        for k in set(params.iterkeys()) - keys_to_ignore:
            for v in params.getlist(k):
                self.rest.add(k, v)

        if self.id:
            self.id = int(self.id)
        if self.vertices and not isinstance(self.vertices, (list, tuple)):
            self.vertices = int(self.vertices)
        if self.edges and not isinstance(self.edges, (list, tuple)):
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

    def download(self, network_id=None):
        """Retrieves the actual dataset from Nexus.

        @param network_id: if the dataset contains multiple networks, the ID
          of the network to be retrieved. C{None} returns a single network if
          the dataset contains a single network, or a dictionary of networks
          if the dataset contains more than one network. C{"*"} retrieves
          a dictionary even if the dataset contains a single network only.

        @return: a L{Graph} instance or a dictionary mapping network names to
          L{Graph} instances.
        """
        if self.id is None:
            raise ValueError("dataset ID is empty")
        conn = self._conn or Nexus
        if network_id is None:
            return conn.get(self.id)
        return conn.get("%s.%s" % (self.id, network_id))

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
            key = key.lower()

            if key == "totalsize":
                # Total number of items in the search result
                self._length = int(value)
            elif key == "id":
                # Starting a new dataset
                if current_dataset:
                    self._datasets[offset] = current_dataset
                    offset += 1
                current_dataset = NexusDatasetInfo(id=int(value))
                current_dataset._conn = self._conn
            elif key == "sid":
                current_dataset.sid = value
            elif key == "name":
                current_dataset.name = value
            elif key == "vertices":
                current_dataset.vertices = int(value)
            elif key == "edges":
                current_dataset.edges = int(value)
            elif key == "vertices/edges":
                current_dataset.vertices_edges = value
            elif key == "tags":
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

    def __str__(self):
        """Converts the Nexus result list into a nice human-readable format."""
        max_index_length = len(str(len(self))) + 2
        indent = "\n" + " " * (max_index_length+1)

        result = []
        for index, item in enumerate(self):
            formatted_item = ("[%d]" % index).rjust(max_index_length) + " " + \
                str(item).replace("\n", indent)
            result.append(formatted_item)
        return "\n".join(result)

Nexus = NexusConnection()