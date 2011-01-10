docs='''
<span>
# Contents
</span>

1. [Atom and RSS feeds](#feeds)
1. [The Nexus API](#api)
    1. [Basics](#basics)
    1. [Data set queries](#dsq)
    1. [Data set downloads](#dsd)
    1. [Data formats](#dformat)
    1. [Licences](#licence)
    1. [XML schemas](#schema)

<a name="feeds"></a>
# Atom and RSS feeds

Nexus provides [Atom][6] and [RSS 2.0][7] feeds, these are handy if
you want to get notified about new data sets. There is a feed for all
data sets, but you can also use feeds for data set tags, or even
combination of tags. To subscribe to the non-specific feed, click on
the yellow feed icon at the top right corner of any Nexus page, in the
'What's new?' box.

Not all browsers support Atom and RSS feeds out of the box, e.g. for
the current version of Google Chrome, you need to install an extension
first. 

For subscribing to feeds of data set tags, e.g. if you are only
interested in weighted networks, click on the tag first in the tag
cloud at the top right corner of any Nexus page, in the 'Data tags'
box. Then click on the big yellow feed icon at the top of the page.

If you want to use RSS feeds instead of Atom feeds, then copy the feed
URL and replace the '<code>format=atom</code>' part with
'<code>format=rss</code>'.

<a name="api"></a>
# The Nexus API

Nexus features a web service an an API that allows programatic queries
and data downloads. The API is supposed to be [RESTful][1].

<a name="basics"></a>
## Basics

All queries can be supplied as HTTP GET methods, i.e. simple HTTP
downloads. A Nexus query has the following parts:

<table>
  <tr><td></td><td>Example</td></tr>
  <tr><td>Base URI</td><td><code>http://nexus.igraph.org/api/</code></td></tr>
  <tr><td>Service</td><td><code>dataset_info?</code></td></tr>
  <tr><td>Format</td><td><code>format=xml</code></td></tr>
  <tr><td>Arguments</td><td><code>&id=1</code></td></tr>
</table>

This gives the URI:
<code>http://nexus.igraph.org/api/dataset_info?format=xml&id=1</code>.

Currently Nexus supports four services: 

<table>
  <tr><td><code>dataset_info?</code></td><td>Information about one or
    more data sets.</td></tr>
  <tr><td><code>dataset?</code></td><td>Download a data set.</td></tr>
  <tr><td><code>format?</code></td><td>Information about data formats.
    </td></tr>
  <tr><td><code>licence?</code></td><td>Information about data set
    licences.</td></tr>
</table>

and two data formats, XML and plain text.

<a name="dsq"></a>
## Data set queries

### Parameters

Data set information can be queried using the
'<code>dataset_info?</code>' service. It has the following parameters:

<table>
  <tr><td>id</td><td>The id of the data set to query. If not present,
    then a list is returned, with basic information about all data
    sets. (Or a subset of them, if the '<code>tag</code>' parameter is
    present.) If '<code>id</code>' is given, then detailed information
    is returned about a single dataset. See examples below.</td></tr>
  <tr><td>tag</td><td>Tag or tags separated by a '<code>|</code>'
    character. Only data sets with these tag(s) are listed. You cannot
    give both '<code>id</code>', and '<code>tag</code>', this results an
    error.</td></tr>
  <tr><td>operator</td><td>Specifies whether only one or all tags are
    required. If it is not given or '<code>or</code>' then all data sets
    are listed that have at least one of the given tags. If it is
    '<code>and</code>', then all tags are required for a data set to
    include it in the result.</td></tr>
  <tr><td>order</td><td>Specifies the order of the data sets, in the
    result. Possible values are in the table below.</td></tr>
</table>

Possible values for the '<code>order</code>' parameter:

<table>
  <tr><td>date</td><td>This is the default, newest data sets are
    first. If two data sets were added on the same day, then the one
    with the highest id is first. (This is most likely the newer of
    the two.)</td></tr>
  <tr><td>name</td><td>Alphabetical ordering, according to data sets
    name.</td></tr>
  <tr><td>popularity</td><td>Data sets that were downloaded many times
    are first.</td></tr>
</table>

### Examples

We give a couple of example Nexus queries here. You can click on the
URIs to see the actual response of the Nexus server. Note that some
browsers tend not to show the actual XML code of the response, in this
case try saving the response to a file. Please refer to the XML
schemas below for the exact syntax.

#### List all data sets, XML format

<code>[http://nexus.igraph.org/api/dataset_info?format=xml](http://nexus.igraph.org/api/dataset_info?format=xml)</code>

#### Information about data set #2, XML format

<code>[http://nexus.igraph.org/api/dataset_info?format=xml&id=2](http://nexus.igraph.org/api/dataset_info?format=xml&id=2)</code>

#### List all data sets tagged as '<code>weighted</code>' and '<code>directed</code>', XML format

<code>[http://nexus.igraph.org/api/dataset_info?format=xml&tag=weighted|directed&operator=and](http://nexus.igraph.org/api/dataset_info?format=xml&tag=weighted|directed&operator=and)</code>

#### List all data sets, plain text format, order by popularity

<code>[http://nexus.igraph.org/api/dataset_info?format=text&order=popularity](http://nexus.igraph.org/api/dataset_info?format=text&order=popularity)</code>

#### Information about data set #2, plain text format

<code>[http://nexus.igraph.org/api/dataset_info?format=text&id=2](http://nexus.igraph.org/api/dataset_info?format=text&id=2)</code>

#### List all data sets tagged as '<code>weighted</code>' and '<code>directed</code>', plain text format

<code>[http://nexus.igraph.org/api/dataset_info?format=text&tag=weighted|directed&operator=and](http://nexus.igraph.org/api/dataset_info?format=text&tag=weighted|directed&operator=and)</code>

<a name="dsd"></a>
## Data set downloads

### Parameters

For downloading a data set, one has to use the '<code>dataset?</code>'
service. Parameters:

<table>
  <tr><td>id</td><td>Numeric id of the data set to download. This must
    be given.</td></tr>
  <tr><td>format</td><td>The data format. Nexus supports downloads in
    different formats. For an up to date list, see the
    '<code>format?</code>' service in the next Section.</td></tr>
</table>

### Examples

#### Download data set #1 in Excel format

<code>[http://nexus.igraph.org/api/dataset?format=Excel&id=1](http://nexus.igraph.org/api/dataset?format=Excel&id=1)</code>

<a name="dformat"></a>
## Data formats

### Supported formats

See the current list of supported formats [here][8].

### Parameters

This service returns the data formats supported by Nexus. It has one
optional parameter:

<table>
  <tr><td>dataformat</td><td>The data format to query. If this is not
    given then information is returned about all data formats.</td></tr>
</table>

### Examples

#### Information about all supported data formats

<code>[http://nexus.igraph.org/api/format?format=xml](http://nexus.igraph.org/api/format?format=xml)</code>

#### Information about a single format, in plain text

<code>[http://nexus.igraph.org/api/format?format=text&dataformat=Excel](http://nexus.igraph.org/api/format?format=text&dataformat=Excel)</code>

<a name="licence"></a>
## Licences

### List of licences

The actual list of licences can be found [here][9].

### Parameters

This service returns information about the different data set
licences. It has one optional parameter:

<table>
  <tr><td>id</td><td>The id of the licence to query. If not given,
    then all licences are listed.</td></tr>
</table>

### Examples

#### All licences, plain text format

<code>[http://nexus.igraph.org/api/licence?format=text](http://nexus.igraph.org/api/licence?format=text)</code>

#### A single licence, XML format

<code>[http://nexus.igraph.org/api/licence?format=xml&id=1](http://nexus.igraph.org/api/licence?format=xml&id=1)</code>

<a name="schema"></a>
## XML schemas

### Basic data set information

Basic data set information, for many data sets. The
'<code>dataset_info?</code>' service returns this format, if the
'<code>id</code>' parameter is not given. The schema is [here][2].

### Detailed data set information
 
All available information about a single data set, including a long
description. The '<code>dataset_info?</code>' service returns this
format, if the '<code>id</code>' parameter is given. The schema is
[here][3].

### Data formats

The '<code>format?</code>' service returns this format. The schema is
[here][4].

### Licences

The '<code>licence?</code>' service returns this format. The schema is
[here][5].

[1]: http://en.wikipedia.org/wiki/Representational_State_Transfer 
[2]: /static/index_schema.xsl
[3]: /static/dataset_schema.xsl
[4]: /static/format_schema.xsl
[5]: /static/licence_schema.xsl
[6]: http://en.wikipedia.org/wiki/Atom_(standard)
[7]: http://en.wikipedia.org/wiki/RSS
[8]: /api/format
[9]: /api/licence
'''
