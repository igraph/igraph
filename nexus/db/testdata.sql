-- Test data

INSERT INTO licence (id, name, text, link) VALUES(1, 'CC-ASA 2.0 UK', 
       'Creative Commons Attribution-ShareAlike 2.0 UK: England & Wales',
       'http://creativecommons.org/licenses/by-sa/2.0/uk/');
INSERT INTO licence (id, name, text, link) VALUES(2, 'Public Domain',
       'In the Public Domain', 'http://en.wikipedia.org/wiki/Public_domain');

INSERT INTO dataset (id, name, shortdescription, description, licence, 
                     source) 
       VALUES(1, 'Zachary''s karate club',
       'Social network of a karate club',
       'Social network of a karate club', 2, 'Publication');
INSERT INTO dataset (id, name, shortdescription, description, licence, 
                     source) 
       VALUES(2, 'UK faculty social network',
       'Social network among faculty members at a UK university',
       'Social network among faculty members at a UK university', 1,
       'http://www.cs.rhul.ac.uk/home/tamas/research/datasets/');

INSERT INTO network (dataset, id, description, vertices, edges, filename)
       VALUES(1, 1, NULL, 34, 78, 'karate');
INSERT INTO network (dataset, id, description, vertices, edges, filename)
       VALUES(2, 1, NULL, 81, 817, 'UKfaculty');

INSERT INTO citation VALUES(1, 'W. W. Zachary, An information flow model for conflict and fission in small groups, Journal of Anthropological Research 33, 452-473 (1977)');

INSERT INTO citation VALUES(2, 'Nepusz T., Petróczi A., Négyessy L., Bazsó F.: Fuzzy communities and the concept of bridgeness in complex networks. Physical Review E 77:016107, (2008)');

INSERT INTO dataset_citation VALUES(1, 1);
INSERT INTO dataset_citation VALUES(2, 2);

INSERT INTO tag VALUES(1, 'weighted', 'The network has weighted edges.');
INSERT INTO tag VALUES(2, 'bipartite', 
       'Bipartite, in other words, two-mode network.');
INSERT INTO tag VALUES(3, 'undirected', 'Edges do not have directions.');
INSERT INTO tag VALUES(4, 'directed', 'Edges have directions.');

INSERT INTO dataset_tag VALUES(1, 1);
INSERT INTO dataset_tag VALUES(1, 3);
INSERT INTO dataset_tag VALUES(2, 1);
INSERT INTO dataset_tag VALUES(2, 4);

INSERT INTO metadata values (1, 1, 'vertex', 'string', 'name', 
       'Names of the actors. These are pseudonyms and only serve to identify the two distinguished actors of the network: Mr Hi and John A.');
INSERT INTO metadata values (1, 1, 'vertex', 'numeric', 'Faction',
       'Either 1 or 2, giving the faction of the actor after the split of the club.');
INSERT INTO metadata values (1, 1, 'edge', 'numeric', 'weight', 
       'Edge weights, the number of common activities the two club members took part of.');

INSERT INTO metadata values (2, 1, 'vertex', 'numeric', 'Group',
       'The numeric ID of the school affiliation');
INSERT INTO metadata values (2, 1, 'edge', 'numeric', 'weight',
       'Edge weights, based on the questionnaire');

INSERT INTO format VALUES('R-igraph', 
       'Rdata file for use with the igraph R package',
       'igraph is an R package for network analysis, with support for large graphs. This data file can be loaded into an R session via the ''load'' function, after loading the igraph package itself using the ''library'' function. See the igraph homepage for details.', 
       '.Rdata',
       'http://igraph.sourceforge.net');

INSERT INTO format VALUES('GraphML',
       'XML-based graph description language',
       'GraphML is an XML-based file format (an XML application in the XML terminology) to describe graphs. It is a modern format, and can store graphs with an extensible set of vertex and edge attributes', 
       '.graphml',
       'http://graphml.graphdrawing.org');

INSERT INTO format VALUES('Pajek',
       'File format of the program Pajek',
       'Pajek it a popular network analysis program for Windows. Note that the Pajek data format does not support arbitrary attributes, so some of the metadata might be missing from the data file if you choose this format.',
       '.net',
       'http://vlado.fmf.uni-lj.si/pub/networks/pajek/');

INSERT INTO format VALUES('Excel',
       'Excel workbook containing all network and metadata',
       'This is a standard file for MS Excel, and other tools capable of reading it, e.g. LibreOffice. The first worksheet includes the network itself and the edge metadata. The second worksheet contains the vertex metadata, if available. The third worksheet contains graph-level data.',
       '.xls',
       'http://office.microsoft.com/en-us/excel/');

INSERT INTO format VALUES('Python-igraph',
       'Pickled Python object for use with the igraph Python module',
       'igraph is a Python module for network analysis, with support for large graphs. This file can be loaded from a Python script using the ''load'' function of the ''pickle'' module if you have installed the ''igraph'' module before. See the igraph homepage for details.',
       'http://igraph.sourceforge.net');

INSERT INTO user VALUES('gabor.csardi', 
       'https://launchpad.net/~gabor.csardi', 1);
INSERT INTO user VALUES('ntamas', 'https://launchpad.net/~ntamas', 1);
