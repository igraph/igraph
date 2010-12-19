-- Test data

INSERT INTO licence (id, name, text) VALUES(1, 'CC-ASA 2.0 UK', 
       'Creative Commons Attribution-ShareAlike 2.0 UK: England & Wales');
INSERT INTO licence (id, name, text) VALUES(2, 'Public Domain',
       'In the Public Domain');

INSERT INTO dataset (id, name, description, licence, source) 
       VALUES(1, 'Zachary''s karate club',
       'Social network of a karate club', 2, 'Publication');
INSERT INTO dataset (id, name, description, licence, source) 
       VALUES(2, 'UK faculty social network',
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

INSERT INTO format VALUES('R-igraph', 
       'Rdata file for use with the igraph R package',
       'igraph is an R package for network analysis, with support for large graphs. This data file can be loaded into an R session via the ''load'' function, after loading the igraph package itself using the ''library'' function. See the igraph homepage for details.', 
       'http://igraph.sourceforge.net');
	
INSERT INTO user VALUES('gabor.csardi', 
       'https://launchpad.net/~gabor.csardi', 1);
INSERT INTO user VALUES('ntamas', 'https://launchpad.net/~ntamas', 1);
