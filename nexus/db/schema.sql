CREATE TABLE licence (
       id INTEGER,
       name TEXT,
       text TEXT,	
       fulltext TEXT,
       link TEXT,
       PRIMARY KEY(id)
);

CREATE TABLE dataset (
       id INTEGER,
       name TEXT,
       description TEXT,
       licence INTEGER,
       source TEXT,
       date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
       PRIMARY KEY(id),
       FOREIGN KEY(licence) REFERENCES licence(id)
);

CREATE TABLE network (
       dataset INTEGER,
       id INTEGER,       
       description TEXT,
       vertices INTEGER,
       edges INTEGER,
       filename TEXT,
       date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
       PRIMARY KEY(dataset, id)
);

CREATE TABLE citation (
       id INTEGER,
       citation TEXT,
       PRIMARY KEY (id)
);

CREATE TABLE dataset_citation (
       dataset INTEGER,
       citation INTEGER,
       PRIMARY KEY(dataset, citation),
       FOREIGN KEY(dataset) REFERENCES dataset(id),
       FOREIGN KEY(citation) REFERENCES citation(id)
);

CREATE TABLE tag (
       id INTEGER,
       tag TEXT,
       description TEXT,
       PRIMARY KEY(id)
);

CREATE TABLE dataset_tag (
       dataset INTEGER,
       tag INTEGER,
       PRIMARY KEY(dataset, tag),
       FOREIGN KEY(dataset) REFERENCES dataset(id),
       FOREIGN KEY(tag) REFERENCES tag(id)
);

CREATE TABLE format (
       name VARCHAR(20),
       shortdesc VARCHAR(100),
       description TEXT,
       link VARCHAR(200),
       PRIMARY KEY(name)
);

CREATE TABLE metadata (
       dataset INTEGER,
       network INTEGER,
       type VARCHAR(10),
       datatype VARCHAR(10),
       name VARCHAR(30),
       description TEXT,
       PRIMARY KEY(dataset, network, type, name),
       FOREIGN KEY(dataset, network) REFERENCES network(dataset, network)
);
CREATE TABLE user (
       name TEXT,
       openid TEXT,
       admin BOOLEAN,
       PRIMARY KEY(openid),
       UNIQUE(name)
);
