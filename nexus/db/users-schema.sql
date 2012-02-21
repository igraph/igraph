CREATE TABLE session (
       session_id CHAR(128) UNIQUE NOT NULL,
       atime TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
       data TEXT
);

CREATE TABLE users (
       username TEXT UNIQUE NOT NULL,
       password CHAR(32) NOT NULL,
       admin BOOLEAN NOT NULL DEFAULT False
);

INSERT INTO users (username, password, admin) 
       VALUES ('csardi', '', 1);

INSERT INTO users (username, password, admin)
       VALUES ('testuser', '', 0);

