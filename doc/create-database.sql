
# First create database igraph & user igraph and grant 
# all rights for igraph user to database igraph
# Then execute the commands below as user igraph

create table graph(gname SERIAL PRIMARY KEY);

create table gal (gname       integer not null references graph(gname) on delete cascade,
       	     	  name        varchar(30),
		  value       text,
		  primary key (gname, name));

create index gal_gname_idx on gal (gname);
create index gal_name_idx  on gal (name);

create table data (gname       integer not null references graph(gname) on delete cascade, 
       	     	   v1          integer not null, 
                   v2          integer not null,
		   eid	       SERIAL,
		   PRIMARY KEY (gname, eid));
		   
create index data_gname_idx on data (gname);
create index data_v1_idx on data (v1);
create index data_v2_idx on data (v2);
create index data_v1v2_idx on data (v1, v2);

create table val (gname	       integer not null references graph(gname) on delete cascade, 
       	          v            integer not null, 
       	          name 	       varchar(30) not null, 
		  value	       text,
		  PRIMARY KEY(gname, v, name));
		  
create index val_gname_idx on val (gname);		  
create index val_v_idx     on val (v);
create index val_name_idx  on val (name);		  

create table eal (gname	       integer,
		  eid	       integer,
		  name 	       varchar(30) not null, 
		  value	       text,
		  FOREIGN KEY (gname,eid) REFERENCES data(gname, eid) on delete cascade,
		  PRIMARY KEY (gname,eid,name));

create index eal_gname_idx on eal (gname);
create index eal_eid_idx   on eal (eid);		  

# Erasing the database

drop table graph cascade;
drop table gal cascade;
drop table data cascade;
drop table val cascade;
drop table eal cascade;
