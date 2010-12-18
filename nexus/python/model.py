import web
import os.path

if os.path.exists("../db/nexus.db"):
    dbfile='../db/nexus.db'
elif os.path.exists("../db/test.db"):
    dbfile='../db/test.db'
db = web.database(dbn='sqlite', db=dbfile)

def get_list_tagged_as(tagname):
    return db.query('''SELECT d.id id, d.name name,
                              d.description description, 
                              d.licence licence, d.date date,
                              n.id netid, n.description netdescription, 
                              n.vertices vertices, n.edges edges
                       FROM dataset d, dataset_tag dt, tag t, network n
                       WHERE d.id=dt.dataset
                         AND dt.tag=t.id
                         AND t.tag=$tagname
                         AND d.id=n.dataset
                       ORDER BY d.date DESC''',
                    vars={'tagname': tagname})

def get_list_of_datasets():
    return db.query('''SELECT d.id id, d.name name, 
                              d.description description, 
                              d.licence licence, d.date date, 
                              n.id netid, n.description netdescription, 
                              n.vertices vertices, n.edges edges
                       FROM dataset d, network n 
                       WHERE d.id=n.dataset
                       ORDER BY d.date DESC''')

def get_tags(id):
    return db.query('''SELECT t.id, t.tag 
                       FROM tag t, dataset_tag d
                       WHERE t.id=d.tag AND dataset=$id''', 
                    vars={'id': id})

def get_dataset(id):
    return db.query('''SELECT d.id id, d.name name, 
                              d.description description, l.name licence,
                              d.date date, n.id netid, 
                              n.description netdescription, 
                              n.vertices vertices, n.edges edges
                       FROM dataset d, licence l, network n
                       WHERE d.licence=l.id
                         AND d.id=$id AND d.id=n.dataset''',
                    vars={'id': id})

def get_dataset_file(id, nid):
    return db.query('''SELECT filename FROM network 
                         WHERE dataset=$id AND id=$nid''',
                    vars={'id': id, 'nid': nid})

def whatsnew():
    return db.select('dataset', limit=5, order="date DESC")

def datatags():
    return db.query('''SELECT t.tag tag, COUNT(*) count
                       FROM dataset_tag dt, tag t 
                       WHERE dt.tag=t.id 
                       GROUP BY dt.tag
                       ORDER BY count DESC;''')

def get_formats():
    return db.select('format')

def get_papers(id):
    return db.query('''SELECT c.citation citation
                       FROM citation c, dataset_citation d
                       WHERE c.id=d.citation 
                       AND d.dataset=%s''' % web.safestr(id))

def get_format(name):
    return db.select('format', where="name='%s'" % web.safestr(name))

def get_licences():
    return db.select('licence')

def new_dataset(**args):
    return db.insert('dataset', seqname='id', **args)

def new_network(**args):
    db.insert('network', seqname=None, **args)

def new_citation(dataset, citation):
    citid=db.insert('citation', seqname='id', citation=citation)
    db.insert('dataset_citation', seqname=None, 
              dataset=dataset, citation=citid)

def new_dataset_tag(dataset, tag):
    ex=db.select('tag', what='id', where='tag=$tag', 
                 vars={ 'tag': tag })
    ex=[e for e in ex]
    if not ex:
        tagid=db.insert('tag', seqname='id', tag=tag)
    else:
        tagid=ex[0].id
    db.insert('dataset_tag', seqname=None, dataset=dataset, tag=tagid)

def list_data_formats():
    return db.select('format')

def new_licence(**args):
    return db.insert('licence', seqname="id", **args)
