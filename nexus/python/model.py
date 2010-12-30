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
                       WHERE t.id=d.tag AND dataset=$id
                       ORDER BY t.tag''', 
                    vars={'id': id})

def get_dataset(id):
    return db.query('''SELECT d.id id, d.name name, 
                              d.description description, l.id licence,
                              l.name licence_name,
                              d.date date, n.id netid, 
                              n.description netdescription, 
                              n.vertices vertices, n.edges edges,
                              n.filename, d.source source
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

def get_licences(what=None):
    return db.select('licence', what=what)

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

def update_dataset(id, name, description, tags, licence, vertices,
                   edges, filename, source, papers):

    db.update('dataset', where='id=%s' % id, name=name, 
              description=description, licence=int(licence), 
              source=source)
    db.update('network', where='dataset=%s AND id=1' % id, 
              vertices=vertices, edges=edges, filename=filename)

    db.delete('dataset_tag', where='dataset=%s' % id)
    for t in tags:
        new_dataset_tag(id, t)

    db.delete('dataset_citation', where='dataset=%s' % id)
    for p in papers:
        new_citation(id, p)
    
    return True

def get_username(openid=None):
    if not openid:
        openid=web.webopenid.status()
    user=db.select('user', what='name', 
                   where="openid='%s'" % openid)
    if user:
        return [u for u in user][0].name
    else:
        return None

def get_licence(id):
    return db.select('licence', where="id='%s'" % id)

def update_licence(id, **args):
    db.update('licence', where='id=%s' % id, **args)
    return True

def get_metadata(id):
    res=db.select('metadata', where='dataset=%s' % int(id), 
                  order="type")
    return list(res)

def get_blog(ids=None, unpublished=False):
    if ids is None:
        limit=5
        order='date DESC'
        where='1=1'
    else:
        limit=None
        order='id'
        where="id in (" + ",".join(ids) + ")"

    if not unpublished:
        where=where + " AND published=1"

    return list(db.select('blog', limit=limit, order=order, where=where))

def new_blog_entry(**args):
    return db.insert('blog', seqname='id', **args)

def update_blog(id, **args):
    db.update('blog', where='id=%s' % id, **args)
