import web
import os.path
from web import websafe

if os.path.exists("../db/nexus.db"):
    dbfile='../db/nexus.db'
elif os.path.exists("../db/test.db"):
    dbfile='../db/test.db'
db = web.database(dbn='sqlite', db=dbfile)

def get_list_tagged_as(tagname, operator="or", order="date", limit=10, 
                       offset=0):
    ord={ 'date': 'd.date DESC, d.id DESC',
          'name': 'd.name',
          'popularity': 'd.downloads DESC, d.date DESC, d.id DESC' }
    tagstr=["'%s'" % websafe(t) for t in tagname]
    tagstr='(' + ",".join(tagstr) + ")"
    if operator=="or":
        table='dataset d, dataset_tag dt, tag t, network n'
        what='''DISTINCT d.id id, d.name name,
                         d.description description,
                         d.shortdescription shortdescription,
                         d.licence licence, d.date date,
                         n.id netid, n.description netdescription, 
                         n.vertices vertices, n.edges edges'''
        where='''d.id=dt.dataset 
                 AND dt.tag=t.id 
                 AND t.tag IN %s 
                 AND d.id=n.dataset''' % tagstr
        res=db.select(table, what=what, where=where, order=ord[order],
                      offset=offset, limit=limit)
        count=list(db.select(table, what='COUNT(DISTINCT d.id) c',
                             where=where))[0].c
    else:
        table="dataset d, network n"
        what='''d.id id, d.name name,
                d.description description, 
                d.shortdescription shortdescription,
                d.licence licence,
                d.date date, n.id netid, 
                n.description netdescription, n.vertices vertices,
                n.edges edges'''
        where='''d.id IN (SELECT id FROM 
                                    (SELECT dt.dataset id, COUNT(*) count
                                     FROM dataset_tag dt, tag t
                                     WHERE dt.tag=t.id AND t.tag IN %s
                                     GROUP BY dt.dataset)
                                 WHERE count=%s)
                 AND d.id=n.dataset''' % (tagstr, len(tagname))
        res=db.select(table, what=what, where=where, order=ord[order],
                      offset=offset, limit=limit)
        count=list(db.select(table, what='COUNT(*) c', where=where))[0].c
    
    return list(res), count

def get_list_of_datasets(order="date", limit=10, offset=0):
    ord={ 'date': 'd.date DESC, d.id DESC',
          'name': 'd.name',
          'popularity': 'd.downloads DESC, d.date DESC, d.id DESC' }
    count=list(db.select('dataset', what='COUNT(id) count'))[0].count
    res=db.query('''SELECT d.id id, d.name name, 
                           d.description description, 
                           d.shortdescription shortdescription,
                           d.licence licence, d.date date, 
                           n.id netid, n.description netdescription, 
                           n.vertices vertices, n.edges edges
                    FROM dataset d, network n 
                    WHERE d.id=n.dataset
                    ORDER BY %s
                    LIMIT %s 
                    OFFSET %s''' % (ord[order], int(limit), int(offset)))
    res=list(res)
    return res, count

def get_tags(id):
    return db.query('''SELECT t.id, t.tag 
                       FROM tag t, dataset_tag d
                       WHERE t.id=d.tag AND dataset=$id
                       ORDER BY t.tag''', 
                    vars={'id': id})

def get_dataset(id):
    return db.query('''SELECT d.id id, d.name name, 
                              d.description description, 
                              d.shortdescription shortdescription,
                              l.id licence,
                              l.name licence_name,
                              l.link licence_url,
                              d.date date, n.id netid, 
                              n.description netdescription, 
                              n.vertices vertices, n.edges edges,
                              n.filename filename, d.source source
                       FROM dataset d, licence l, network n
                       WHERE d.licence=l.id
                         AND d.id=$id AND d.id=n.dataset''',
                    vars={'id': id})

def delete_dataset(id):
    db.delete('dataset', where='id=%s' % int(id))
    db.delete('network', where='dataset=%s' % int(id))
    db.delete('dataset_citation', where='dataset=%s' % int(id))
    db.delete('dataset_tag', where='dataset=%s' % int(id))
    db.delete('metadata', where='dataset=%s' % int(id))

def get_dataset_file(id, nid):
    return db.query('''SELECT filename FROM network 
                         WHERE dataset=$id AND id=$nid''',
                    vars={'id': id, 'nid': nid})

def whatsnew():
    return db.select('dataset', limit=5, order="date DESC, id DESC")

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
                       AND d.dataset=%s''' % websafe(id))

def get_format(name):
    return db.select('format', where="name='%s'" % websafe(name))

def get_licences(what=None):
    return db.select('licence', what=what, order="name")

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

def update_dataset(id, name, shortdescription,
                   description, tags, licence, vertices,
                   edges, filename, source, papers):

    db.update('dataset', where='id=%s' % int(id), name=name,
              shortdescription=shortdescription,
              description=description, licence=int(licence), 
              source=source)
    db.update('network', where='dataset=%s AND id=1' % int(id), 
              vertices=vertices, edges=edges, filename=filename)

    db.delete('dataset_tag', where='dataset=%s' % int(id))
    for t in tags:
        new_dataset_tag(id, t)

    db.delete('dataset_citation', where='dataset=%s' % int(id))
    for p in papers:
        new_citation(id, p)
    
    return True

def get_username(openid=None):
    if not openid:
        openid=web.webopenid.status()
    user=db.select('user', what='name', 
                   where="openid='%s'" % websafe(openid))
    if user:
        return [u for u in user][0].name
    else:
        return None

def get_licence(id):
    return db.select('licence', where="id='%s'" % int(id))

def update_licence(id, **args):
    db.update('licence', where='id=%s' % int(id), **args)
    return True

def get_metadata(id):
    res=db.select('metadata', where='dataset=%s' % int(id), 
                  order="type")
    return list(res)

def delete_meta(id, type, name):
    res=db.delete('metadata', where="type='%s' AND name='%s'" % \
                      (websafe(type), websafe(name)))

def update_meta(id, type, name, **args):
    res=db.update('metadata', where="type='%s' AND name='%s'" % 
                  (websafe(type), websafe(name)),
                  type=type, name=name, **args)

def add_meta(id, type, name, **args):
    res=db.insert('metadata', dataset=id, network=1, type=type, name=name,
                  **args)

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
    db.update('blog', where='id=%s' % int(id), **args)

def increase_downloads(id):
    db.query('''UPDATE dataset SET downloads=downloads+1 WHERE id=%s''' % \
                 int(id))

def get_totals():
    res=db.select('dataset', 
                  what='COUNT(*) AS count, SUM(downloads) AS downloads')
    return list(res)[0]

def get_format_extensions():
    res=db.select('format', what='name,extension')
    return dict( (f.name,f.extension) for f in res )

def get_format_extension(format):
    res=db.select('format', what='extension', 
                  where="name='%s'" % websafe(format))
    return list(res)[0].extension
