#! /bin/env /home/csardi/software/bin/mypython
# vim:set ts=4 sw=4 sts=4 et

import sys
sys.path.append("../python")

import web
import model
import math
import random
import openid.store.filestore
import os
import url_helper
import odict
import subprocess
import markdown
import igraph
import xlwt
import xlrd

from datetime import datetime
from functools import wraps
from itertools import izip
from operator import attrgetter
from recaptcha.client import captcha
from textwrap import dedent

# web.config.debug = False
web.config.smtp_server = '173.192.111.8'
web.config.smtp_port = 26
web.config.smtp_username = 'csardi@mail.igraph.org'
web.config.smtp_password = file("../config/emailpass").read().strip()
web.config.smtp_starttls = True

urls = (
    '/?',                                  'About',
    '/api/dataset_info',                   'Index',
    '/api/dataset',                        'Dataset',
    '/api/format',                         'Format',
    '/api/licence',                        'Licence',
    '/web/about',                          'About',
    '/web/addblog',                        'AddBlog',
    '/web/addlicence',                     'AddLicence',
    '/web/add',                            'Add',
    '/web/admin',                          'Admin',
    '/web/blog',                           'Blog',
    '/web/check/(\d+)',                    'Check',
    '/web/delete/(\d+)',                   'Delete',
    '/web/docs',                           'Docs',
    '/web/donate',                         'Donate',
    '/web/editblog/(\d+)',                 'EditBlog',
    '/web/editlicence/(\d+)',              'EditLicence',
    '/web/edit/(\d+)',                     'Edit',
    '/web/feedback',                       'Feedback',
    '/web/loginfailed',                    'LoginFailed',
    '/web/logout',                         'Logout',
    '/web/openid',                         'OpenID',
    '/web/recreate/(\d+)',                 'Recreate',
    '.*',                                  'NotFound'    
    )

# reCAPTCHA keys
recaptcha_pubkey = "6Lfqjb8SAAAAAJyGZQrvqgs7JWjpNp_Vf9dpTMxy"

# OpenIDs of the site administrators
admins_openid=('https://launchpad.net/~gabor.csardi',
               'https://launchpad.net/~ntamas')

# Text to be included in the HTML output whenever a reCAPTCHA is needed
recaptcha_text = """
<div class="recaptcha">
<script type="text/javascript">
   var RecaptchaOptions = {
      theme : 'white'
   };
</script>
<script type="text/javascript"
        src="http://www.google.com/recaptcha/api/challenge?k=%s">
</script>
<noscript>
  <iframe src="http://www.google.com/recaptcha/api/noscript?k=%s"
          height="300" width="500" frameborder="0"></iframe><br>
  <textarea name="recaptcha_challenge_field" rows="3" cols="40">
  </textarea>
  <input type="hidden" name="recaptcha_response_field"
         value="manual_challenge">
</noscript>
</div>
""" % (recaptcha_pubkey, recaptcha_pubkey)

# List of supported webpage formats
formats = ('html', 'xml', 'text', 'rss', 'atom')

# Dictionary mapping supported dataset formats to extensions
dataformats = { 'R-igraph': '.Rdata', 
                'GraphML': '.graphml',
                'Pajek': '.net',
                'Excel': '.xls'        }

def get_current_url():
    """Returns the URL of the current page being produced"""
    return web.ctx.fullpath

def get_current_not_logout_url():
    """Returns the URL of the current page being produced, except for the
    logout page.

    For the logout and login failure pages, returns the root page."""
    fp=web.ctx.fullpath
    if fp in ('/web/logout', '/web/loginfailed'):
        fp='/'
    return fp

def get_whatsnew():
    """Retrieves the new datasets from the model and renders it."""
    w=model.whatsnew()
    return render_plain.whatsnew(w)

def get_datatags():
    """Retrieves the list of tags and tag counts in the model and renders a
    nice tag cloud."""
    tags = [tag for tag in model.datatags()]
    counts = [tag.count for tag in tags]
    max_count = max(counts)
    sizes = [ int(math.log(count+1) / math.log(max_count+1) * 5) for count in counts]
    for tag, size in izip(tags, sizes):
        tag.count = size
    tags.sort(key = attrgetter("tag"))
    return render_plain.datatags(tags)

def mymarkdown(text, limit=4):
    if text.count("\n") <= limit:
        return markdown.markdown(text)
    else:
        res="\n".join(text.split("\n")[0:limit])
        return markdown.markdown(res) + '&hellip;'

def get_motto():
    tot=model.get_totals()
    return str(tot.count) + " data sets, downloaded " + \
        str(tot.downloads) + " times."

tempglob = { 'dataformats': dataformats,
             'openid': web.webopenid,
             'getusername': model.get_username,
             'currenturl': get_current_url,
             'currentnotlogouturl': get_current_not_logout_url,
             'get_whatsnew': get_whatsnew,
             'get_datatags': get_datatags,
             'get_motto': get_motto,
             'mymarkdown': mymarkdown,
             'markdown': markdown.markdown,
             'type': type }
for name in url_helper.__all__:
    tempglob[name] = getattr(url_helper, name)

render = web.template.render('templates', base='base', globals=tempglob)
render_plain = web.template.render('templates', globals=tempglob)

class NotFound:
    """Handler class for URLs that do not encode a known resource."""

    def GET(self):
        """Responds with a HTTP error 404 for GET requests."""
        return web.notfound()

class Donate:
    """Handler class for the "Donate data" page."""

    donate_form=web.form.Form(
        web.form.Textbox("name", description="Your name:", id="focused"),
        web.form.Textbox("email", description="Your email address:"),
        web.form.Textbox("url", description="FTP or Web URL:"),
        web.form.Checkbox("directed", description="Directed:", 
                          value="True"),
        web.form.Checkbox("weighted", description="Weighted:",
                          value="True"),
        web.form.Checkbox("bipartite", description="Two-mode:",
                          value="True"),
        web.form.Checkbox("dynamic", description="Dynamic:", 
                          value="True"),
        web.form.Textbox("tags", description="Tags:"),
        web.form.Dropdown("licence", description="Licence:", args=[]),
        web.form.Textarea("description", 
                          description="Data format description:", 
                          cols=50, rows=10),
        web.form.Textarea("papers", description="Publication(s):", 
                          cols=50, rows=10),
        web.form.Button("Donate!", pre=recaptcha_text)
        )
    
    def GET(self):
        form=self.donate_form()
        lic=model.get_licences('id,name')
        form.licence.args=[(l.id,  l.name) for l in lic]        
        return render.donate(form, False, False, False)

    def POST(self):
        form=self.donate_form()
        if not form.validates():
            ## TODO
            None
        
        user_input=web.input()
        recaptcha_private_key = \
            file("../config/recaptcha_private_key").read().strip()
        valid=captcha.submit(user_input.recaptcha_challenge_field,
                             user_input.recaptcha_response_field,
                             recaptcha_private_key,
                             web.ctx.ip)

        if not valid.is_valid:
            return render.donate(form, True, False, False)

        try:
            web.sendmail('nexus@igraph.org', 'nexus@igraph.org', 
                         'Nexus donation', 
                         'Name:        ' + form.d.name        + '\n'   +
                         'Email:       ' + form.d.email       + '\n'   +
                         'URL:         ' + form.d.url         + '\n'   +
                         'Directed:    ' + str(form.d.directed) + '\n' + 
                         'Weighted:    ' + str(form.d.weighted) + '\n' +
                         'Bipartite:   ' + str(form.d.bipartite)+ '\n' + 
                         'Dynamic:     ' + str(form.d.dynamic)  + '\n' + 
                         'Tags:        ' + form.d.tags        + '\n'   + 
                         'Licence:     ' + form.d.licence     + '\n\n' + 
                         'Description: ' + form.d.description + '\n\n' +
                         'Publication: ' + form.d.papers      + '\n\n')
            return render.donate(form, True, True, True)
        except:
            return render.donate(form, True, True, False)

class About:
    """Renders the contents of the About page."""

    def GET(self):
        return render.about()

class Feedback:
    """Renders the feedback form."""

    feedback_form = web.form.Form(
        web.form.Textbox("name", description="Your name (optional):"),
        web.form.Textbox("email", description="Your email (optional):"), 
        web.form.Textarea("message", description="Your message:", 
                          cols=50, rows=10),
        web.form.Button("Send message", pre=recaptcha_text)
        )

    def GET(self):
        form=self.feedback_form()
        return render.feedback(form)

    def POST(self):
        form=self.feedback_form()
        if not form.validates():
            ## TODO
            pass

        user_input=web.input()
        recaptcha_private_key = \
            file("../config/recaptcha_private_key").read().strip()
        valid=captcha.submit(user_input.recaptcha_challenge_field,
                             user_input.recaptcha_response_field,
                             recaptcha_private_key,
                             web.ctx.ip)

        if not valid.is_valid:
            return render.feedback_ok(form, False, False)

        try:
            web.sendmail('nexus@igraph.org', 'nexus@igraph.org', 
                         'Nexus feedback', 
                         form.d.message + "\n\nEmail:" + form.d.email)
            return render.feedback_ok(form, True, True)
        except Exception, x:
            return render.feedback_ok(form, True, False)

class Index:
    """Renders the index page."""

    def list_datasets(self, user_input):
        format=user_input.format

        datasets=list(model.get_list_of_datasets(order=user_input.order))
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = list(model.get_tags(i))

        if format=='html':
            feed='/api/dataset_info?format=atom'
            return render.index(datasets, tags, "All Nexus data sets", feed)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_index(datasets, tags, 
                                          'All Nexus data sets')
        elif format=='text':
            for k, t in tags.iteritems():
                tags[k]=";".join(x.tag for x in t)
            web.header('Content-Type', 'text/plain')
            return render_plain.text_index(datasets, tags, 
                                           'All Nexus data sets')
        elif format=='rss':
            date=datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/rss+xml')
            return render_plain.rss_index(datasets, tags, 
                                          "All Nexus data sets", 
                                          date, web.ctx.homedomain, '')
        elif format=='atom':
            date=datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/atom+xml')
            return render_plain.atom_index(datasets, tags, 
                                           "All Nexus data sets", 
                                           date, web.ctx.homedomain, '')

    def format_text(self, dataset, tags, meta, papers):

        def format_attr(attr):
            return "Attribute: %s %s %s\n  %s" % \
                (attr.type, attr.datatype, attr.name, attr.description)

        tags=";".join([x.tag for x in tags])
        papers=[p.citation.replace("\n", "\n  ").strip() for p in papers]
        papers="\n  .\n".join(papers)
        main=dedent("""\
                Id: %i
                Name: %s
                Vertices: %s
                Edges: %s
                Tags: %s
                Date: %s
                Licence: %s
                Licence URL: %s
                Short description: %s
                Description: %s
                Citation: %s
        """) % (dataset.id, dataset.name, dataset.vertices, dataset.edges,
                tags, dataset.date[0:10], dataset.licence_name, 
                dataset.licence_url,
                dataset.shortdescription.replace("\n", "\n  ").strip(),
                dataset.description.replace("\n\n", "\n.\n").replace("\n", "\n  ").strip(), papers)
        meta="\n".join(format_attr(m) for m in meta)
        return main + meta

    def formats_for_dataset(self, id, filename):
        def good(base, f):
            fn=os.path.join("..", "data", id, base + dataformats[f])
            return os.path.exists(fn)
        base=list(model.get_dataset_file(id, 1))[0].filename
        form=[(f.name, f) for f in model.get_formats()]
        form=[ (n,f) for (n,f) in form if good(base, n) ]
        return dict(form)

    def dataset(self, user_input):
        id=user_input.id
        format=user_input.format
        dataset=[d for d in model.get_dataset(id)][0]
        if not dataset:
            return web.notfound()        
        tags=list(model.get_tags(dataset.id))
        meta=list(model.get_metadata(dataset.id))
        papers=model.get_papers(id)

        if format=='html':
            formats=self.formats_for_dataset(id, dataset.filename)
            return render.dataset(dataset, tags, meta, formats, papers)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_dataset(dataset, tags, meta, papers)
        elif format=='text':
            formatted=self.format_text(dataset, tags, meta, papers)
            web.header('Content-Type', 'text/plain')
            return render_plain.text_dataset(formatted)

    def list_tagged(self, user_input):
        tagname=user_input.tag.split("|")
        format=user_input.format
        if user_input.operator not in ("and", "or"):
            return web.notfound()
        datasets=list(model.get_list_tagged_as(tagname, user_input.operator))
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = list(model.get_tags(i))

        tagname=["'%s'" % t for t in tagname]
        tagname=(" " + user_input.operator + " ").join(tagname)

        if format=='html':
            feed='/api/dataset_info?format=atom&tag=%s&operator=%s' % \
                (user_input.tag, user_input.operator)
            return render.index(datasets, tags, 
                                "Data sets tagged %s" % tagname, feed)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_index(datasets, tags, 
                                          "Data sets tagged %s" 
                                          % tagname)
        elif format=='text':
            for k,t in tags.items():
                tags[k]=";".join([x.tag for x in t])
            web.header('Content-Type', 'text/plain')
            return render_plain.text_index(datasets, tags, 
                                           "Data sets tagged '%s'" 
                                           % tagname)
        elif format=="rss":
            date=datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/rss+xml')
            return render_plain.rss_index(datasets, tags, 
                                          "Nexus data sets tagged %s" 
                                          % tagname, date,
                                          web.ctx.homedomain, 'tagged/%s' 
                                          % tagname)
        elif format=="atom":
            date=datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/atom+xml')
            return render_plain.atom_index(datasets, tags, 
                                           "Nexus data sets tagged %s" 
                                           % tagname, date,
                                           web.ctx.homedomain, 'tagged/%s' 
                                           % tagname)        

    def GET(self):
        user_input=web.input(format="html", id=None, tag=None,
                             operator="or", order="date")
        if user_input.order not in ('date', 'name', 'popularity'):
            return web.notfound()
        if user_input.tag is not None and user_input.id is not None:
            return web.notfound()
        if user_input.id is not None:
            return self.dataset(user_input)
        elif user_input.tag is not None:
            return self.list_tagged(user_input)
        else:
            return self.list_datasets(user_input)

class Dataset:
    
    def GET(self):
        user_input=web.input(id=None, networkid=1)
        format=user_input.format
        if user_input.id is None:
            return web.notfound()
        
        datafile=model.get_dataset_file(user_input.id, user_input.networkid)
        if not datafile:
            return web.notfound()
        else:
            basename=[ d.filename for d in datafile][0]
            ext=dataformats[format]
        filename=os.path.join('..', 'data', user_input.id, basename + ext)
        try:
            f=open(filename)
            data=f.read()
            web.header('Content-Type', 'application/octet-stream')
            web.header('Content-Disposition', 
                       'attachment; filename="%s%s"' % (basename,ext))
            model.increase_downloads(user_input.id)
            return data
        except:
            return web.internalerror()

class Format:

    def format_one(self, format):
        return """Name: %s
Short description: %s
Description: %s
URL: %s""" % (format.name, format.shortdesc, 
              format.description.replace("\n", "\n  .\n").strip(),
              format.link)

    def format_text(self, formats):
        return "\n\n".join([ self.format_one(f) for f in formats ])

    def GET(self):
        user_input=web.input(format="html", dataformat=None)
        format=user_input.format
        dataformat=user_input.dataformat

        if dataformat:
            data=[d for d in model.get_format(dataformat)]
            if not data:
                return web.notfound
        else:
            data=[d for d in model.list_data_formats()]

        if format=='html':
            return render.format(data)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_formats(data)
        elif format=='text':
            formatted=self.format_text(data)
            web.header('Content-Type', 'text/plain')
            return render_plain.text_formats(formatted)

class MyForm(web.form.Form):

    def insert(self, where, *what):
        
        def find(f, seq):
            """Return first item in sequence where f(item) == True."""
            for idx, item in enumerate(seq):
                if f(item): 
                    return idx

        self.inputs=list(self.inputs)
        pos=find(lambda inp: inp.name==where, self.inputs)
        self.inputs[pos:pos] = what

class Spacer(web.form.Input):
    def get_type(self):
        return 'spacer'

    def render(self):
        return '<hr/>'
    
    def get_value(self):
        return None

class Label(web.form.Input):
    def get_type(self):
        return 'label'

    def render(self):
        return self.name

    def get_value(self):
        return None

add_form=MyForm(
    web.form.Textbox("name", description="Name:", id="focused", size=50),
    web.form.Textbox("shortdescription", description="Short description",
                     size=50),
    web.form.Textarea("description", description="Description:",
                      rows=10, cols=50),
    web.form.Textbox("tags", description="Tags:",
                     post="<br/> (comma separated)", size=50),
    web.form.Checkbox("directed", description="Directed:",
                      value="True"),
    web.form.Checkbox("weighted", description="Weighted:",
                      value="True"),
    web.form.Checkbox("bipartite", description="Two-mode:",
                      value="True"),
    web.form.Checkbox("dynamic", description="Dynamic:",
                      value="True"),
    web.form.Dropdown("licence", description="Licence:", args=[]),
    web.form.Textbox("vertices", description="Vertices:"),
    web.form.Textbox("edges", description="Edges:"),
    web.form.Textbox("filename", description="File name:", size=50),
    web.form.Textbox("source", description="Source:", size=50),
    web.form.Textarea("papers", description="Publication(s):", 
                      rows=5, cols=50),
    Spacer('spacer', description=''),
    web.form.Button("Add", description="Add")
    )

def add_meta_form(form, id):
    meta=model.get_metadata(id)
    for n, m in enumerate(meta):
        mkeep=web.form.Dropdown('meta%skeep' %n, 
                                description='Keep or delete:',
                                args=['keep', 'delete'])
        mtype=web.form.Dropdown("meta%stype" %n, 
                                description='Type:', 
                                args=['deleted', 'graph', 'vertex',
                                      'edge'])
        mtype.set_value(m.type)
        mdatatype=web.form.Dropdown('meta%sdatatype' %n,
                                    description='Data type:',
                                    args=['numeric', 'string'])
        mdatatype.set_value(m.datatype)
        mname=web.form.Textbox('meta%sname' %n,
                               description='Name:')
        mname.set_value(m.name)
        mdescription=web.form.Textarea('meta%sdescription' %n,
                                       description='Description:', 
                                       rows=5, cols=50)
        mdescription.set_value(m.description)
        form.insert("spacer", 
                    Spacer(''), 
                    Label('<h4>Meta data #%s</h4>' % str(int(n)+1), 
                          description=''),
                    mkeep, mtype, mdatatype, mname, mdescription)
    form.insert("spacer",
                Spacer(''),
                Label('<h4>Additional meta data</h4>', description=""),
                web.form.Dropdown('metaxadd', description='Add or not:',
                                  args=['do not add', 'add']),
                web.form.Dropdown("metaxtype", 
                                description='Type:', 
                                args=['graph', 'vertex', 'edge']),
                web.form.Dropdown('metaxdatatype',
                                    description='Data type:',
                                    args=['numeric', 'string']),
                web.form.Textbox('metaxname',
                               description='Name:'),
                web.form.Textarea('metaxdescription',
                                       description='Description:', 
                                       rows=5, cols=50))
    return form

def add_meta_form_add(form):
    for n in range(0,5):
        form.insert('spacer',
                    Spacer(''),
                    Label('<h4>Meta data #%s</h4>' % str(int(n)+1), 
                          description=""),
                    web.form.Dropdown('meta%sadd' % n, 
                                      description='Add:',
                                      args=['do not add', 'add']),
                    web.form.Dropdown("meta%stype" % n, 
                                      description='Type:', 
                                      args=['graph', 'vertex', 'edge']),
                    web.form.Dropdown('meta%sdatatype' % n,
                                      description='Data type:',
                                      args=['numeric', 'string']),
                    web.form.Textbox('meta%sname' % n,
                                     description='Name:'),
                    web.form.Textarea('meta%sdescription' % n,
                                      description='Description:', 
                                      rows=5, cols=50))
    return form

class Add:

    def GET(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_form()
        form=add_meta_form_add(form)
        lic=model.get_licences('id,name')
        form.licence.args=[(l.id,  l.name) for l in lic]
        return render.add(form, False, False, None)

    def POST(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_form()
        form=add_meta_form_add(form)
        if not form.validates():
            ## TODO
            None

        did=model.new_dataset(name=form.d.name, 
                              shortdescription=form.d.shortdescription,
                              description=form.d.description,
                              licence=int(form.d.licence),
                              source=form.d.source,
                              date=web.SQLLiteral("CURRENT_DATE"))

        model.new_network(dataset=did, id=1, description="",
                          vertices=form.d.vertices,
                          edges=form.d.edges,
                          filename=form.d.filename, 
                          date=web.SQLLiteral("CURRENT_DATE"))

        cits=form.d.papers.split("\r\n\r\n")
        for cit in cits:
            model.new_citation(dataset=did, citation=cit.strip())

        if form.d.tags.strip()=="":
            tags=[]
        else:
            tags=[t.strip() for t in form.d.tags.split(",")]
        for tag in tags:
            model.new_dataset_tag(dataset=did, tag=tag.strip())

        if form.d.directed and 'directed' not in tags:
            model.new_dataset_tag(dataset=did, tag="directed")
        elif 'undirected' not in tags:
            model.new_dataset_tag(dataset=did, tag="undirected")

        if form.d.weighted and 'weighted' not in tags:
            model.new_dataset_tag(dataset=did, tag="weighted")

        if form.d.bipartite and 'bipartite' not in tags:
            model.new_dataset_tag(dataset=did, tag="bipartite")
            
        if form.d.dynamic and 'dynamic' not in tags:
            model.new_dataset_tag(dataset=did, tag="dynamic")

        for m in range(0,5):
            if form.d['meta%sadd' % m] == "add":
                model.add_meta(id=did, type=form.d['meta%stype' % m],
                               name=form.d['meta%sname' % m],
                               datatype=form.d['meta%sdatatype' % m],
                               description=form.d['meta%sdescription' % m])
                               
        return render.add(form, True, False, did)

class Edit:

    def GET(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        
        form=add_form()
        form.Add.name='Submit'
        lic=model.get_licences('id,name')
        form.licence.args=[(l.id,  l.name) for l in lic]
        ds=model.get_dataset(id)
        ds=[d for d in ds][0]
        form.fill(ds)
        papers=[p.citation for p in model.get_papers(id)]
        form.papers.value="\n\n".join(papers)
        tags=set(t.tag for t in model.get_tags(id))
        form.directed.set_value('directed' in tags)
        form.weighted.set_value('weighted' in tags)
        form.bipartite.set_value('bipartite' in tags)
        form.dynamic.set_value('dynamic' in tags)
        form.directed.value=form.weighted.value=form.bipartite.value= \
            form.dynamic.value="True"
        tags -= set(("weighted", "bipartite", "directed", "undirected", 
                     "dynamic"))
        form.tags.value=", ".join(sorted(tags))
        
        form=add_meta_form(form, id)
        return render.add(form, False, True, int(id))

    def POST(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        
        form=add_form()
        form=add_meta_form(form, id)
        if not form.validates():
            ## TODO
            None

        if form.d.tags.strip()=="":
            tags=[]
        else:
            tags=[t.strip() for t in form.d.tags.split(',')]
        if form.d.weighted and 'weighted' not in tags:
            tags.append('weighted')
        if form.d.directed and 'directed' not in tags:
            tags.append('directed')
        if not form.d.directed and 'undirected' not in tags:
            tags.append('undirected')
        if form.d.bipartite and 'bipartite' not in tags:
            tags.append('bipartite')
        if form.d.dynamic and 'dynamic' not in tags:
            tags.append('dynamic')

        papers=[p.strip() for p in form.d.papers.split("\r\n\r\n")]

        model.update_dataset(id=id, name=form.d.name, 
                             shortdescription=form.d.shortdescription,
                             description=form.d.description,
                             tags=tags, licence=form.d.licence, 
                             vertices=form.d.vertices,
                             edges=form.d.edges,
                             filename=form.d.filename,
                             source=form.d.source,
                             papers=papers)

        m=0
        while ('meta%skeep' % m ) in form.d:
            if form.d['meta%skeep' % m] == "delete":
                model.delete_meta(id, type=form.d['meta%stype' % m], 
                                  name=form.d['meta%sname' % m])
            else:
                model.update_meta(id, type=form.d['meta%stype' % m],
                                  name=form.d['meta%sname' % m],
                                  datatype=form.d['meta%sdatatype' % m],
                                  description=form.d['meta%sdescription' % m])
            m=m+1

        if form.d['metaxadd']=='add':
            ## TODO: check that it does not exist already
            model.add_meta(id=id, type=form.d['metaxtype'], 
                           name=form.d['metaxname'],
                           datatype=form.d['metaxdatatype'],
                           description=form.d['metaxdescription'])
    
        return render.add(form, True, True, int(id))
    
add_licence_form=web.form.Form(
    web.form.Textbox("name", description="Name:", id="focused", size=50),
    web.form.Textarea("text", description="Text:", cols=50, rows=2),
    web.form.Textarea("fulltext", description="Full text:", cols=50, 
                      rows=10),
    web.form.Textbox("link", description="URL:", size=50),
    web.form.Button("Add")
    )

class AddLicence:
    
    def GET(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        form=add_licence_form()
        return render.addlicence(form, False, False, False)

    def POST(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_licence_form()
        if not form.validates():
            # TODO
            None
            
        lid=model.new_licence(name=form.d.name, text=form.d.text, 
                              fulltext=form.d.fulltext, link=form.d.link)

        return render.addlicence(form, True, False, lid)

class EditLicence:
    def GET(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_licence_form()
        form.Add.name='Submit'
        lic=model.get_licence(id)
        lic=list(lic)[0]
        form.fill(lic)
        return render.addlicence(form, False, False, False)

    def POST(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        
        form=add_licence_form()
        if not form.validates():
            ## TODO
            pass

        model.update_licence(id=id, name=form.d.name, text=form.d.text,
                             fulltext=form.d.fulltext, link=form.d.link)

        return render.addlicence(form, True, True, int(id))

class Logout:
    def GET(self):
        web.webopenid.logout()
        return render.logout()

    def POST(self):
        user_input=web.input()
        web.webopenid.logout()
        return render.logout()

class OpenID(web.webopenid.host):
    def POST(self):
        i=web.input(return_to="/")
        user=model.get_username(i.openid)
        if user:
            web.webopenid.host.POST(self)
        else:
            return web.seeother("/loginfailed")

class LoginFailed:
    def GET(self):
        return render.loginfailed()

class Admin:
    def GET(self):
        return render.admin()

def run_r(cmd):
    Rscript=file("../config/Rscript").read().strip()
    args=[Rscript, "-e", cmd]
    p=subprocess.Popen(args, stdout=subprocess.PIPE)
    ret=p.wait()
    out=p.stdout.read()
    return [ret, out]

class Check:

    def check_r_igraph(self, ds, filename, tags, meta):

        res=odict.odict()
        
        ## File exists
        ex=os.path.exists(filename)
        res['Data file exists'] = ex
        if not ex:
            return res

        ## File can be loaded
        loadcode = 'library(igraph) ; load("%s") ;' %filename
        try:
            run_r(loadcode)
            res['Data file can be loaded'] = True
        except:
            res['Data file can be loaded'] = False
            return res
        
        ## File contains a single igraph graph
        code = loadcode + 'g <- get(ls()[1]) ; cat(is.igraph(g), " ") ;'
        try:
            ret, out=run_r(code)
            if out[0:4] == 'TRUE':
                res['File contains proper data'] = True
            else:
                res['File contains proper data'] = False
                return res
        except:
            res['File contains proper data'] = False
            return res

        ## Number of vertices and edges 
        code = code + 'cat(vcount(g), ecount(g), is.directed(g), ' + \
            'is.weighted(g), is.bipartite(g), " ") ; '
        ret, out=run_r(code)
        nm=out.split()[1:]
        res['Number of vertices'] = int(nm[0])==ds.vertices
        res['Number of edges'] = int(nm[1])==ds.edges

        ## Tags
        tags=[t.tag for t in tags]
        if 'directed' in tags:
            res['Tags, directed'] = nm[2]=='TRUE'
        if 'undirected' in tags:
            res['Tags, undirected'] = nm[2]=='FALSE'
        if 'weighted' in tags:
            res['Tags, weighted'] = nm[3]=='TRUE'
        if 'bipartite' in tags:
            res['Tags, bipartite'] = nm[4]=='TRUE'

        ## Metadata
        code=loadcode + 'g <- get(ls()[1]) ; ' + \
            'cat(list.vertex.attributes(g), sep="", "\\n") ;' + \
            'cat(list.edge.attributes(g), sep="", "\\n") ;'
        ret, out=run_r(code)
        attr=out.split("\n")[0:2]
        for m in meta:
            if m.type == "vertex":
                res["Metadata, vertex, '%s'" % m.name]=m.name in attr[0]
            elif m.type == "edge":
                res["Metadata, edge, '%s'" % m.name]=m.name in attr[1]

        return res

    def check_graphml(self, ds, filename, tags, meta):
        res=odict.odict()

        ## File exists
        ex=os.path.exists(filename)
        res['Data file exists'] = ex
        if not ex:
            return res

        ## File can be loaded
        try:
            g=igraph.Graph.Read_GraphML(filename)
            res['Data file can be loaded'] = True
        except Exception, x:
            print str(x)
            res['Data file can be loaded'] = False
            return res

        ## Number of vertices and edges
        res['Number of vertices'] = ds.vertices == g.vcount()
        res['Number of edges'] = ds.edges == g.ecount()

        vattr=g.vs.attribute_names()
        eattr=g.es.attribute_names()

        ## Tags 
        tags=[t.tag for t in tags]
        if 'directed' in tags:
            res['Tags, directed'] = g.is_directed()
        if 'undirected' in tags:
            res['Tags, undirected'] = not g.is_directed()
        if 'weighted' in tags:
            res['Tags, weighted'] = 'weight' in eattr
        if 'bipartite' in tags:
            res['Tags, bipartite'] = 'type' in vattr
            
        ## Metadata
        for m in meta:
            if m.type == "vertex":
                res["Metadata, vertex, '%s'" % m.name] = m.name in vattr
            elif m.type == "edge":
                res["Metadata, edge, '%s'" % m.name] = m.name in eattr

        return res

    def check_pajek(self, ds, filename, tags, meta):
        res=odict.odict()

        ## File exists
        ex=os.path.exists(filename)
        res['Data file exists'] = ex
        if not ex:
            return res

        ## File can be loaded
        try:
            g=igraph.Graph.Read_Pajek(filename)
            res['Data file can be loaded'] = True
        except Exception, x:
            print str(x)
            res['Data file can be loaded'] = False
            return res

        ## Number of vertices and edges 
        res['Number of vertices'] = ds.vertices == g.vcount()
        res['Number of edges'] = ds.edges == g.ecount()

        vattr=g.vs.attribute_names()
        eattr=g.es.attribute_names()

        ## Tags 
        tags=[t.tag for t in tags]
        if 'directed' in tags:
            res['Tags, directed'] = g.is_directed()
        if 'undirected' in tags:
            res['Tags, undirected'] = not g.is_directed()
        if 'weighted' in tags:
            res['Tags, weighted'] = 'weight' in eattr

        pajek_meta=(('edge', 'weight'), ('vertex', 'id') )
        meta=[ m for m in meta if (m.type, m.name) in pajek_meta ]
        for m in meta:
            if m.type == "vertex":
                res["Metadata, vertex, '%s'" % m.name] = m.name in vattr
            elif m.type == "edge":
                res["Metadata, edge, '%s'" % m.name] = m.name in eattr        
            
        return res

    def check_excel(self, ds, filename, tags, meta):
        res=odict.odict()

        ## File exists
        ex=os.path.exists(filename)
        res['Data file exists'] = ex
        if not ex:
            return res

        ## File can be loaded
        try:
            wb=xlrd.open_workbook(filename)
            res['Data file can be loaded'] = True
        except Exception, x:
            res['Data file can be loaded'] = False
            return res

        ## Number of sheets
        res['Number of sheets'] = wb.nsheets in (2,3)        

        ## Sheet names
        shnames=wb.sheet_names()        
        res['Sheet names'] = shnames[0] == 'Network' and \
            shnames[1] == 'Edges' and shnames[2] == 'Vertices'

        ## Number of vertices and edges
        sh1=wb.sheet_by_index(0)
        res['Number of vertices'] = int(sh1.cell_value(1,1)) == ds.vertices
        res['Number of edges'] = int(sh1.cell_value(2,1)) == ds.edges

        ## Tags
        mytags=sh1.row(4)[1:]
        mytags=[ t.value for t in mytags ]
        tags=[ t.tag for t in tags ]
        res['Tags'] = len(set(mytags) & set(tags)) == len(mytags)

        ## Metadata
        sh2=wb.sheet_by_index(1)
        eattr=[ c.value for c in sh2.row(0)[2:] ]
        if wb.nsheets == 3: 
            sh3=wb.sheet_by_index(2)
            vattr=[ c.value for c in sh3.row(0) ]
        else:
            vattr=[]
        for m in meta:
            if m.type == "vertex":
                res["Metadata, vertex, '%s'" % m.name] = m.name in vattr
            elif m.type == "edge":
                res["Metadata, edge, '%s'" % m.name] = m.name in eattr

        return res

    def check(self, ds, format, tags, meta, filename):
        if format=="R-igraph":
            return self.check_r_igraph(ds, filename, tags, meta)
        elif format=="GraphML":
            return self.check_graphml(ds, filename, tags, meta)
        elif format=="Pajek":
            return self.check_pajek(ds, filename, tags, meta)
        elif format=="Excel":
            return self.check_excel(ds, filename, tags, meta)
        return None

    def GET(self, id):
        ds=model.get_dataset(id)
        ds=[d for d in ds][0]
        checkres=odict.odict()
        for k,v in dataformats.items():
            tags=model.get_tags(ds.id)
            meta=model.get_metadata(ds.id)
            try: 
                checkres[k] = self.check(ds, k, tags, meta,
                                         "../data/" + str(id) + "/" + 
                                         ds.filename + v)
            except Exception, x:
                checkres[k] = { "Cannot run check": str(x) }
                
        return render.check(id, ds, checkres)      

class Licence:

    def format_text(self, licences):
        
        def format_one(lic):
            return dedent('''\
                     Id: %s
                     Name: %s
                     Short Description: %s
                     URL: %s''' % (lic.id, lic.name, lic.text, lic.link))
        
        return "\n".join(format_one(l) for l in licences)

    def GET(self):
        user_input=web.input(id=None, format="html")
        format=user_input.format

        if user_input.id is None:
            lic=model.get_licences("*")
        else:
            lic=model.get_licence(user_input.id)
        
        if format=='html':
            return render.licence(lic)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_licence(lic)
        elif format=='text':
            web.header('Content-Type', 'text/plain')
            return self.format_text(lic)

class Blog:
    
    def GET(self):
        user_input=web.input(id=None, admin=False)
        unpublished = web.webopenid.status() in admins_openid and \
            user_input.admin
        if user_input.id is None:
            entries=model.get_blog(unpublished=unpublished)
        else:
            entries=model.get_blog(ids=[user_input.id], 
                                   unpublished=unpublished)
        return render.blog(entries, user_input.id)

add_blog_form=web.form.Form(
    web.form.Textbox('title', description='Title:', id='focused', size=50),
    web.form.Textbox('date', description='Date:', size=50),
    web.form.Textarea('entry', description='Text:', rows=20, cols=50),
    web.form.Checkbox('published', description="Published?", value="True"),
    web.form.Button("Add")
    )

class AddBlog:
    
    def GET(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_blog_form()
        return render.addblog(form, None, added=False, edit=False)

    def POST(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        form=add_blog_form()
        if not form.validates():
            pass
        
        bid=model.new_blog_entry(title=form.d.title, 
                                 date=form.d.date,
                                 entry=form.d.entry,
                                 published=form.d.published)
        
        return render.addblog(form, bid, added=True, edit=False)

class EditBlog:
    
    def GET(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        form=add_blog_form()        
        form.Add.name='Submit'
        blog=model.get_blog(ids=[id], unpublished=True)
        form.fill(blog[0])
        return render.addblog(form, id, added=False, edit=True)

    def POST(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        form=add_blog_form()
        if not form.validates():
            # TODO
            pass
        
        model.update_blog(id=id, title=form.d.title, date=form.d.date,
                          entry=form.d.entry, published=form.d.published)
        
        return render.addblog(form, id, added=True, edit=True)

class Docs:
    
    def GET(self):
        import docs
        cc=markdown.markdown(docs.docs)
        return render.docs(cc)

class Delete:

    def GET(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        model.delete_dataset(int(id))
        return render.delete()

class Recreate:

    def to_excel(self, inputfile, id, ds):
        g=igraph.Graph.Read_GraphML(inputfile + ".graphml")
        wb=xlwt.Workbook()

        ## The network
        tags=list(model.get_tags(id))
        papers=list(model.get_papers(id))
        sh1=wb.add_sheet('Network')
        sh1.write(0, 0, 'Name')
        sh1.write(0, 1, ds.name)
        sh1.write(1, 0, 'Vertices')
        sh1.write(1, 1, ds.vertices)
        sh1.write(2, 0, 'Edges')
        sh1.write(2, 1, ds.edges)
        sh1.write(3, 0, 'Type')
        if 'directed' in tags:
            sh1.write(3, 1, 'directed')
        else:
            sh1.write(3, 1, 'undirected')
        sh1.write(4, 0, 'Tags')
        for n, t in enumerate(tags):
            sh1.write(4, n+1, t.tag)
        sh1.write(5, 0, 'Licence')
        sh1.write(5, 1, ds.licence_name)
        sh1.write(6, 0, 'Licence URL')
        sh1.write(6, 1, ds.licence_url)
        sh1.write(7, 0, 'Short description')
        sh1.write(7, 1, ds.shortdescription)
        sh1.write(8, 0, 'Description')
        sh1.write(8, 1, ds.description)
        for n, p in enumerate(papers):
            sh1.write(9+n, 0, 'Paper')
            sh1.write(9+n, 1, p.citation)

        ## The edges
        el=g.get_edgelist()
        if 'name' in g.vs.attribute_names():
            n=g.vs['name']
            el=[ (n[e[0]], n[e[1]]) for e in el]
        sh2=wb.add_sheet('Edges')
        sh2.write(0, 0, 'From')
        sh2.write(0, 1, 'To')
        for n, e in enumerate(el):
            sh2.write(n+1, 0, e[0])
            sh2.write(n+1, 1, e[1])

        ## Edge attributes
        ea=g.es.attribute_names()
        for n, attr in enumerate(ea):
            sh2.write(0, n+2, attr)
            for r, v in enumerate(g.es[attr]):
                sh2.write(r+1, n+2, v)

        ## Vertex attributes
        va=g.vs.attribute_names()
        if len(va) != 0:
            sh3=wb.add_sheet('Vertices')
            for n, attr in enumerate(va):
                sh3.write(0, n, attr)
                for r, v in enumerate(g.vs[attr]):
                    sh3.write(r+1, n, v)
            
        wb.save(inputfile + ".xls")
    
    def GET(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        
        inputfile = list(model.get_dataset_file(id, 1))[0].filename
        inputfile = "../data/%s/%s" % (id, os.path.basename(inputfile))
        rcode = dedent('library(igraph) ; \
                  load("%s.Rdata") ; \
                  g <- get(ls()[1]) ; \
                  write.graph(g, file="%s.graphml", format="graphml") ; \
                  write.graph(g, file="%s.net", format="pajek") ; \
                  ' % ((inputfile,) * 3))
        ret, out=run_r(rcode)

        ds=list(model.get_dataset(id))[0]
        if ds.vertices < 2**26 and ds.edges < 2**16:
            self.to_excel(inputfile, id, ds)

        return render.recreate(id)
        
app = web.application(urls, globals())

web.webopenid.sessions = \
    web.session.Session(app, web.session.DiskStore('../sessions'),
                        initializer={})
web.webopenid.store=openid.store.filestore.FileOpenIDStore('../sessions')

if __name__ == '__main__':
    try:
        app.run()
    except Exception, x:
        web.header('Content-Type', 'text/plain')
        print str(x)
