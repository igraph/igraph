#! /usr/bin/python2.6
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

from datetime import datetime
from functools import wraps
from itertools import izip
from operator import attrgetter
from recaptcha.client import captcha
from textwrap import dedent

web.config.debug = True

# URL mappings used in the Nexus web application
urls = (
    '/?',                                  'About',
    '/openid',                             'OpenID',
    '/loginfailed',                        'LoginFailed',
    '/logout',                             'Logout',
    '/admin',                              'Admin',
#    '/advanced',                           'AdvancedSearch',
    '/donate',                             'Donate',
#    '/blog',                               'Blog',
    '/about',                              'About',
    '/feedback',                           'Feedback',
    '/addlicence',                         'AddLicence',
    '/add',                                'Add',
    '/edit/(\d+)',                         'Edit',
    '/(\w+)/?',                            'Index',
    '/(\w+)/dataset/?',                      'Index',
    '/(\w+)/dataset/(\d+)',                'Dataset',
    '/([\w-]+)/getdata/(\d+)(?:/(\d+))?',  'GetData',
    '/(\w+)/tagged/(\w+)',                 'Tagged',
    '/(\w+)/format/?',                     'Format',
    '/(\w+)/format/([\w-]+)',              'Format',
    '.*',                                  'NotFound'
    )

# reCAPTCHA keys
recaptcha_pubkey = "6Lfqjb8SAAAAAJyGZQrvqgs7JWjpNp_Vf9dpTMxy"
recaptcha_private_key = "6Lfqjb8SAAAAAO_ElXNZyzVXbP5xffMs6IVypJbB"

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
dataformats = { 'R-igraph': '.Rdata' }

def get_current_url():
    """Returns the URL of the current page being produced"""
    return web.ctx.fullpath

def get_current_not_logout_url():
    """Returns the URL of the current page being produced, except for the
    logout page.

    For the logout and login failure pages, returns the root page."""
    fp=web.ctx.fullpath
    if fp in ('/logout', '/loginfailed'):
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


tempglob = { 'whatsnew': 'Nothing',
             'datatags': 'None',
             'dataformats': dataformats,
             'openid': web.webopenid,
             'getusername': model.get_username,
             'currenturl': get_current_url,
             'currentnotlogouturl': get_current_not_logout_url,
             'get_whatsnew': get_whatsnew,
             'get_datatags': get_datatags}
for name in url_helper.__all__:
    tempglob[name] = getattr(url_helper, name)

render = web.template.render('templates', base='base', globals=tempglob)
render_plain = web.template.render('templates', globals=tempglob)

def knownformat(fn):
    """Decorator that checks whether the format argument of a function is
    one of the known page formats."""
    @wraps(fn)
    def new(*args, **kwds):
        if args[1] not in formats:
            return web.badrequest()
        return fn(*args)
    return new

def knowndataformat(fn):
    """Decorator that checks whether the format argument of a function is
    one of the known data formats."""
    @wraps(fn)
    def new(*args, **kwds):
        if args[1] not in dataformats:
            return web.badrequest()
        return fn(*args, **kwds)
    return new

class NotFound:
    """Handler class for URLs that do not encode a known resource."""

    def GET(self):
        """Responds with a HTTP error 404 for GET requests."""
        return web.notfound()

class Home:
    """Handler class for the home URL."""

    def GET(self):
        """Renders the homepage."""
        return render.home()

class AdvancedSearch:
    """Handler class for the advanced search page."""

    def GET(self):
        """Yet to be implemented."""
        return render.home()
    
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
        valid=captcha.submit(user_input.recaptcha_challenge_field,
                             user_input.recaptcha_response_field,
                             recaptcha_private_key,
                             web.ctx.ip)

        if not valid.is_valid:
            return render.donate(form, True, False, False)

        web.config.smtp_server = 'smtp.gmail.com'
        web.config.smtp_port = 587
        web.config.smtp_username = 'nexus.repository@gmail.com'
        web.config.smtp_password = 'bhu8nji9'
        web.config.smtp_starttls = True
        try:
            web.sendmail(form.d.name, 'nexus.repository@gmail.com', 
                         'Donation', 
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

class Blog:
    
    def GET(self):
        return render.blog()

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
        valid=captcha.submit(user_input.recaptcha_challenge_field,
                             user_input.recaptcha_response_field,
                             recaptcha_private_key,
                             web.ctx.ip)

        if not valid.is_valid:
            return render.feedback_ok(form, False, False)

        web.config.smtp_server = 'smtp.gmail.com'
        web.config.smtp_port = 587
        web.config.smtp_username = 'nexus.repository@gmail.com'
        web.config.smtp_password = 'bhu8nji9'
        web.config.smtp_starttls = True
        try:
            web.sendmail(form.d.name, 'nexus.repository@gmail.com', 
                         'Feedback', 
                         form.d.message + "\n\nEmail:" + form.d.email)
            return render.feedback_ok(form, True, True)
        except:
            return render.feedback_ok(form, True, False)

class Index:
    """Renders the index page."""

    @knownformat
    def GET(self, format='html'):
        datasets=list(model.get_list_of_datasets())
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = list(model.get_tags(i))

        if format=='html':
            feed='/atom'
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
            

class Dataset:
    """Renders the page of a dataset."""

    def format_text(self, dataset, tags, papers):
        tags=";".join([x.tag for x in tags])
        papers=[p.citation.replace("\n", "\n  ").strip() for p in papers]
        papers="\n  .\n".join(papers)
        return dedent("""\
                Id: %i
                Name: %s
                Vertices: %s
                Edges: %s
                Tags: %s
                Date: %s
                Licence: %s
                Description: %s
                Citation: %s
        """ % (dataset.id, dataset.name, dataset.vertices, dataset.edges,
               tags, dataset.date, dataset.licence, 
               dataset.description.replace("\n", "\n  ").strip(), papers))

    @knownformat
    def GET(self, format, id):
        dataset=[d for d in model.get_dataset(id)][0]
        if not dataset:
            return web.notfound()
        tags=list(model.get_tags(dataset.id))
        formats=dict((f.name, f) for f in model.get_formats())
        papers=model.get_papers(id)

        if format=='html':
            return render.dataset(dataset, tags, formats, papers)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_dataset(dataset, tags, papers)
        elif format=='text':
            formatted=self.format_text(dataset, tags, papers)
            web.header('Content-Type', 'text/plain')
            return render_plain.text_dataset(formatted)

class Tagged:
    @knownformat
    def GET(self, format, tagname=None):
        datasets=list(model.get_list_tagged_as(tagname))
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = list(model.get_tags(i))

        if format=='html':
            feed='/atom/tagged/%s' % tagname
            return render.index(datasets, tags, 
                                "Data sets tagged '%s'" % tagname, feed)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_index(datasets, tags, 
                                          "Data sets tagged '%s'" 
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
            
            
        
class GetData:
    
    @knowndataformat
    def GET(self, format, id, nid):
        if nid is None:
            nid=1
        datafile=model.get_dataset_file(id, nid)
        if not datafile:
            return web.notfound()
        else:
            basename=[ d.filename for d in datafile][0]
            ext=dataformats[format]
        filename=os.path.join('..', 'data', id, basename + ext)
        try:
            f=open(filename)
            data=f.read()
            web.header('Content-Type', 'application/octet-stream')
            web.header('Content-Disposition', 
                       'attachment; filename="%s%s"' % (basename,ext))
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

    @knownformat
    def GET(self, format, dataformat=None):
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

add_form=web.form.Form(
    web.form.Textbox("name", description="Name:", id="focused", size=50),
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
    web.form.Button("Add")
    )

class Add:

    def GET(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_form()
        lic=model.get_licences('id,name')
        form.licence.args=[(l.id,  l.name) for l in lic]
        return render.add(form, False, False, None)

    def POST(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=add_form()
        if not form.validates():
            ## TODO
            None

        did=model.new_dataset(name=form.d.name, 
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
        tags -= ("weighted", "bipartite", "directed", "undirected", "dynamic")
        form.tags.value=", ".join(sorted(tags))
        return render.add(form, False, True, id)

    def POST(self, id):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        
        form=add_form()
        if not form.validates():
            ## TODO
            None

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
                             description=form.d.description,
                             tags=tags, licence=form.d.licence, 
                             vertices=form.d.vertices,
                             edges=form.d.edges,
                             filename=form.d.filename,
                             source=form.d.source,
                             papers=papers)
    
        return render.add(form, True, True, id)
    
class AddLicence:
    
    add_licence_form=web.form.Form(
        web.form.Textbox("name", description="Name:", id="focused", size=50),
        web.form.Textarea("text", description="Text:", cols=50, rows=2),
        web.form.Textarea("fulltext", description="Full text:", cols=50, 
                          rows=10),
        web.form.Textbox("link", description="URL:", size=50),
        web.form.Button("Add")
        )

    def GET(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")
        form=self.add_licence_form()
        return render.addlicence(form, False)

    def POST(self):
        if web.webopenid.status() not in admins_openid:
            return web.seeother("/login")

        form=self.add_licence_form()
        if not form.validates():
            # TODO
            None
            
        lid=model.new_licence(name=form.d.name, text=form.d.text, 
                              fulltext=form.d.fulltext, link=form.d.link)

        return render.addlicence(form, True)

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
        
app = web.application(urls, globals())
web.webopenid.sessions = \
    web.session.Session(app, web.session.DiskStore('../sessions'),
                        initializer={})
web.webopenid.store=openid.store.filestore.FileOpenIDStore('../sessions')

if __name__ == '__main__':
    app.run()
