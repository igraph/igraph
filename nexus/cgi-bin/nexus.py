#! /usr/bin/python2.6

import sys
sys.path.append("../python")

import web
import model
from recaptcha.client import captcha
import math
import random
import openid.store.filestore
import datetime

web.config.debug = True

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
    '/(\w+)/dataset/(\d+)',                'Dataset',
    '/([\w-]+)/getdata/(\d+)(?:/(\d+))?',  'GetData',
    '/(\w+)/tagged/(\w+)',                 'Tagged',
    '/(\w+)/format/?',                     'Format',
    '/(\w+)/format/([\w-]+)',              'Format',
    '.*',                                  'NotFound'
    )

recaptcha_pubkey = "6Lfqjb8SAAAAAJyGZQrvqgs7JWjpNp_Vf9dpTMxy"
recaptcha_private_key = "6Lfqjb8SAAAAAO_ElXNZyzVXbP5xffMs6IVypJbB"

admins_openid=('https://launchpad.net/~gabor.csardi',
               'https://launchpad.net/~ntamas')

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

formats = ('html', 'xml', 'text', 'rss')
dataformats = { 'R-igraph': '.Rdata' }

def get_current_url():
    return web.ctx.fullpath

def get_current_not_logout_url():
    fp=web.ctx.fullpath
    if fp in ('/logout', '/loginfailed'):
        fp='/'
    return fp

def get_whatsnew():
    w=model.whatsnew()
    return render_plain.whatsnew(w)

def get_datatags():
        t=[tag for tag in model.datatags()]
        c=[tag.count for tag in model.datatags()]
        maxc=max(c)
        c=[ int(math.log(cc+1) / math.log(maxc+1) * 5) for cc in c]
        for i in range(len(c)):
            t[i].count=c[i]
        random.shuffle(t)
        return render_plain.datatags(t)


tempglob = { 'whatsnew': 'Nothing',
             'datatags': 'None',
             'dataformats': dataformats,
             'openid': web.webopenid,
             'getusername': model.get_username,
             'currenturl': get_current_url,
             'currentnotlogouturl': get_current_not_logout_url,
             'get_whatsnew': get_whatsnew,
             'get_datatags': get_datatags }

render = web.template.render('templates', base='base', globals=tempglob)
render_plain = web.template.render('templates', globals=tempglob)

def knownformat(fn):
    def new(*args):
        if args[1] not in formats:
            return web.badrequest()
        return fn(*args)
    return new

def knowndataformat(fn):
    def new(*args):
        if args[1] not in dataformats:
            return web.badrequest()
        return fn(*args)
    return new

class NotFound:

    def GET(self):
        return web.notfound()

class Home:
    
    def GET(self):
        return render.home()

class AdvancedSearch:
    
    def GET(self):
        ## TODO
        return render.home()
    
class Donate:

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
    
    def GET(self):
        return render.about()

class Feedback:

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
            None

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
    
    @knownformat
    def GET(self, format='html'):
        datasets=model.get_list_of_datasets()
        datasets=[d for d in datasets]
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:            
            tags[i] = [t for t in model.get_tags(i)]

        if format=='html':
            feed='/rss'
            return render.index(datasets, tags, "All Nexus data sets", feed)
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_index(datasets, tags, 
                                          'All Nexus data sets')
        elif format=='text':
            for k,t in tags.items():
                tags[k]=";".join([x.tag for x in t])
            web.header('Content-Type', 'text/plain')
            return render_plain.text_index(datasets, tags, 
                                           'All Nexus data sets')
        elif format=='rss':
            date=datetime.datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/rss+xml')
            return render_plain.rss_index(datasets, tags, 
                                          "All Nexus data sets", 
                                          date, web.ctx.homedomain, '')

class Dataset:

    def format_text(self, dataset, tags, papers):
        tags=";".join([x.tag for x in tags])
        papers=[p.citation.replace("\n", "\n  ").strip() for p in papers]
        papers="\n  .\n".join(papers)
        return """Id: %i
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
       dataset.description.replace("\n", "\n  ").strip(), papers)

    @knownformat
    def GET(self, format, id):
        dataset=[d for d in model.get_dataset(id)][0]
        if not dataset:
            return web.notfound()
        tags=[t for t in model.get_tags(dataset.id)]
        formats={}
        for f in model.get_formats():
            formats[f.name] = f
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
        datasets=model.get_list_tagged_as(tagname)
        datasets=[d for d in datasets]
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = [t for t in model.get_tags(i)]

        if format=='html':
            feed='/rss/tagged/%s' % tagname
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
            date=datetime.datetime.today().strftime("%a, %d %b %Y %H:%M:%S +0200")
            web.header('Content-Type', 'application/rss+xml')
            return render_plain.rss_index(datasets, tags, 
                                          "Nexus data sets tagged %s" 
                                          % tagname, date,
                                          web.ctx.homedomain, 'tagged/%s' 
                                          % tagname)
            
            
        
class GetData:
    
    @knowndataformat
    def GET(self, format, id, nid):
        if nid==None: nid=1
        datafile=model.get_dataset_file(id, nid)
        if not datafile:
            return web.notfound()
        else:
            basename=[ d.filename for d in datafile][0]
            ext=dataformats[format]
        filename='../data/' + id + '/' + basename + ext
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
        tags=[t.tag for t in model.get_tags(id)]
        form.directed.set_value('directed' in tags)
        form.weighted.set_value('weighted' in tags)
        form.bipartite.set_value('bipartite' in tags)
        form.dynamic.set_value('dynamic' in tags)
        form.directed.value=form.weighted.value=form.bipartite.value= \
            form.dynamic.value="True"
        tags=[t for t in tags if t not in ("weighted", "bipartite", 
                                           "directed", "undirected", 
                                           "dynamic")]
        form.tags.value=", ".join(tags)
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
