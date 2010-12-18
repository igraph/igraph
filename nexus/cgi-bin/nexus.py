#! /usr/bin/python2.6

import sys
sys.path.append("../python")

import web
import model
from recaptcha.client import captcha
import math
import random

web.config.debug=True

urls = (
    '/?',                                  'About',
    '/openid',                             'web.webopenid.host',
    '/login',                              'Login',
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

formats = ('html', 'xml', 'text')
dataformats = { 'R-igraph': '.Rdata' }

tempglob = { 'whatsnew': 'Nothing',
             'datatags': 'None',
             'dataformats': dataformats,
             'openid': web.webopenid      }

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

class Base:
    def __init__(self): 
        global tempglob
        w=model.whatsnew()
        t=[tag for tag in model.datatags()]
        print t
        c=[tag.count for tag in model.datatags()]
        maxc=max(c)
        c=[ int(math.log(cc+1) / math.log(maxc+1) * 5) for cc in c]
        for i in range(len(c)):
            t[i].count=c[i]
        random.shuffle(t)
        tempglob['whatsnew']=render_plain.whatsnew(w)
        tempglob['datatags']=render_plain.datatags(t)

class NotFound:

    def GET(self):
        return web.notfound()

class Home(Base):
    
    def GET(self):
        return render.home()

class AdvancedSearch(Base):
    
    def GET(self):
        ## TODO
        return render.home()
    
class Donate(Base):

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

class Blog(Base):
    
    def GET(self):
        return render.blog()

class About(Base):
    
    def GET(self):
        return render.about()

class Feedback(Base):

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

class Index(Base):
    
    @knownformat
    def GET(self, format='html'):
        datasets=model.get_list_of_datasets()
        datasets=[d for d in datasets]
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:            
            tags[i] = [t for t in model.get_tags(i)]
        if format=='html':
            return render.index(datasets, tags, "All data sets")
        elif format=='xml':
            web.header('Content-Type', 'text/xml')
            return render_plain.xml_index(datasets, tags, 
                                          'All data sets')
        elif format=='text':
            for k,t in tags.items():
                tags[k]=";".join([x.tag for x in t])
            web.header('Content-Type', 'text/plain')
            return render_plain.text_index(datasets, tags, 
                                           'All data sets')

class Dataset(Base):

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

class Tagged(Base):

    @knownformat
    def GET(self, format, tagname=None):
        datasets=model.get_list_tagged_as(tagname)
        datasets=[d for d in datasets]
        ids=[d.id for d in datasets]
        tags={}
        for i in ids:
            tags[i] = [t for t in model.get_tags(i)]

        if format=='html':
            return render.index(datasets, tags, 
                                "Data sets tagged '%s'" % tagname)
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
            
        
class GetData(Base):
    
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

class Format(Base):

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

class Add(Base):

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

class Edit(Base):

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
        print form.d
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
    
class AddLicence(Base):
    
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

class Login(Base):
    def GET(self):
        return render.login(web.webopenid.form('/openid'))

class Admin(Base):
    def GET(self):
        return render.admin()
        
app = web.application(urls, globals())

if __name__ == '__main__':
    app.run()
