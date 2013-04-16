#!/usr/bin/python
__author__ = 'xia'
from jinja2 import Template,Environment,PackageLoader
import cgi,cgitb
cgitb.enable()
form=cgi.FieldStorage(keep_blank_values=True)
env=Environment(loader=PackageLoader(__name__,"templates"))
print "Content-type: text/html\n\n"
if not form:
    template=env.get_template("index.html")
    print template.render()
if form.has_key("goto_modification1"):
    template=env.get_template("modification.html")
    print template.render()
elif form.has_key("goto_report"):
    template=env.get_template("report.html")
    print template.render()
elif form.has_key("goto_muta"):
    template=env.get_template("mutagenesis.html")
    print template.render()



