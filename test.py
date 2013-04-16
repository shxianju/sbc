#!/usr/bin/python
__author__ = 'xia'
print "Content-type: text/html\n\n"
from jinja2 import Template,Environment,PackageLoader


env=Environment(loader=PackageLoader(__name__))
template=env.get_template("hw.html")
print template.render(name="world")