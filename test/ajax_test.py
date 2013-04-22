#!/usr/bin/python
__author__ = 'xia'
from jinja2 import Template,Environment,PackageLoader
import cgi,cgitb,re,random


cgitb.enable()
form=cgi.FieldStorage(keep_blank_values=True)
print "Content-type: text/html\n\n"
print form["a"].value,form["b"].value,random.random()