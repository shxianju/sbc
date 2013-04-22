#!/usr/bin/python
__author__ = 'xia'
from jinja2 import Template,Environment,PackageLoader
import cgi,cgitb,re

cgitb.enable()
form=cgi.FieldStorage(keep_blank_values=True)
env=Environment(loader=PackageLoader(__name__,"templates"))
print "Content-type: text/html\n\n"
def check_seq(ins_seq):
    if ins_seq.strip()=="":
        return "Please input your insert sequence."
    elif re.search(r"[^atgcu0-9\s]",ins_seq,re.I):
        return "Please only input a DNA/mRNA sequence"
    true_seq=re.sub(r"[\s0-9]+","",ins_seq)
    if len(true_seq)%3 != 0:
        return "Please input your CDS region, which means the length of your sequence should be multiples of 3."



if not form:
    template=env.get_template("index.html")
    print template.render()
elif form.has_key("goto_modification1"):
    template=env.get_template("modification.html")
    print template.render()
elif form.has_key("goto_report"):
    template=env.get_template("report.html")
    print template.render()
elif form.has_key("goto_muta"):
    template=env.get_template("mutagenesis.html")
    print template.render()
elif form.has_key("ins_seq"):
    print check_seq(form["ins_seq"].value)


