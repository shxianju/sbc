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
    if len(true_seq)<=60:
        return "The input sequence should be longer than 60bp"
    elif len(true_seq)%3 != 0:
        return "Please input your CDS region, which means the length of your sequence should be multiple of 3."
    else:
        return "ok result"#pass


if not form:
    import mysql.connector
    config = {
        'user': 'root',
        'password': '123',
        'host': '127.0.0.1',
        'database': 'test',
        'raise_on_warnings': True,
        }
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()
    query = ("SELECT exp_sys,name,tag,cp_num,resist,path FROM test_vector ")
    cursor.execute(query)
    records=[]
    for (exp_sys,name,tag,cp_num,resist,path) in cursor:
        records.append({"exp_sys": exp_sys, "name": name, "tag": tag,
                        "cp_num": cp_num, "resist": resist, "path": path})
    cursor.close()
    cnx.close()
    template=env.get_template("index.html")
    print template.render(records=records)
elif form.has_key("goto_modification1"):
    template=env.get_template("modification.html")
    print template.render()
elif form.has_key("goto_report"):
    template=env.get_template("report.html")
    print template.render()
elif form.has_key("goto_muta"):
    template=env.get_template("mutagenesis.html")
    print template.render()
elif form.has_key("ins_seq_check"):
    print check_seq(form["ins_seq_check"].value)


