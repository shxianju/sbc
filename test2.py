from jinja2 import Template,Environment,PackageLoader
print "Content-type: text/html\n\n"

env=Environment(loader=PackageLoader(__name__))
template=env.get_template("hw.html")
print template.render(name="world")