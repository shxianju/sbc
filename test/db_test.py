__author__ = 'xia'
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
for (exp_sys,name,tag,cp_num,resist,path) in cursor:
    print "exp_sys: %s, name: %s, tag: %s, cp_num: %s, resist: %s, path: %s"\
          %(exp_sys,name,tag,cp_num,resist,path)
cursor.close()
cnx.close()