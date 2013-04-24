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
records=[]
for (exp_sys,name,tag,cp_num,resist,path) in cursor:
    records.append({"exp_sys": exp_sys, "name": name, "tag": tag,
           "cp_num": cp_num, "resist": resist, "path": path})
cursor.close()
cnx.close()