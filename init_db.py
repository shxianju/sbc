import mysql.connector,os,re
from mysql.connector import errors
config = {
    'user': 'root',
    'password': '123',
    'host': '127.0.0.1',
    'database': 'test',
    'raise_on_warnings': True,
    }
DB_NAME = 'test'
TABLES = {}
TABLES['test_vector'] = (
    "CREATE TABLE `test_vector` ("
    "  `seq_no` int(11) NOT NULL AUTO_INCREMENT,"
    "  `exp_sys` varchar(20) DEFAULT NULL,"
    "  `name` varchar(100) NOT NULL,"
    "  `tag` varchar(10) NOT NULL,"
    "  `cp_num` varchar(4) NOT NULL,"
    "  `resist` varchar(10) DEFAULT NULL,"
    "  `path` varchar(200) NOT NULL,"
    "  PRIMARY KEY (`seq_no`)"
    ") ENGINE=InnoDB")

cnx = mysql.connector.connect(**config)
cursor = cnx.cursor()
for name, ddl in TABLES.items():
    try:
        print("Creating table {}: ".format(name))
        cursor.execute(ddl)
    except mysql.connector.Error as err:
        print "Error",
        print(str(err))
    else:
        print("OK")
add_seq = ("INSERT INTO test_vector "
                "(exp_sys, name, tag, cp_num, resist, path)"
                "VALUES (%s, %s, %s, %s, %s,%s)")
files_path="./PlasmidData/tmp/"
for path in os.listdir(files_path):
    name= path.strip(".gb")
    path= files_path+path
    tag_match=re.search(r"(myc|his)",name,re.I)
    if tag_match:
        tag=tag_match.group(1)
        print tag
    else:tag="no_tag"
    resist_match=re.search(r"(amp|neo|kan)",name,re.I)
    if resist_match:
        resist=resist_match.group(1)
        print resist
    else: resist="no_resist"
    data=("e_coli",name,tag, "high",resist,path)
    cursor = cnx.cursor()
    cursor.execute(add_seq, data)
    print("Inserting "+name)
cnx.commit()
cursor.close()
cnx.close()
