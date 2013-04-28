#coding=utf-8
import os
import re
path="./commercial_vectors/genbank/"
path_w="./tmp/"
file_names=os.listdir(path)
date_dict={u"一月": "Jan",u"二月":"Feb",u"三月":"Mar",u"四月":"Apr",
			u"五月":"May",u"六月":"Jun",u"七月":"Jul",u"八月":"Aug",
			u"九月":"Sep",u"十月":"Oct",u"十一月":"Nov",u"十二月":"Dec"}
def change_month(match):
	return date_dict[match.group()]
for file_name in file_names: 
	file_content=unicode(open(path+file_name).read(),"gbk")
	file_content= re.sub("|".join(date_dict.keys()),change_month,file_content)
	file_content=re.sub(u"月|年","-",file_content)
	file_content=re.sub(u"日","",file_content)
	line_1 = file_content.split("\n")[0]
	line_1=re.sub("Exported File","Exported_File",line_1)
	if line_1[55:63]=="ircular"or "inear":
		line_1=line_1[:54] + " " + line_1[54:]
	file_content="\n".join([line_1]+file_content.split("\n")[1:])
	open(path_w+file_name[:-1],'w').write(file_content)
