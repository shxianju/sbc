#!/usr/bin/python
__author__ = 'xia'
from jinja2 import Template,Environment,PackageLoader
import cgi,cgitb,re,os,time,json
DEBUG=1
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
table_codon={}
for i in range(len(amino_acids)):
    aa =amino_acids[i]
    if aa in table_codon.keys():
        table_codon[aa].append(codons[i])
    else:
        table_codon[aa]=[codons[i],]
hydrophobicity_table=dict(
        I="#ff0000",V="#f60009",
        L="#ea0015",F="#cb0034",
        C="#c2003d",M="#b0004f",
        A="#ad0052",G="#6a0095",
        X="#680097",T="#61009e",
        S="#5e00a1",W="#5b00a4",
        Y="#4f00b0",P="#4600b9",
        H="#1500ea",E="#0c00f3",
        Z="#0c00f3",Q="#0c00f3",
        D="#0c00f3",B="#0c00f3",
        N="#0c00f3",K="#0000ff",
        R="#0000ff")
hydrophobicity_table["*"]="#000000"
cgitb.enable()
form=cgi.FieldStorage(keep_blank_values=True)
env=Environment(loader=PackageLoader(__name__,"templates"))
print "Content-type: text/html\n\n"
def translate_format(nt_seq):
    if len(nt_seq)%3 != 0:
        raise ValueError("length of the input sequence %s is not multiples of 3."%nt_seq)
    format_dna=""
    format_protein=""
    for i in range(len(nt_seq)/3):
        codon = nt_seq[3*i:3*i+3].lower()
        aa=codon_table[codon].upper()
        aa_color=hydrophobicity_table[aa]
        format_dna += codon+"&nbsp;"
        format_protein +="<span style='color:%s;'>&nbsp;%s&nbsp;&nbsp;</span>"%(aa_color,aa)
        if (i+1)%20==0:
            format_dna += "<br>"
            format_protein += "<br>"
    return format_protein,format_dna



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
def vector_in_db():
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
        records.append({"exp_sys": exp_sys, "name": name, "tag": tag.lower(),
                        "cp_num": cp_num, "resist": resist.lower(), "path": path})
    cursor.close()
    cnx.close()
    return records
def get_sticky_end(re_pattern):#re_pattern in the form of "a|tcgat"
    if re_pattern:
        cut_pos=re_pattern.index("|")
        sticky_end=re.sub(r"\|","",re_pattern)[min(cut_pos,len(re_pattern)-1-cut_pos):
                                               max(cut_pos,len(re_pattern)-1-cut_pos)]
        if cut_pos<(len(re_pattern)-1)/2.: return "+"+sticky_end
        elif cut_pos>(len(re_pattern)-1)/2.:return "-"+sticky_end
        else: return ""
    else: return ""
def insert_insertion(form):
    vec_sites=form.getlist("ez_pos")
    vec_res=[re.search(r"\((.*)\)",x).group(1) for x in form.getlist("ez_txt")]
    vec_re_names=[re.search(r"(.+)\(.*\)",x).group(1) for x in form.getlist("ez_txt")]
    sign=cmp(int(vec_sites[0])+vec_res[0].index("|"),int(vec_sites[1])+vec_res[1].index("|"))
    vec_5_site=vec_sites[int(.5+sign*.5)]
    vec_5_re=vec_res[int(.5+sign*.5)]
    vec_3_site=vec_sites[1-int(.5+sign*.5)]
    vec_3_re=vec_res[1-int(.5+sign*.5)]
    ins_dict=eval(form["ins_dict"].value)
    ins_dict["vec_re_names"]=vec_re_names
    utr_5=ins_dict["utr_5"]
    utr_5_re=ins_dict["utr_5_re"]
    ins_seq=ins_dict["ins_seq"]
    utr_3_re=ins_dict["utr_3_re"]
    utr_3=ins_dict["utr_3"]
    vec_5_site,utr_5_seq,ins_seq,utr_3_seq,vec_3_site=ins_organize(vec_5_site, vec_5_re, utr_5, utr_5_re,
        ins_seq, utr_3_re, utr_3, vec_3_re, vec_3_site)
    print "<br>here goes test ",form.getvalue("path_dict")
    path_dict=eval(form.getvalue("path_dict"))
    vector_path_gb=path_dict["vector_path_gb"]
    out_path=path_dict["out_path"]
    old_vector=dp.plasmid(vector_path_gb,out_path=out_path,step_name="after_modification")
    old_vector.ins_insert(vec_5_site,utr_5_seq,ins_seq,utr_3_seq,vec_3_site,"ins_name")
    modified_vector_path=old_vector.write_gb_file()#add insert to old vector and update it
    return modified_vector_path,out_path,ins_dict
def insert_tag(form):
    mcs_dict=eval(form.getvalue("mcs_dict"))
    ins_dict=eval(form.getvalue("ins_dict"))
    ins_start=mcs_dict["ins_list"][0][0]
    ins_end=mcs_dict["ins_list"][0][1]
    tag_5_val=form.getvalue("tag_5_val")
    tag_5_protease=form.getvalue("tag_5_protease")
    tag_3_val=form.getvalue("tag_3_val")
    tag_3_protease=form.getvalue("tag_3_protease")
    path_dict=eval(form.getvalue("path_dict"))
    vector_path_gb=path_dict["vector_path_gb"]
    out_path=path_dict["out_path"]
    old_vector=dp.plasmid(vector_path_gb,out_path=out_path,step_name="getting_advanced_modification")
    if tag_5_val!="none" and tag_5_protease!="none":
        old_vector.ins_tag(tag_5_val.split(":")[1],tag_5_protease,tag_5_val.split(":")[0],
        [ins_start,ins_start],side=5)
        ins_end += len(tag_5_val.split(":")[1]+tag_5_protease)
        ins_dict["tag_5_seq"]=tag_5_val.split(":")[1]
        ins_dict["tag_5_name"]=tag_5_val.split(":")[0]
        ins_dict["tag_5_protease"]=tag_5_protease

    if tag_3_val!="none" and tag_3_protease!="none":
        old_vector.ins_tag(tag_3_val.split(":")[1], tag_3_protease,tag_3_val.split(":")[0],
        [ins_end,ins_end],side=3)
        ins_dict["tag_3_seq"]=tag_3_val.split(":")[1]
        ins_dict["tag_3_name"]=tag_3_val.split(":")[0]
        ins_dict["tag_3_protease"]=tag_3_protease
    modified_vector_path=old_vector.write_gb_file()
    return modified_vector_path,out_path,ins_dict

def reverse_complement_dna(seq):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    return str(Seq(seq,IUPAC.ambiguous_dna).reverse_complement())

def ins_organize(vec_5_site,vec_5_re,utr_5,utr_5_re,ins_seq,
                 utr_3_re,utr_3,vec_3_re,vec_3_site):
    vec_5_site=int(vec_5_site)+vec_5_re.index("|")
    vec_5_sticky=get_sticky_end(vec_5_re)
    utr_5_seq=re.search(
        re.sub(r"\|","(",dp.deg_pattern(utr_5_re,change_slash=False)+r".*)"),
        utr_5).group(1)
    utr_5_sticky=get_sticky_end(utr_5_re)

    utr_3_sticky=get_sticky_end(utr_3_re)
    utr_3_seq=re.search(
        re.sub(r"\|",")",r"(.*"+dp.deg_pattern(utr_3_re,change_slash=False)),
        utr_3).group(1)
    vec_3_sticky=get_sticky_end(vec_3_re)
    vec_3_site=int(vec_3_site)+vec_3_re.index("|")
    return vec_5_site,utr_5_seq,ins_seq,utr_3_seq,vec_3_site

def report_check(form):
    path_dict=eval(form.getvalue("path_dict"))
    ins_dict=eval(form.getvalue("ins_dict"))
    mcs_dict=eval(form.getvalue("mcs_dict"))
    entries=[]
    icon,item,content,comment=report_check_start_stop_codon(ins_dict["ins_seq"])
    entries.append((icon,item,content,comment))
    if ins_dict.has_key("tag_5_name"):
        icon,item,content,comment=report_check_5_tag(ins_dict)
        entries.append((icon,item,content,comment))
    if ins_dict.has_key("tag_3_name"):
        icon,item,content,comment=report_check_3_tag(ins_dict)
        entries.append((icon,item,content,comment))
    icon,item,content,comment=report_check_res(ins_dict,path_dict)
    entries.append((icon,item,content,comment))
    return entries

def report_check_res(ins_dict,path_dict):
    import draw_plasmid as dp
    vector_path_gb=path_dict["vector_path_gb"]
    vector_rep_seq=str(dp.plasmid(vector_path_gb).record.seq)
    re_names=ins_dict["vec_re_names"]+\
             [ins_dict["utr_5_re_name"]]+[ins_dict["utr_3_re_name"]]
    RE_dict=dict(Acc65I="G'GTACC",AflIII="A'CRYGT",ApaI="GGGCC'C",
        AseI="AT'TAAT",AvaI="C'YCGRG",AvrII="C'CTAGG",
        BamHI="G'GATCC",BbeI="GGCGC'C",BlpI="GC'TNAGC",
        BmeT110I="CY'CGRG",BmtI="GCTAG'C",Bpu10I="CC'TNAGC",
        BsoBI="C'YCGRG",BspDI="AT'CGAT",BspEI="T'CCGGA",
        BspQI="",BsrBI="CCG'CTC",Bsu36I="CC'TNAGG",
        BtgZI="",ClaI="AT'CGAT",Eco53kI="GAG'CTC",
        EcoNI="CCTNN'NNNAGG",EcoRI="G'AATTC",
        EcoRV="GAT'ATC",FspI="TGC'GCA",
        HpaI="GTT'AAC",KasI="G'GCGCC",KpnI="GGTAC'C",
        NarI="GG'CGCC",NheI="G'CTAGC",NotI="GC'GGCCGC",
        PaeR7I="C'TCGAG",PasI="CC'CWGGG",PciI="A'CATGT",
        PflFI="GACN'NNGTC",PspOMI="G'GGCCC",PstI="CTGCA'G",
        PvuI="CGAT'CG",SacI="GAGCT'C",SapI="",
        ScaI="AGT'ACT",SfoI="GGC'GCC",SphI="GCATG'C",
        SspI="AAT'ATT",StuI="AGG'CCT",TliI="C'TGGAG",
        Tth111I="GACN'NNGTC",XcmI="CCANNNNN'NNNNTGG",XhoI="C'TCGAG")
    repeat_list=[]#record the repeat enzymes
    for re_name in set(re_names):
        re_seq=re.sub(r"'","",RE_dict[re_name])
        re_pattern=dp.deg_pattern(re_seq)
        if re.search(re_pattern,vector_rep_seq,re.I):
            repeat_list.append(re_name)
    if repeat_list:
        return "warning","The restriction enzymes used to build the construct.",\
                "The REs you have used are %s."%", ".join(set(re_names)),\
                "%s still exist in your construct, please be careful."%", ".join(repeat_list)
    else:
        return "check","The restriction enzymes used to build the construct.",\
               "The REs you have used are %s."%", ".join(set(re_names)),\
                "None of them will cut your construct, it is safe to use."

def report_check_5_tag(ins_dict):
    if ins_dict.has_key("tag_5_name"):
        return "check","5' tag added by customer.",\
                "You have added %s to 5' end near the insertion"%ins_dict["tag_5_name"],\
                "none"
def report_check_3_tag(ins_dict):
    if ins_dict.has_key("tag_3_name"):
        return "check","3' tag added by customer.",\
               "You have added %s to 3' end near the insertion"%ins_dict["tag_3_name"],\
               "none"
def report_check_start_stop_codon(ins_seq):
    start_codon=ins_seq[:3]
    end_codon=ins_seq[-3:]
    if start_codon.upper()=="ATG" and end_codon.upper() in ("TAG","TGA","TAA"):
        return "check","Start and stop codon",\
               "Your start and stop codons are %s and %s"\
               %(start_codon.upper(),end_codon.upper()),\
                "none"
    elif start_codon.upper()!="ATG" and not end_codon.upper() in ("TAG","TGA","TAA"):
        return "error","Start and stop codon",\
               "Your start and stop codons are %s and %s"\
               %(start_codon.upper(),end_codon.upper()),\
               "None of your start/end codon is a valid codon for protein expression."
    else:
        return "warning","Start and stop codon",\
               "Your start and stop codons are %s and %s"\
               %(start_codon.upper(),end_codon.upper()),\
                "One of your start/end codons may not be valid codon for protein expression."


if not form:
    records=vector_in_db()
    template=env.get_template("index.html")
    print template.render(records=records)
elif form.has_key("goto_modification"):
    import draw_plasmid as dp
    if DEBUG==1:print form
    for i in form:
        print i,form[i].value,"<br>"
    user_name="test"
    if form.has_key("advance"): button_name="goto_advance"
    else:button_name="goto_review" #decide whether to go to advanced mode
    vector_path=form["vector"].value
    out_path="tmp/"+user_name+"___"+re.sub(r"\.","",str(time.time()))+"/"
    os.mkdir(out_path)
    vector_um=dp.plasmid(vector_path,out_path=out_path,step_name="goto_modification")#vector_um:unmodified vector
#    open("./"+"test.txt","w").write(str(vector_um.get_mcs()))
    mcs_dict=vector_um.get_mcs()
    mcs_path=vector_um.draw_mcs()
    vector_um.organize()
    mcs_ac=[]
    for start_pos,end_pos,RE_name,RE_seq in mcs_dict["RE_list"]:
        RE_seq=re.sub("'","|",vector_um.RE_dict[RE_name])
        mcs_ac.append(dict(value=str(start_pos),
                    label=str(start_pos)+": "+RE_name+"("+RE_seq+")",
                    desc=RE_name+"("+RE_seq+")"))
    mcs_dict["mcs_ac"]=sorted(mcs_ac,key=lambda x:int(x["value"]))
    vector_um_path= vector_um.draw()#this function returns the path of figure
    utr_5_re=form["5re"].value.split(":")[1]
    utr_3_re=form["3re"].value.split(":")[1]
    utr_5=form["5utr"].value
    utr_3=form["3utr"].value
    ins_seq=form["ins_seq"].value
    ins_dict=dict(utr_5_re=utr_5_re,utr_5_re_name=form["5re"].value.split(":")[0],
        utr_3_re=utr_3_re,utr_3_re_name=form["3re"].value.split(":")[0],
        utr_5=utr_5,utr_3=utr_3,ins_seq=ins_seq)

    path_dict=dict(mcs_path=mcs_path,vector_path=vector_um_path,
        vector_path_gb=form["vector"].value,out_path=out_path,
        original_vector_path_gb=form["vector"].value)
    template=env.get_template("modification.html")
    print template.render(path_dict=path_dict, button_name=button_name,
        mcs_dict=mcs_dict,ins_dict=ins_dict)

elif form.has_key("goto_advance"):
    import draw_plasmid as dp
    if DEBUG==1:print form
    modified_vector_path,out_path,ins_dict=insert_insertion(form)
    vector_path_gb=modified_vector_path
    vector_ins=dp.plasmid(vector_path_gb,out_path=out_path)
    mcs_dict=vector_ins.get_mcs()
    path_dict=eval(form.getvalue("path_dict"))
    mcs_path=vector_ins.draw_mcs()
    vector_ins.organize()
    vector_ins_path= vector_ins.draw()#this function returns the path of figure
    path_dict["mcs_path"]=mcs_path
    path_dict["vector_path"]=vector_ins_path
    path_dict["vector_path_gb"]=vector_path_gb
    path_dict["out_path"]=out_path
    template=env.get_template("advance_step.html")
    print template.render(ins_dict=ins_dict,path_dict=path_dict,
        button_name="goto_review",mcs_dict=mcs_dict)

elif form.has_key("goto_review"):
    import draw_plasmid as dp
    if DEBUG==1:print form
    if form.has_key("tag_5_val"):#this form come from tag page,deal with add tag method
        modified_vector_path,out_path,ins_dict=insert_tag(form)
    else:#this form come from modification page,deal with cut and paste method
        modified_vector_path,out_path,ins_dict=insert_insertion(form)
    vector_path_gb=modified_vector_path
    vector_rev=dp.plasmid(vector_path_gb,out_path=out_path)
    mcs_dict=vector_rev.get_mcs()
    mcs_path=vector_rev.draw_mcs()
    vector_rev.organize()
    vector_rev_path= vector_rev.draw()#this function returns the path of figure
    path_dict=dict(mcs_path=mcs_path,vector_path=vector_rev_path,
        vector_path_gb=vector_path_gb,out_path=out_path,
        original_vector_path_gb=eval(form.getvalue("path_dict"))["original_vector_path_gb"])
    template=env.get_template("review_step.html")
    print template.render(path_dict=path_dict, button_name="goto_report",
        mcs_dict=mcs_dict,ins_dict=ins_dict)


elif form.has_key("goto_report"):
    if DEBUG==1:print form
    path_dict=eval(form.getvalue("path_dict"))
    ins_dict=eval(form.getvalue("ins_dict"))
    check_result=report_check(form)
    vector_path_gb=path_dict["vector_path_gb"]
    format_protein,format_dna=translate_format(ins_dict["ins_seq"])
    format_seqs=format_protein,format_dna
    for entry in check_result:
        if entry[0]=="error":
            template=env.get_template("report_error.html")
            print template.render(path_dict=path_dict,
                mcs_dict=eval(form.getvalue("mcs_dict")),
                ins_dict=ins_dict,check_result=check_result,
                format_seqs=format_seqs)
            exit()
    template=env.get_template("report.html")
    print template.render(path_dict=path_dict,
        mcs_dict=eval(form.getvalue("mcs_dict")),
        ins_dict=ins_dict,check_result=check_result,format_seqs=format_seqs)


elif form.has_key("goto_muta"):
    if DEBUG==1:print form
    path_dict=eval(form.getvalue("path_dict"))
    vector_path_gb=path_dict["vector_path_gb"]
    ins_dict=eval(form.getvalue("ins_dict"))
    format_protein,format_dna=translate_format(ins_dict["ins_seq"])
    format_seqs=format_protein,format_dna
    template=env.get_template("mutagenesis.html")
    print template.render(format_seqs=format_seqs,protein_seq="".join(re.findall(r"[A-Z*]",format_protein)))
elif form.has_key("in_seq_mut"):
    in_seq=form["in_seq_mut"].value
    trans_list=[str(i+1)+": "+in_seq[i] for i in range(len(in_seq))]
    target_list=[]
    for aa in sorted(table_codon.keys()):
        for codon in sorted(table_codon[aa]):
            target_list.append("%s: %s"%(aa,codon))
    print json.dumps({"source_list":trans_list,"target_list":target_list})
elif form.has_key("ins_seq_check"):
    print check_seq(form["ins_seq_check"].value)

