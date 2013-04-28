import re
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Circle
import matplotlib.lines as mlines
import numpy as np
from math import pi,sin,cos
from matplotlib.collections import PatchCollection
from Bio import SeqIO
from Bio.SeqFeature import *
import time
from os import remove
#sss=time.time()
__version__="v_1_2"
class plasmid():
    n=0
    annotate_size=6.#font size of annotation
    big_size=annotate_size*18./7.
    fig_size=annotate_size*9./7.
    coord_diff=annotate_size*0.02/7.
    wh_ratio=1.3#width/height ratio
    c_center=(0.5,0.5)
    c_radius=0.21
    block_width=0.04
    arrow_width=0.03
    text_Y=1.-coord_diff
    anno_xpad=.1
    anno_ypad=.1
    anno_ratio=1.35
    Arrow=dict(pro_ter=10.,gene=5.)
    etype_list=["gene", "insert", "pro_ter", "tag", "RE"]
    etype_colors=["red","lightgreen","lightblue","magenta","grey"]
    Deg_nt=dict(A="A",T="T",G="G",C="C",
                R="[AG]",Y="[CT]",M="[AC]",K="[GT]",
                S="[GC]",W="[AT]",H="[ATC]",B="[GTC]",
                V="[GAC]",D="[GAT]",N="[AGCT]")
    RE_dict=dict(Acc65I="G'GTACC",AflIII="A'CRYGT",ApaI="GGGCC'C",
                 AseI="AT'TAAT",AvaI="C'YCGRG",AvrII="C'CTAGG",
                 BamHI="G'GATCC",BbeI="GGCGC'C",BlpI="GC'TNAGC",
                 BmeT110I="CY'CGRG",BmtI="GCTAG'C",Bpu10I="CC'TNAGC",
                 BsoBI="C'YCGRG",BspDI="AT'CGAT",BspEI="T'CCGGA",
                 BspQI="",BsrBI="CCG'CTC",Bsu36I="CC'TNAGG",
                 BtgZI="",ClaI="AT'CGAT",Eco53kI="GAG'CTC",
                 EcoNI="CCTNN'NNNAGG",EcoRV="GAT'ATC",FspI="TGC'GCA",
                 HpaI="GTT'AAC",KasI="G'GCGCC",KpnI="GGTAC'C",
                 NarI="GG'CGCC",NheI="G'CTAGC",NotI="GC'GGCCGC",
                 PaeR7I="C'TCGAG",PasI="CC'CWGGG",PciI="A'CATGT",
                 PflFI="GACN'NNGTC",PspOMI="G'GGCCC",PstI="CTGCA'G",
                 PvuI="CGAT'CG",SacI="GAGCT'C",SapI="",
                 ScaI="AGT'ACT",SfoI="GGC'GCC",SphI="GCATG'C",
                 SspI="AAT'ATT",StuI="AGG'CCT",TliI="C'TGGAG",
                 Tth111I="GACN'NNGTC",XcmI="CCANNNNN'NNNNTGG",XhoI="C'TCGAG")
    #{"BamHI": "ggatcc", "EcoRI":"gaattc","KpnI":"ggtacc", "XhoI":"ctcgag"}
    #coord_dict={0:[1-anno_xpad,1-anno_ypad-coord_diff,'left',-1],
    #            1:[anno_xpad,1-anno_ypad-coord_diff,'right',-1],
    #            2:[anno_xpad,0+anno_ypad+coord_diff,'right',1],
    #            3:[1-anno_xpad,0+anno_ypad+coord_diff,'left',1]} #coordinate list of annotation in each sector    

    def __init__(self,file_name="./test.gb",construct_name="",out_path=""):
        self.file_name=file_name
        self.record = SeqIO.read(file_name, "genbank")
        self.whole_len=len(self.record)
        self.out_path=out_path
        for etype in self.etype_list:
            exec("self.%sList=[]"%etype)
        self.used_RE=[]
        if self.whole_len<100:
            print "Sequence too short! There is no plasmid like this!"
            exit()
        if construct_name=="":
            if re.search(r'/',file_name):
                self.construct_name=re.search(r'.+/(.+?)\.gb$',file_name).group(1)
            else:
                self.construct_name=re.search(r'^(.+)\.gb$',file_name).group(1)
			
        
    def append_block(self,start_theta,end_theta,etype,color="blue"):
        if etype=="RE":
            alpha=0.6
        else:
            alpha=1.
        block=Wedge(self.c_center,self.c_radius+self.block_width/2.,
        theta1=end_theta,theta2=start_theta,width=self.block_width,color=color,alpha=alpha)
        exec("self."+etype+"List.append(block)")
    def append_arrow(self,arrow_point,arrow_root,etype,color="blue"):
        arrow_dots=np.array([self.polar(arrow_point,radius=self.c_radius),
                             np.array(self.polar(arrow_root,radius=self.c_radius+self.arrow_width)),
                             np.array(self.polar(arrow_root,radius=self.c_radius-self.arrow_width))])
        arrow=Polygon(arrow_dots,color=color)
        exec("self."+etype+"List.append(arrow)")
#    def append_line(self,start_theta,end_theta,etype,color="blue"):
#        mid_theta=self.mid(start_theta, end_theta)
#        line = mlines.Line2D([self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[0],1.] ,
#                             [self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[1],0], lw=1., alpha=0.6,color=color)
#        exec("self."+etype+"List.append(line)")
    def append_annotation(self,mid_theta,text,etype):
        if mid_theta in self.anno_dict:
            if (text,etype) in self.anno_dict[mid_theta]:
                pass
            else:
                self.anno_dict[mid_theta].append((text,etype))
        else:
            self.anno_dict[mid_theta]=[(text,etype)]
    def draw_annotations(self):
        anno_sorted_list=[x for x in sorted(self.anno_dict.keys(),reverse=True) if int(x/90)==0]+\
                          [x for x in sorted(self.anno_dict.keys(),reverse=False) if int(x/90)==1]+\
                          [x for x in sorted(self.anno_dict.keys(),reverse=True) if int(x/90)==2]+\
                          [x for x in sorted(self.anno_dict.keys(),reverse=False) if int(x/90)==3]
        anno_sorted_list.reverse()
        for mid_theta in anno_sorted_list:
            sector=int(mid_theta/90)#sector of the angle,0:upper_right,1:upper_left,2:lower_left,3:lower_right
            indicator=self.coord_dict[sector][3]
            for anno in self.anno_dict[mid_theta]:
                (text,etype)=anno
                color=self.etype_colors[self.etype_list.index(etype)]
                preferred_x,preferred_y=self.polar(mid_theta,radius=self.c_radius*self.anno_ratio)
                if -1*indicator*(preferred_y-self.coord_dict[sector][1])>=self.coord_diff:#if there is enough space between the two annotations
                    this_x=preferred_x
                    this_y=preferred_y
                else:
                    this_x=preferred_x#self.coord_dict[sector][0]
                    this_y=self.coord_dict[sector][1] -self.coord_diff*indicator
                plt.text(this_x,this_y,text,fontsize=self.annotate_size,
                         horizontalalignment=self.coord_dict[sector][2],
                         verticalalignment='center',color=color)
                line = mlines.Line2D([self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[0],this_x] ,
                                     [self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[1],this_y],
                                     lw=.3, alpha=1.,color=color)
                self.ax.add_line(line)
                self.coord_dict[sector][0] =this_x
                self.coord_dict[sector][1] =this_y
                if this_y<0 or this_y>1:
                    self.n+=1
                    print "Warning: Too many features!",self.n
                    
                #plt.text(self.coord_dict[sector][0],self.coord_dict[sector][1],text,fontsize=self.annotate_size,
                #     horizontalalignment=self.coord_dict[sector][2],verticalalignment='center')
                #line = mlines.Line2D([self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[0],self.coord_dict[sector][0]] ,
                #                     [self.polar(mid_theta,radius=self.c_radius+.5*self.block_width)[1],self.coord_dict[sector][1]],
                #                     lw=.3, alpha=1.,color="black")
                #self.ax.add_line(line)
                #self.coord_dict[sector][1] += self.coord_diff*self.coord_dict[sector][3]
        
    def add_RE(self,start,end,name):
        etype="RE"
        self.add_insert(start,end,name,etype)
    def add_tag(self,start,end,name):
        etype="tag"
        self.add_insert(start,end,name,etype)
    def add_insert(self,start,end,name,etype="insert"):
        color=self.etype_colors[self.etype_list.index(etype)]
        start_theta,end_theta=self.bp2theta(start,end)
        self.append_block(start_theta,end_theta,etype,color=color)
        self.append_annotation(self.mid(start_theta,end_theta),name+"("+str(start)+".."+str(end)+")",etype)
    def add_pro_ter(self,start,end,name,complement=False,):#promoter/terminator
        etype="pro_ter"
        self.add_gene(start,end,name,complement,etype)
        #~ color=self.etype_colors[self.etype_list.index(etype)]
        #~ start_theta,end_theta=self.bp2theta(start,end)
        #~ arrow_theta=start_theta-end_theta
        #~ arrow_point=complement*start_theta+(1-complement)*end_theta
        #~ arrow_root=(1-complement)*start_theta+complement*end_theta
        #~ self.append_arrow(arrow_point,arrow_root,etype,color=color)
        #~ self.append_annotation(self.mid(start_theta,end_theta),name+"("+str(start)+".."+str(end)+"c"*complement+")",etype)
    def add_gene(self,start,end,name,complement=False,etype="gene"): #theta: angle,0~360
        color=self.etype_colors[self.etype_list.index(etype)]
        start_theta,end_theta=self.bp2theta(start,end)
        arrow_theta=self.Arrow[etype] #in theta
        if (start_theta-end_theta)%360>arrow_theta:
            arrow_point=(1-complement)*end_theta+complement*start_theta
            start_theta=(start_theta-complement*arrow_theta)%360
            end_theta=(end_theta+(1-complement)*arrow_theta)%360
            arrow_root=(1-complement)*end_theta+complement*start_theta
            self.append_arrow(arrow_point,arrow_root,etype,color=color)
            self.append_block(start_theta,end_theta,etype,color=color)
        else: #if gene too short
            arrow_theta=start_theta-end_theta
            arrow_point=complement*start_theta+(1-complement)*end_theta
            arrow_root=(1-complement)*start_theta+complement*end_theta
            self.append_arrow(arrow_point,arrow_root,etype,color=color)
        self.append_annotation(self.mid(start_theta*(1-complement)+arrow_point*complement,
                                        arrow_point*(1-complement)+end_theta*complement),
                               name+"("+str(start)+".."+str(end)+"c"*complement+")",etype)
    def bp2theta(self,start,end):
        return (90 - float(start-1)/self.whole_len * 360)%360, (90 - float(end)/self.whole_len * 360)%360
    def polar(self,theta,radius=0.25):
        return self.c_center[0]+radius*cos(theta/180.*pi),self.c_center[1]+radius*sin(theta/180.*pi)
    def mid(self,start_theta,end_theta):
        return (end_theta+(start_theta-end_theta)%360/2.)%360		
    def organize(self):#rules
        self.anno_dict= {}
        self.coord_dict={0:[self.c_center[0]+self.c_radius* self.anno_ratio,self.c_center[1]-.5*self.coord_diff,'left',-1],
                    1:[self.c_center[0]-self.c_radius* self.anno_ratio,self.c_center[1]-.5*self.coord_diff,'right',-1],
                    2:[self.c_center[0]-self.c_radius* self.anno_ratio,self.c_center[1]+.5*self.coord_diff,'right',1],
                    3:[self.c_center[0]+self.c_radius* self.anno_ratio,self.c_center[1]+.5*self.coord_diff,'left',1]}
                    #{mid_theta:[(name,etype),(name,etype)],mid_theta:[(name,etype)]}   |etype=element type|
        for feature in self.record.features:
            start=int(str(feature.location.start))
            end=int(str(feature.location.end))
            if feature.type=='gene': #gene
                self.add_gene(start,end,feature.qualifiers['gene'][0],complement=(1-feature.strand)/2)
            elif feature.type=='CDS' and feature.qualifiers.has_key("gene"): #gene
                self.add_gene(start,end,feature.qualifiers['gene'][0],complement=(1-feature.strand)/2)
            elif feature.type=='promoter': #pro_ter
                self.add_pro_ter(start,end,feature.qualifiers['note'][0],complement=(1-feature.strand)/2)
            elif feature.type=='terminator': #pro_ter
                self.add_pro_ter(start,end,feature.qualifiers['note'][0],complement=(1-feature.strand)/2)
            elif feature.type=='promoter': #pro_ter
                self.add_pro_ter(start,end,feature.qualifiers['note'][0],complement=(1-feature.strand)/2)
            elif feature.type=='CDS' and feature.qualifiers.has_key("product")\
                     and feature.qualifiers.has_key("note"): #tag
                if re.search(r"tag",feature.qualifiers["product"][0],re.I) \
                       or re.search(r"tag",feature.qualifiers["note"][0],re.I):
                    self.add_tag(start,end,feature.qualifiers['note'][0])
            elif feature.type=='insert':#insert
                self.add_insert(start,end,feature.qualifiers['note'][0])
            elif feature.type=='tag':#tag (customized tag)
                self.add_tag(start,end,feature.qualifiers['note'][0])
            elif feature.type=='CDS' and feature.qualifiers.has_key("translation"):
                if feature.qualifiers.has_key("note") and \
                (re.search(r"his$|^his",feature.qualifiers['note'][0],re.I) or
                re.search(r"myc$|^myc",feature.qualifiers['note'][0],re.I) or
                re.search(r"HA",feature.qualifiers['note'][0]) or
                re.search(r"tag",feature.qualifiers['note'][0],re.I)):#tag
                    self.add_tag(start,end,feature.qualifiers['note'][0])
                elif feature.qualifiers.has_key("gene"):#gene
                    self.add_gene(start,end,feature.qualifiers['gene'][0],complement=(1-feature.strand)/2)
                elif feature.qualifiers.has_key("note") and \
                not re.search(r"site$|atg",feature.qualifiers['note'][0],re.I):#gene
                    self.add_gene(start,end,feature.qualifiers['note'][0],complement=(1-feature.strand)/2)
                elif feature.qualifiers.has_key("product"):#gene
                    self.add_gene(start,end,feature.qualifiers['product'][0],complement=(1-feature.strand)/2)
            elif feature.type=='CDS' and not feature.qualifiers.has_key("translation") and feature.qualifiers.has_key("note"):
                print feature.qualifiers['note']
            for RE in ["ApaI","AvaI","BamHI","ClaI","EcoRV","HpaI","KpnI","NheI",
            "NotI","PasI","SacI","SphI","StuI","XhoI"]:#self.RE_dict:
                if self.RE_dict[RE]=="": continue
                RE_site=re.sub("'","",self.RE_dict[RE])
                if RE_site in self.used_RE: continue
                self.used_RE.append(RE_site)
                RE_pattern="".join([self.Deg_nt[x] for x in RE_site])
                for m in re.finditer(RE_pattern,str(self.record.seq)+str(self.record.seq)[:20],re.I):#+str(self.record.seq)[:10] add some to complete the RE search
                    start=(m.start()+1)%self.whole_len
                    end=m.end()%self.whole_len
                    self.add_RE(start,end,RE)
    def draw(self):
        self.fig = plt.figure(1, figsize=(self.fig_size*self.wh_ratio, self.fig_size))
        self.ax = self.fig.add_subplot(111)
        #ax=plt.axes()#[0,0,1,1])
        plt.axis("scaled")
        self.ax.xaxis.set_ticks_position('none')
        self.ax.set_xticklabels([])
        self.ax.yaxis.set_ticks_position('none')
        self.ax.set_yticklabels([])
        plt.axis('off')
        self.c=Circle(self.c_center,radius=self.c_radius,fill=False,edgecolor="grey")
        self.ax.add_patch(self.c)
        if len(self.construct_name)<=15:
            name_size=self.big_size
        elif len(self.construct_name)>15 and len(self.construct_name)<=22:
            name_size=self.big_size*.7
        else:
            name_size=self.big_size*.5
            print self.construct_name,len(self.construct_name)
        plt.text(self.c_center[0],self.c_center[1]*1.03,self.construct_name,fontsize=name_size,
                horizontalalignment='center',
                verticalalignment='bottom')
        plt.text(self.c_center[0],self.c_center[1]*.97,str(self.whole_len)+'bp',fontsize=self.big_size,
             horizontalalignment='center',
             verticalalignment='top')
        #self.ax.axhline(y=self.c_center[0])
        #self.ax.axvline(x=self.c_center[1])

        for etype in self.etype_list:
            for patch in eval("self."+etype+"List"):
                self.ax.add_patch(patch)
        self.draw_annotations()

        self.ax.set_xlim(0-(self.wh_ratio-1)/2.,1.+(self.wh_ratio-1)/2.)
        self.ax.set_ylim(0,1)
        fname=self.out_path+self.construct_name+re.sub(r"\.","",str(time.time()))+".png"
        plt.savefig(fname, bbox_inches="tight",dpi=300,pad_inches=0.)
        plt.clf()#plt.show()
        return fname

    def ins_seq(self,sequence,ins_name,ins_sites,etype): #the cutpoint is after ins_sites[0] bp and after ins_sites[1] bp
        from Bio.Alphabet import IUPAC
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        ins_record=SeqRecord(Seq(sequence,IUPAC.ambiguous_dna))
        f=SeqFeature(FeatureLocation(0,len(ins_record)),type=etype)
        f.qualifiers["note"]=[ins_name,]
        ins_record.features=[f]
        old_name=self.record.name
        self.record=self.record[:ins_sites[0]]+ins_record+self.record[ins_sites[1]:]
        self.record.name=old_name
        self.whole_len=len(self.record)
        #record_new.annotations["data_file_division"]=record_test.annotations["data_file_division"]
        #~ tmp_time=str(time.time())
        #~ SeqIO.write(self.record,"test_%s.gb"%tmp_time,"genbank")
        #~ self.__init__("test_%s.gb"%tmp_time,old_name)
        #~ remove("test_%s.gb"%tmp_time)
    def test(self):
        self.whole_len=4000
        self.add_RE(250,255,"RE1")
        self.add_gene(225,500,"gene1",complement=False)
        self.add_RE(595,600,"RE2")
        self.add_gene(2300,2800,"gene3")
        self.add_gene(799,1200,"gene2",complement=True)
        self.add_RE(999,1004,"RE3")
        self.add_tag(1860,1900,"tag3")
        self.add_tag(1610,1650,"tag2")
        self.add_RE(2850,2855,"RE3")
        self.add_tag(540,580,"tag1")
        self.add_tag(3300,3450,"tag4")
        self.add_RE(2850,2855,"RE6")
        self.add_pro_ter(3900,13,"promoter1",complement=True)
        self.add_RE(3800,3805,"RE5")
        self.add_RE(3500,3505,"RE4")
        self.add_insert(2800,3280,"insert1")
        self.add_RE(3998,4,"RE7")
        self.add_RE(3998,4,"RE7")
        self.draw()
def main():
    A=plasmid("test.gb")
    A.ins_seq("AAGCTTA"*500,"test tag",[1000,1025],"insert")
    A.organize()
    A.draw()
    exit()
    #print time.time()-sss

    #plasmid().test()
    #exit()
    #~ A=plasmid("pUC.gb")
    #~ A.organize()
    #~ A.draw()
    #~ exit()
    n=0
    from os import listdir
    for file_name in listdir("./tmp/"):
        if re.search(r"\.gb$",file_name):
            n+=1
            A=plasmid("./tmp/"+file_name)#pUC.gb")pET-28a(+) pbr322
            print A.construct_name,n
            A.out_path="./tmp_fig/"
            A.organize()
            A.draw()
if __name__=="__main__":
    main()
