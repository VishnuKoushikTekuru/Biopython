from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
import pylab
import os


Entrez.email="test@test.com"
print ("process initiated......")
v=input("enter a term you want to search for: ")
handle=Entrez.esearch(db="nucleotide", retmax=5, term= v, idtype="acc")
record=Entrez.read(handle)
print ("file read and ID lists retrieved.....\n")
handle.close()
ID_List=record['IdList']
print (ID_List)
for hj in ID_List:
    Entrez.email="test@test.com"
    handle=Entrez.efetch(db="nucleotide", id=hj, rettype="fasta",retmode="text")
    sr=SeqIO.read(handle,"fasta")
    handle.close()
    print ("sequnce length of ", hj,"is",len(sr.seq))

r=input("enter the desired ID from the above list you want to search for: ")
Entrez.email="test@test.com"
print ("step 2-fetching the data from the above ID's......\n") 
handle=Entrez.efetch(db="nucleotide", id=r, rettype="gb",retmode="text")
sr=SeqIO.read(handle,"genbank")
print ("we read the files and this is their length",len(sr.seq))
r=SeqIO.write(sr,sr.id+".fasta","fasta")
print("we have successfully downloaded the files in fasta format")
p=(sr.seq)
z=sr.id
# print (z)
print ("STEP-3....\n")

print ("Blast Process Started..")

if os.path.isfile(z+".fasta"):
    print("blasting using fasta")
    result_handle = NCBIWWW.qblast("blastn","nt",z+".fasta")
else:
    print ("blasting using sequence")
    result_handle = NCBIWWW.qblast("blastn","nt",p)
save_file=open("my_bla9.xml","w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
print("Done downloading the XML file from Blast...\n")
print ("STEP-4....\n")

print ("Now reopening the file to convert it to fasta")
rh = open("my_bla9.xml","r")
blast_record= NCBIXML.read(rh)
E_VALUE_THRESH=float(input("enter threshold value"))
f_handle=open("testing.fasta","w")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:        
        if hsp.expect< E_VALUE_THRESH:
            sr=SeqRecord(Seq(hsp.sbjct),id=alignment.hit_id)
            SeqIO.write(sr,f_handle,"fasta")
            # print("*****Alignment******")
            # print ("sequemce:", alignment.title)
            # print ("length:", alignment.length)
            # print ("e.value:", hsp.expect)
            # print (hsp.query[0:75]+"....")
            # print (hsp.match[0:75]+"...")
            # print (hsp.sbjct[0:75]+"...")
            break
print ("successfully downloaded in fasta format")

print ("STEP-5....\n")
cl_exe=r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
cl_cline= ClustalwCommandline(cl_exe, infile="testing.fasta")
stdout,stderr= cl_cline()
print ("file is converted")        
print ("STEP-6....\n")



align=AlignIO.read("opuntia.aln","clustal")
print(align)
print ("STEP-7....\n")

Phylo.convert("opuntia.dnd","newick","testing3.xml","phyloxml")
tree= Phylo.read("testing3.xml","phyloxml")
Phylo.draw_ascii(tree)
