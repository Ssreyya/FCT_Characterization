# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:40:35 2022

@author: SSudharsan
"""



from Bio import SeqIO
import functools
import os

class FCT():
    def readgbk(self,filename):
        file="Annotated/" + filename + ".gbk"
        records = list(SeqIO.parse(file, "gb"))
        #accounts for genbank files with more than one record
        for i in range(len(records)):
            #list would include all gene features of that record
            genes=[]
            for feature in records[i].features:
                if feature.type=="gene":
                    genes.append(feature)
        
    
            gene_of_interest = "hslO"
            #look for our gene of interest
            for j in range(len(genes)):
                if 'gene' in genes[j].qualifiers.keys():
                    #print(genes[0].location.strand)
                    #print(genes[0].location.start)
                    #print(genes[0].location.end)
                    if gene_of_interest==genes[j].qualifiers['gene'][0]:
                        locus_hslO=genes[j].qualifiers['locus_tag'][0]
                        direction = genes[j].location.strand #strand 1 or -1
                        avglocPAS = self.getcoordinates(filename)
                        
                        if direction == 1:
                            direction_a=-1
                            if avglocPAS != None:
                                PASgenefeature=self.findlocusPAS1(genes,j,avglocPAS)
                            if PASgenefeature != None:
                                locus_PAS1=PASgenefeature.qualifiers['locus_tag'][0]
                                self.locusList1(genes,j,filename,records,locus_PAS1,i)
                                self.getFCTType(filename,locus_hslO,direction_a)
                        
                        if direction == -1:
                            direction_a=1
                            PASgenefeature = self.findlocusPAS2(genes,j,avglocPAS)
                            if PASgenefeature != None:
                                locus_PAS2=PASgenefeature.qualifiers['locus_tag'][0]
                                self.locusList2(genes,j,filename,records,locus_PAS2,i)
                                self.getFCTType(filename,locus_hslO,direction_a)
     
    def findlocusPAS1(self,genes,j,avglocPAS):
        for k in range(0,j+1):
            if k < 30:
                a=genes[j-k].location.start
                b=genes[j-k].location.end
                if avglocPAS in range(a,b):
                    return genes[j-k]
            
    def findlocusPAS2(self,genes,j,avglocPAS):
        l=len(genes) - j
        for k in range(0,l):
            if k <30:
                a=genes[j+k].location.start
                b=genes[j+k].location.end
                if avglocPAS in range(a,b):
                    return genes[j+k]

    
    def locusList1(self,genes,j,filename,records,locus_PAS1,i):#for postive strand of hslO
        locus_tag=genes[j].qualifiers['locus_tag'][0]
        extracted_sequence=genes[j].extract(records[i])
        extracted_sequence.id=locus_tag + "@" +filename
        extracted_sequence.description=locus_tag
        a="FCTgenes/"+filename+"FCTgenes.fasta"
        with open(a,'a') as f:
            SeqIO.write(extracted_sequence,f,'fasta')
        k=0
        while locus_tag != locus_PAS1:
            k=k+1
            locus_tag=genes[j-k].qualifiers['locus_tag'][0]
            extracted_sequence=genes[j-k].extract(records[i])
            extracted_sequence.id=locus_tag + "@" + filename
            extracted_sequence.description=locus_tag
            with open(a,'a') as f:
                SeqIO.write(extracted_sequence,f,'fasta')

    def locusList2(self,genes,j,filename,records,locus_PAS2,i): #for postive strand of hslO
        i=int(i)
        locus_tag=genes[j].qualifiers['locus_tag'][0]
        extracted_sequence=genes[j].extract(records[i])
        extracted_sequence.id=locus_tag + "@" + filename
        extracted_sequence.description=locus_tag
        a="FCTgenes/"+filename+"FCTgenes.fasta"
        with open(a,'a') as f:
            SeqIO.write(extracted_sequence,f,'fasta')
        k=0
        while locus_tag != locus_PAS2:
            k=k+1
            locus_tag=genes[j+k].qualifiers['locus_tag'][0]
            extracted_sequence=genes[j+k].extract(records[i])
            extracted_sequence.id=locus_tag + "@" + filename
            extracted_sequence.description=locus_tag
            with open(a,'a') as f:
                SeqIO.write(extracted_sequence,f,'fasta')
        
        
                
    def getFCTType(self,filename,locus_hslO,direction_a):
        hslO_locus = int(locus_hslO.split("_")[1])
        q="FCTgenes/"+filename+"FCTgenes.fasta"
        makedb = "makeblastdb -in " + "Ref_Ann_genes.fasta" + " -dbtype nucl -out " + "refblastdb"
        os.system(makedb)
        blast = " blastn -query " + q + " -db refblastdb "+ " -evalue 0.05 -outfmt '6 delim=@ sacc qacc' -out temp/testblast"+filename+".txt"
        os.system(blast)
        
        c="temp/testblast"+filename+".txt"

    
        #opening blastoutput file and reading each line
        with open(c) as f:
            lines=f.read().splitlines()
        genelist=[]
        locusngenes={}
        
        with open(q) as f:
            lines1=f.read().splitlines()
            linelist=[]
            for line in lines1:
                if ">" in line:
                    linelist.append(int(line.split(" ")[1].split("_")[1]))

        #putting all values in dictionary, dictionary also ignores repeated keys(here, the keys will be the locustag)
        for line in lines:
            locusngenes.setdefault(int(line.split("@")[1].split("_")[1]), line.split("@")[0])
        curr_locus=int()
        for key in locusngenes:
            if key == hslO_locus:
                genelist.append(locusngenes.get(key))
                curr_locus = key+(1*direction_a)
            else:
                if key == curr_locus:
        
                    genelist.append(locusngenes.get(key))
                    
                else:
                    a=(key-curr_locus)*(direction_a)
                    for i in range(a):
                        genelist.append("nohits"+str(curr_locus+(i*direction_a)))
                    genelist.append(locusngenes.get(key))
                 
                curr_locus=key+(1*direction_a)

        print(locusngenes)
                

        #to remove locus_numbers that originally werent even present (otherwise they would be reported as nohits)
        for gene in genelist:
            if "nohits" in gene:
                if int(gene[6:]) not in linelist:
                    genelist.remove(gene)

        genelistnew = genelist.copy() #making a copy of genelist in genelistnew to store transposases
        
        for gene in genelistnew:
            if "nohits" in gene:
                cr="temp/"+filename+gene+"hits.fasta"
                with open(q) as f, open(cr,'w') as z:
                    lines1=f.read().splitlines()
                    
                    i=0
                    while  i>=0 and i<len(lines1):
                        if gene[6:] in lines1[i]:
                            z.write(lines1[i])
                            z.write("\n")
                            start=i
                        i=i+1
                    i=start+1
                    while i>start and i<len(lines1) and ">" not in lines1[i]:
                        z.write(lines1[i])
                        i=i+1
                blastxl = "blastx -db nrdb/nr -query " + cr +  " -taxids '1301'   -max_target_seqs 10 -outfmt '6 delim=@  sacc qacc stitle' -out nohit.txt"
                os.system(blastxl)
                with open("nohit.txt") as g:
                    lines2=g.read().splitlines()
                    if len(lines2)>0:
                        foundtranspo="FALSE"
                        for line in lines2:
                            if "transposase" in line.lower():
                                genelistnew = list(map(lambda x: x.replace(gene, 'transposase'), genelistnew))
                                foundtranspo="TRUE"
                                break
                        if foundtranspo=="FALSE":
                            genelistnew = list(map(lambda x: x.replace(gene, lines2[0].split("@)")[0]), genelistnew))
                    


                        
        
        with open("genelist.txt",'a') as f:
            f.write(filename)
            for gene in genelist:
                f.write(","+gene)
            f.write("\n")

        with open("genelistnew.txt",'a') as f:
            f.write(filename)
            for gene in genelistnew:
                f.write(","+gene)
            f.write("\n")
                    
        #removing nohits and sortases from our genelist and creating hitgeneslist
        hitgenes=[]
        for gene in genelist:
            if "transposase" not in gene and "sortase" not in gene:
                hitgenes.append(gene)
        
        hitgenes1=[]
        for gene in genelistnew:
            if "transposase" not in gene and "sortase" not in gene:
                hitgenes1.append(gene)
        
        with open("coregenes.txt",'a') as f:
            f.write(filename)
            for gene in hitgenes1:
                f.write(","+gene)
            f.write("\n")
        

        fct_types={"FCT_1":['hslO', 'rofA', 'bp', 'ap2', 'ap1', 'fbp', 'yheo'],
        "FCT_2":['hslO','rofA','PrtF1/sfbl', 'bp', 'ap2', 'ap1', 'fbp', 'yheo'],
        "FCT_3":['hslO', 'rofA', 'PrtF1/sfbl', 'bp', 'fszE/ap2', 'fbp', 'yheo'],
        "FCT_4":['hslO', 'rofA', 'bp', 'fszE/ap2', 'fbp', 'yheo'],
        "FCT_5":['hslO', 'YSIRK_TR', 'bp', 'fszE/ap2', 'fbp', 'yheo'],
        "FCT_6":['hslO', 'rofA', 'PrtF1/sfbl', 'yheo'],
        "FCT_7":['hslO', 'rofA', 'PrtF1/sfbl', 'collagen_adhesion_protein', 'fimbrial_protein', 'yheo'],
        "FCT_8":['hslO', 'rofA', 'ap1', 'signal_peptidase', 'ap', 'yheo'],
        "FCT_9":['hslO', 'rofA', 'ap1','fbp', 'bp', 'bp', 'ap2', 'yheo'],
        "FCT_10":['hslO', 'YSIRK_TR', 'ap1', 'signal_peptidase', 'bp', 'pilin_related_fszE', 'fbp', 'yheo'],
        "FCT_11":['hslO', 'rofA', 'PrtF1/sfbl', 'fbp', 'bp', 'ap', 'yheo'],
        "FCT_12":['hslO','rofA','PrtF1/sfbl','fszE/ap2', 'fbp','yheO']}
        
    
        
        for key in fct_types:
            if functools.reduce(lambda i, j : i and j, map(lambda m, k: m == k, fct_types[key], hitgenes), True): 
                if len(hitgenes) == len(fct_types[key]):
                    with open("FCTtyping.txt",'a') as f:
                        f.write(filename+","+key)
                        for gene in genelist:
                            f.write(","+gene)
                        f.write("\n")
            
        
        
                
     
                            
    def getcoordinates(self,filename):
        file="Assemblies/" + filename+".fasta"
        blast = "blastn  -query PAS_SDSE.fasta -subject " + file  + " -outfmt '6 delim=@ sstart send qstart qend sseqid' " + " -evalue 0.005 -out temp/" + filename + "PASout.txt"
        os.system(blast)
        c= "temp/" + filename + "PASout.txt"
        with open(c) as f:
            lines=f.read().splitlines()
            if len(lines)>0:
                coordinate1 = int(lines[0].split("@")[0])
                coordinate2 = int(lines[0].split("@")[1])
            
                coordinate= int((coordinate1+coordinate2)/2)
        
                return coordinate
    
    def readfilenames(self):
        with open("ID_SDSE.txt") as f:
            filenames=f.read().splitlines()
            for filename in filenames:
               self.readgbk(filename)
            
ob=FCT()
ob.readfilenames()
