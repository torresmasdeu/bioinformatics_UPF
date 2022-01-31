from Bio import SeqIO
from Bio.Seq import Seq
import os 
import pandas

exonerate_output = 1 

ParentDir = "/run/user/1400982974/gvfs/smb-share:server=fs-aules.ds.upf.edu,share=public/20428/Treball/HB"

HBGenome = "/mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.fa" 

def run_cmd(cmd):
    out_stat = os.system(cmd)
    if out_stat!=0:
        raise ValueError("The command '%s' did not work"%cmd)

PredictedDir = "%s/Predicted_proteins"%(ParentDir)
if not os.path.isdir(PredictedDir):
    run_cmd("mkdir %s"%PredictedDir)

df_total=pandas.DataFrame()

for record_query in SeqIO.parse(ParentDir+"/query.txt", "fasta"): #
    record_query.id = record_query.description.split("(")[1].replace(") # ","_").strip()
    protein_name = record_query.id.split("_")[0].strip()

    ProtDir = "%s/%s"%(ParentDir, record_query.id)
    if not os.path.isdir(ProtDir):
        run_cmd("mkdir %s"%ProtDir)
      
    query_protein = "%s/%s.fasta"%(ProtDir, record_query.id)
    record_query.description = ""
    record_query.name = ""
    
    record_query.seq = Seq(str(record_query.seq).replace("U","X")) 

    SeqIO.write([record_query], query_protein, "fasta")

    blast_output = "%s/%s.blast"%(ProtDir, record_query.id)
    
    if not os.path.isfile(blast_output): 
        run_cmd('tblastn -evalue 0.001 -query %s -db %s -out %s -outfmt "6 qstart qend sseqid sstart send sstrand evalue"'%(query_protein,HBGenome,blast_output))

    df_blast = pandas.read_csv(blast_output, sep="\t", names=["qstart", "qend", "sseqid", "sstart", "send", "sstrand", "evalue"], header=None)
    df_blast["sstrand"] = (df_blast.sstart < df_blast.send).map({True:"+", False:"-"}) 
    for scaffold in set(df_blast.sseqid): 
        for strand in set(df_blast.sstrand):
            df_scaffst = df_blast[(df_blast['sseqid'] == scaffold)&(df_blast['sstrand'] == strand)]

            if len(df_scaffst)==0: 
                continue
            if strand =="+":
                df_scaffst = df_scaffst.sort_values(by="sstart", ascending=True) 
            else:
                df_scaffst = df_scaffst.sort_values(by="sstart", ascending=False) 

            previous_qstart = 0
            previous_send = df_scaffst.send.iloc[0]
            current_sstart = df_scaffst.sstart.iloc[0]
            list_regions = []

            for I,r in df_scaffst.iterrows(): 
                if (r.qstart<previous_qstart):
                    list_regions.append([current_sstart,previous_send])
                    current_sstart = r.sstart
                         
                previous_qstart = r.qstart
                previous_send = r.send
            list_regions.append([current_sstart,previous_send])

            for r1, r2 in list_regions:
                if strand=="+": 
                    start, end = r1, r2
                else: 
                    start, end = r2, r1

                if end<=start:
                    raise ValueError("Start can't be after end")

                start = max([1, (start - 50000)])

                df_len = pandas.read_csv("/mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.lengths", sep=" ", names=["len", "scaffold"], header=None)
                len_scaffold = df_len[df_len.scaffold==scaffold]["len"].iloc[0]
                end = min([(end + 50000), len_scaffold])
                length_protein = end-start

                os.chdir("%s"%(ProtDir))

                HB_id = "%s_HB_%s_%s_%s"%(protein_name,scaffold,start,end)

                if not os.path.isfile("%s_Fastafetch.fa"%(HB_id)): 
                    run_cmd("fastafetch %s /mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.index %s > %s_Fastafetch.fa"%(HBGenome, scaffold,HB_id))

                if not os.path.isfile("%s_genomic.fa"%(HB_id)): 
                    run_cmd("fastasubseq %s_Fastafetch.fa %s %s > %s_genomic.fa" %(HB_id,start,length_protein,HB_id))
                
                if not os.path.isfile("%s_exonerate.nuc"%(HB_id)):
                    os.system("exonerate -m p2g --showtargetgff -q %s.fasta -t %s_genomic.fa | egrep -w exon > %s_exonerate.nuc"%(record_query.id, HB_id,HB_id))
                    exonerate_output = "%s_exonerate.nuc"%(HB_id)
                    nlines_exonerate = len(open(exonerate_output, "r").readlines())
                    if nlines_exonerate==0:
                        if strand=="+":
                            run_cmd("genewise -pep -pretty -cdna -gff %s %s_genomic.fa > %s_genewisepredprot.aa"%(query_protein,HB_id,HB_id))
                            
                        else:
                            run_cmd("genewise -pep -pretty -cdna -gff -trev %s %s_genomic.fa > %s_genewisepredprot.aa"%(query_protein,HB_id,HB_id))
                            
                    else:
                        if not os.path.isfile("%s_fastaseqfromGFF"%(HB_id)):
                            run_cmd("fastaseqfromGFF.pl %s_genomic.fa %s_exonerate.nuc > %s_fastaseqfromGFF"%(HB_id, HB_id, HB_id))

                        if not os.path.isfile("%s_predprot.fasta"%(HB_id)):
                            run_cmd("fastatranslate -F 1 %s_fastaseqfromGFF > %s_predprot.fasta"%(HB_id, HB_id))

                        if not os.path.isfile("%s_tcoffee"%(HB_id)):
                            run_cmd("t_coffee %s.fasta %s_predprot.fasta > %s_tcoffee.fa"%(record_query.id, HB_id, HB_id))
                    
                        if not os.path.isfile("%s/%s_predprot.fasta"%(PredictedDir,HB_id)):
                            run_cmd("cp %s/%s_predprot.fasta %s"%(ProtDir,HB_id,PredictedDir))

                        for record_HB in SeqIO.parse("%s/%s_predprot.fasta"%(PredictedDir,HB_id), "fasta"): 
                            if "*" in str(record_HB.seq):
                                HB_selenocistein = "U"
                            else:
                                HB_selenocistein="C"
                        
                        ZF_selenocistein = ""

                        for record_ZF in SeqIO.parse("%s/zf.fasta"%(ParentDir), "fasta"): 
                            record_ZF.id = record_ZF.description.split("(")[1].split(")")[0]
                            if record_ZF.id == protein_name:
                                if "U" in str(record_ZF.seq):
                                    ZF_selenocistein = "U"
                                else:
                                    ZF_selenocistein = "C"
                            if record_ZF.id != protein_name and ZF_selenocistein != "U" and ZF_selenocistein != "C":
                                ZF_selenocistein = "-"

                        df_HB=pandas.DataFrame([{
                                "Protein":protein_name,
                                "Scaffold":scaffold,
                                "Start":start,
                                "End":end,
                                "Strand":strand,
                                "Query_species":record_query.id.split("_")[1],
                                "HB_selenocistein":HB_selenocistein,
                                "ZF_selenocistein":ZF_selenocistein,
                            }])

                        df_total=df_total.append(df_HB)
                    
                os.chdir("%s"%(ParentDir))

df_total.to_csv('%s/selenocistein_table.csv'%ParentDir)
