from Bio import SeqIO #import SeqIO from Biopython module
from Bio.Seq import Seq #import Seq from Bio.Seq module
import os #import operating system module
import pandas #import pandas module

exonerate_output = 1 #exonerate file is full

ParentDir = "/run/user/1400982974/gvfs/smb-share:server=fs-aules.ds.upf.edu,share=public/20428/Treball/Haplochromis_burtoni_02_12" #variable that contains the main directory path, where query.txt file is found. it has to be changed every time for the directory of interest

HBGenome = "/mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.fa" #variable that contains the path where you can find the fasta file containing the genome of our species, Haplochromis burtoni

def run_cmd(cmd): #function that will print an error when an os.system command does not work
    out_stat = os.system(cmd) #os.system always prints an out status when executed, which can be 0 (has not worked) or different to 0 (has worked)
    if out_stat!=0: #conditional statement that is executed when os.system command does not work
        raise ValueError("The command '%s' did not work"%cmd) #raise an error and print the message in the terminal

PredictedDir = "%s/Predicted_proteins"%(ParentDir) #variable that contains the path to the Predicted_proteins folder
if not os.path.isdir(PredictedDir): #conditonal statement that is only executed when the specified directory does not exist
    run_cmd("mkdir %s"%PredictedDir) #command that executes the making of the specified directory

df_total=pandas.DataFrame() #variable that contains an empty dataframe. It will be used later to append the data for each predicted protein

for record_query in SeqIO.parse(ParentDir+"/query.txt", "fasta"): #treat the query.txt file as a fasta file and parse its contents for each fasta sequence, which starts with a >
    record_query.id = record_query.description.split("(")[1].replace(") # ","_").strip() #variable that contains the name of the protein + "_" + the query species . split allows you to split your text in the specified area, replace allows you to replace the specified text with another, and strip deletes any blank spaces found in the end of the variable.
    protein_name = record_query.id.split("_")[0].strip() #variable that just contains the name of the protein

    ProtDir = "%s/%s"%(ParentDir, record_query.id) #variable that stores the name of the folder that will contain all the documents for the predicted proteins of a single query protein
    if not os.path.isdir(ProtDir): #conditional statement that is only executed when the specified directory does not exist
        run_cmd("mkdir %s"%ProtDir) #command that executes the making of the specified directory
    
    query_protein = "%s/%s.fasta"%(ProtDir, record_query.id) #variable that contains the path and filename of a single query protein
    record_query.description = "" #variable that empties the description of the fasta file
    record_query.name = "" #variable that empties the name of the fasta file
    
    record_query.seq = Seq(str(record_query.seq).replace("U","X")) #variable that stores the sequence of each single query protein where U are replaced by X. IT has to be specified as string to be able to work, if not, replace does not work

    SeqIO.write([record_query], query_protein, "fasta") #create a new fasta file containing all the record modifications stored in the aforementioned record variables

    blast_output = "%s/%s.blast"%(ProtDir, record_query.id) #variable that contains the path and filename of a single query protein blast
    
    if not os.path.isfile(blast_output): #conditional statement that is only executed when the specified file does not exist
        run_cmd('tblastn -evalue 0.001 -query %s -db %s -out %s -outfmt "6 qstart qend sseqid sstart send sstrand evalue"'%(query_protein,HBGenome,blast_output)) #run blast and format it as a table

    df_blast = pandas.read_csv(blast_output, sep="\t", names=["qstart", "qend", "sseqid", "sstart", "send", "sstrand", "evalue"], header=None) #variable that stores, as a table, the results of the tblastn
    df_blast["sstrand"] = (df_blast.sstart < df_blast.send).map({True:"+", False:"-"}) #value that substitutes the preexisting sstrand column (which was seen to be N/A) with a + or - sign depending on whether the beginning of a particular hit is higher or lower than the end (if it is higher it will be found in the - strand, whereas if it is lower, it is found in the + strand)
    for scaffold in set(df_blast.sseqid): #each iteration will go through each value (row) that can be taken by the sseqid column. set makes a list with each non-repeated value; as such, no id that has previously been analysed will be analysed again
        for strand in set(df_blast.sstrand): #same as the previous command, but instead of looking at the sseqid column, looking at the sstrand
            df_scaffst = df_blast[(df_blast['sseqid'] == scaffold)&(df_blast['sstrand'] == strand)] #new dataframe that only contains one of the possible sseqid values and sstrand values

            if len(df_scaffst)==0: #set command is random, so it could take a scaffold that does not have a +/- strand (dataframe is empty), signaling an error and stopping the program
                continue #tells the program to continue when the dataframe is empty
            if strand =="+": #conditional statement that is only executed when the strand is +
                df_scaffst = df_scaffst.sort_values(by="sstart", ascending=True) #substitutes the preexisting values of the dataframe by sorting the sstart column in ascending order 
            else: #conditional statement that is only executed when the strand is NOT + (aka, -)
                df_scaffst = df_scaffst.sort_values(by="sstart", ascending=False) #substitutes the preexisting values of the dataframe by sorting the sstart column in descascending order 

            previous_qstart = 0 
            previous_send = df_scaffst.send.iloc[0] #variable that stores the value of the first row from the send column
            current_sstart = df_scaffst.sstart.iloc[0] #variable that stores the value of the first row from the sstart column
            list_regions = [] #new list variable

            for I,r in df_scaffst.iterrows(): #each iteration will look at the different values of each column for each row.
                if (r.qstart<previous_qstart): #conditional statement that is only executed when the value of the qstart of the current row is smaller than the value of the qstart of the previous row (which means that the current hit belongs to a different protein)
                    list_regions.append([current_sstart,previous_send]) #add current start and end values for the subject to the list_regions variable (list of lists)
                    current_sstart = r.sstart #change the preexisting subject start to the current one
                         
                previous_qstart = r.qstart #for every iteration, change the variable value
                previous_send = r.send #for every iteration, change the variable value
            list_regions.append([current_sstart,previous_send]) #add current start and end values for the subject to the list_regions variable (if this command was not here, the last prediction would not be recorded)

            for r1, r2 in list_regions: #each interation will go through the start and end values of each position of the list_regions vector
                if strand=="+": #conditional statement that is only executed when the strand is +
                    start, end = r1, r2 #start = r1, end = r2
                else: #conditional statement that is only executed when the strand is NOT + (aka, -)
                    start, end = r2, r1 #start = r2, end = r1

                if end<=start: #conditional statement that is only executed when end value is smaller or equal to start
                    raise ValueError("Start can't be after end") #raise and print error in the terminal (stops the program)

                start = max([1, (start - 50000)]) #start value that will be used for the fastasubseq. we are choosing the maximum value between 1 and the substraction to the previous start (these will be genomic coordinates, which cannot be negative)

                df_len = pandas.read_csv("/mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.lengths", sep=" ", names=["len", "scaffold"], header=None) #variable that contains the lengths and scaffold id of all the Haplocgromis burtoni genome
                len_scaffold = df_len[df_len.scaffold==scaffold]["len"].iloc[0] #variable that stores the length value of the scaffold of interest
                end = min([(end + 50000), len_scaffold]) #end value that will be used for the fastasubseq. we are choosing the minimum value between the total length of the scaffold and the addition  to the previous end (these will be genomic coordinates, which cannot be higher than the total length of the scaffold)
                length_protein = end-start #variable that contains the length of the start and end values for the fastasubseq (highest length than the real one)

                os.chdir("%s"%(ProtDir)) #change directory

                HB_id = "%s_HB_%s_%s_%s"%(protein_name,scaffold,start,end) #variable that contains the name of all the documents that will now be created

                if not os.path.isfile("%s_Fastafetch.fa"%(HB_id)): #conditional statement that is only executed when the specified file does not exist 
                    run_cmd("fastafetch %s /mnt/NFS_UPF/soft/genomes/2021/Haplochromis_burtoni/genome.index %s > %s_Fastafetch.fa"%(HBGenome, scaffold,HB_id)) #run fastafetch

                if not os.path.isfile("%s_genomic.fa"%(HB_id)): #conditional statement that is only executed when the specified file does not exist
                    run_cmd("fastasubseq %s_Fastafetch.fa %s %s > %s_genomic.fa" %(HB_id,start,length_protein,HB_id)) #run fastasubseq
                
                if not os.path.isfile("%s_exonerate.nuc"%(HB_id)): #conditional statement that is only executed when the specified file does not exist
                    os.system("exonerate -m p2g --showtargetgff -q %s.fasta -t %s_genomic.fa | egrep -w exon > %s_exonerate.nuc"%(record_query.id, HB_id,HB_id)) #run exonerate
                    exonerate_output = "%s_exonerate.nuc"%(HB_id) #variable that stores the exonerate filename
                    nlines_exonerate = len(open(exonerate_output, "r").readlines()) #variable that stores the number of lines found in the file
                    if nlines_exonerate==0: #conditional statement that is only executed when there are no lines in the variable
                        if strand=="+": #conditional statement that is only executed when the strand is +
                            run_cmd("genewise -pep -pretty -cdna -gff %s %s_genomic.fa > %s_genewisepredprot.aa"%(query_protein,HB_id,HB_id)) #run genewise
                            
                        else: #conditional statement that is only executed when the strand is NOT + (aka, -)
                            run_cmd("genewise -pep -pretty -cdna -gff -trev %s %s_genomic.fa > %s_genewisepredprot.aa"%(query_protein,HB_id,HB_id)) #run genewise
                            
                    else: #conditional statement that is only executed when there are  lines in the variable
                        if not os.path.isfile("%s_fastaseqfromGFF"%(HB_id)): #conditional statement that is only executed when the specified file does not exist
                            run_cmd("fastaseqfromGFF.pl %s_genomic.fa %s_exonerate.nuc > %s_fastaseqfromGFF"%(HB_id, HB_id, HB_id)) #run fastaseqfromGFF

                        if not os.path.isfile("%s_predprot.fasta"%(HB_id)): #conditional statement that is only executed when the specified file does not exist
                            run_cmd("fastatranslate -F 1 %s_fastaseqfromGFF > %s_predprot.fasta"%(HB_id, HB_id)) #run fastatranslate

                        if not os.path.isfile("%s_tcoffee"%(HB_id)): #conditional statement that is only executed when the specified file does not exist
                            run_cmd("t_coffee %s.fasta %s_predprot.fasta > %s_tcoffee.fa"%(record_query.id, HB_id, HB_id)) #run t-coffee
                    
                        if not os.path.isfile("%s/%s_predprot.fasta"%(PredictedDir,HB_id)): #conditional statement that is only executed when the specified file does not exist
                            run_cmd("cp %s/%s_predprot.fasta %s"%(ProtDir,HB_id,PredictedDir)) #copy the file containing the aminoacid sequence of the predicted protein to the Predicted_proteins folder

                        for record_HB in SeqIO.parse("%s/%s_predprot.fasta"%(PredictedDir,HB_id), "fasta"): #treat the predprot.fasta file as a fasta file and parse its contents for each fasta sequence, which starts with a >
                            if "*" in str(record_HB.seq): #conditional statement that is only executed when the sequence of the fasta file contains a "*"
                                HB_selenocistein = "U" #variable that stores a "U"
                            else: #conditional statement that is only executed when the sequence of the fasta does NOT file contains a "*"
                                HB_selenocistein="C" #variable that stores a "C"
                        
                        ZF_selenocistein = "" #variable that initialises in every iteration

                        for record_ZF in SeqIO.parse("%s/zf.fasta"%(ParentDir), "fasta"): #treat the query.txt file as a fasta file and parse its contents for each fasta sequence, which starts with a >
                            record_ZF.id = record_ZF.description.split("(")[1].split(")")[0] #substitute the preexisting id for a simplified version that only contains the name of the protein
                            if record_ZF.id == protein_name: #conditional statement that is only executed when the name of the zf protein is the same as the name of the query protein (aka, it is found in zebrafish)
                                if "U" in str(record_ZF.seq): #conditional statement that is only executed when the sequence of the fasta file contains a "U"
                                    ZF_selenocistein = "U" #variable that stores a "U"
                                else: #conditional statement that is only executed when the sequence of the fasta file does NOT contains a "U"
                                    ZF_selenocistein = "C" #variable that stores a "C"
                            if record_ZF.id != protein_name and ZF_selenocistein != "U" and ZF_selenocistein != "C": #conditional formating that is only executed when the name of the zf protein is not the same as the name of the query protein (aka, it is not found in zebrafish) AND when ZF_selenocistein is not a "U" and a "C" (if this last part were not here, the program would overwrite the U/C written after going into line 133)
                                ZF_selenocistein = "-" #variable that stores a "-"

                        df_HB=pandas.DataFrame([{ #make a new dataframe with the following column names (in orange) and values (blue)
                                "Protein":protein_name,
                                "Scaffold":scaffold,
                                "Start":start,
                                "End":end,
                                "Strand":strand,
                                "Query_species":record_query.id.split("_")[1],
                                "HB_selenocistein":HB_selenocistein,
                                "ZF_selenocistein":ZF_selenocistein,
                            }])

                        df_total=df_total.append(df_HB) #append this dataframe containing information for each single predicted protein to the previous existing data
                    
                os.chdir("%s"%(ParentDir)) #change directory

df_total.to_csv('%s/selenocistein_table.csv'%ParentDir) #save the dataframe containing the data for all predicted proteins as a csv file
