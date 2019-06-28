#!/usr/bin/python
import argparse, sys, subprocess
from Bio import SeqIO
import multiprocessing as mp

def filterFasta(fastaFile, minLength, outFile):
    newFasta = []
    for record in SeqIO.parse(fastaFile,"fasta"):
        if len(record.seq) >= minLength :
            newFasta.append(record)
    SeqIO.write(newFasta,outFile,"fasta")


def createGFF(fastaFile, outGFF):
    out_gff = open(outGFF,"w")

    for rec in SeqIO.parse(fastaFile,"fasta"):
        line = [""]*9
        line[1] = "triniUmap"
        line[2]="CDS"
        line[5]= "."
        line[6] = "+"
        line[7]="."

        # specific for each gene
        line[0]=rec.name
        line[3]= "1"
        line[4]= str(len(rec.seq))
        line[8]= "ID="+rec.name+";Name="+rec.name
        out_gff.write( "\t".join(line) + "\n")

    out_gff.close()

def runHtSeq(cmdHtseq):
    subprocess.call(cmdHtseq,shell=True)
    return 0
    

def main():

    '''arguments parsing '''

    parser = argparse.ArgumentParser(description="Run full analysis")

    ### general arguments ###
    parser.add_argument('--out', help="Output name for the analysis, required", required=True)
    parser.add_argument('--umap',help="Set this option if you want to run Trinity on unmap reads", action='store_true')
    parser.add_argument('--single',help="Set this option if you have single-end reads", action='store_true')
    
    ### arguments bowtie2 ###        
    parser.add_argument('--ref', help="Reference PanGenome or Genome, required", required=True)
    parser.add_argument('--fq', help="Tab delimited file with paired fastq file names", required=True)
    parser.add_argument('--p', help="Number of threads to run bowtie2 and samtools, default 1", default=1)
    parser.add_argument('--index',help="Set this option if you do NOT want to index reference", action='store_false')

    ### arguments htseq-count ###
    parser.add_argument('--gff', help="Gff file describing pangenome or genome", required=True)
    parser.add_argument('--a', help="Minimum quality of the reads [0-10], default 10", default="10")
    parser.add_argument('--idattr', help="Feature in gff defining gene name, default Name", default="Name")
    parser.add_argument('--t', help="Feature to consider from gff, default CDS", default="CDS")
    parser.add_argument('--m', help="htseq-count sensibility, default intersection-nonempty", default="intersection-nonempty")
    parser.add_argument('--Proc', help="Number of processes to run htseq in parallell, default 1", default="1")
    parser.add_argument('--stranded', help="whether the data is from a strand-specific assay <yes/no/reverse>, default: no", default="no")

    ### arguments Trinity ###
    parser.add_argument('--max_memory', help="Maximum RAM allowed for Trinity, default 10G", default="10G")
    parser.add_argument('--CPU', help="Maximum CPU allowed for Trinity, default 5", default="5")
    parser.add_argument('--filt', help="Minimal length for Trinity scaffold filter, default 800", default="800")

    ### arguments ORFfinder ###
    parser.add_argument('--ml', help="Minimal length of the ORF (nt) that will detect ORFfinder, default 600", default="600")

    ### arguments cd-hit-est-2d ###
    parser.add_argument('--c', help="sequence identity threshold, this is the default cd-hit's global sequence identity calculated as: \
    number of identical amino acids or bases in alignment divided by the full length of the shorter sequence, default 0.8", default="0.8")
    parser.add_argument('--s2', help="length difference cutoff by default, seqs in db1 >= seqs in db2 in a same cluster\
    if set to 0.9, seqs in db1 may just >= 90% seqs in db2, default 0.8", default="0.8")
    parser.add_argument('--i', help="Fasta file with all CDS sequences annotated in REF (in case of using complete genome), default REF", default="NO")

    
    parser.set_defaults()
    args = parser.parse_args()

    print args

    ##################################################################################
    ################################# BASIC RUN ######################################
    
    ''' Index Reference Genome '''
    if args.index :
        subprocess.call(["bowtie2-build", args.ref, args.ref])

    ''' Create Directory for the full analysis --out and
        subdirectories structure'''
    subprocess.call(["mkdir", args.out])
    subprocess.call(["mkdir", args.out+"/alignments"])
    subprocess.call(["mkdir", args.out+"/counts"])
    if args.umap :
        subprocess.call(["mkdir", args.out+"/Trinity_umap"])
    

    ''' Parse the input fq list file '''

    fqFilesList = []
    inputFq = open(args.fq,"r")    
    for line in inputFq:
        fqFilesList.append(line.rstrip().split("\t"))


    ''' Run analysis on each sample '''
    
    print "PANRNA : BEGIN PIPELINE.... BASE"
    for sample in fqFilesList:

        # run Name is whatever is before the first point and after last / in first read #
        runName = sample[0].split("/")[-1].split(".")[0]        

        ### map with bowtie2 and sort bam by name ###
        outMap = args.out+"/alignments/"+runName

        if args.single:
            print "PANRNA : BASE RUN BOWTIE2 ON SINGLE END"
            subprocess.call(["bowtie2", "-x", args.ref, "-U", sample[0],
                             "-N", "1", "-q","--local", "-p", args.p,
                             "-S", outMap+".sam"])
            
        else :
            print "PANRNA : BASE RUN BOWTIE2 ON PAIRED END"
            subprocess.call(["bowtie2", "-x", args.ref, "-1", sample[0], "-2", sample[1],
                             "-N", "1", "-q","--local", "-p", args.p,
                             "-S", outMap+".sam"])
            
            

        subprocess.call(["samtools", "view", "-b","-@", args.p, "-S",
                         "-o", outMap+".bam",
                         outMap+".sam"])

        subprocess.call(["samtools", "sort", "-n","-@", args.p,
                         "-o", outMap+".sorted.bname.bam",
                         outMap+".bam"])

        # delete .sam and unsorted .bam files #
        subprocess.call(["rm", outMap+".sam"])
        subprocess.call(["rm", outMap+".bam"])

        


    ### parallell counting with htseq-count###
    htseqCmds = []
    for sample in fqFilesList:
        # run Name is whatever is before the first point and after last / in first read #
        runName = sample[0].split("/")[-1].split(".")[0]      
        
        outMap = args.out+"/alignments/"+runName
        cmdHtseq = "samtools view " + outMap +".sorted.bname.bam" + \
                   " | htseq-count -r name -a " + args.a + \
                   " --idattr "+ args.idattr +" -t "+ args.t +" -m "+ args.m +" --stranded=" + args.stranded +" - "+ \
                   args.gff +" > " + args.out+"/counts/"+runName+".counts.table"
        htseqCmds.append(cmdHtseq)

    pool = mp.Pool(processes=int(args.Proc))
    pool.map(runHtSeq, htseqCmds)

        
    print "PANRNA : END BASE PIPELINE"
    ##################################################################################
    ############################## Trinity Unmapped ##################################

    if args.umap :
        print "PANRNA : BEGIN UNMAPPED ANALYSIS"
        ''' Create Directories structure'''
        subprocess.call(["mkdir", args.out+"/Trinity_umap/alignments_Umap"])
        subprocess.call(["mkdir", args.out+"/Trinity_umap/reads_Umap"])
        subprocess.call(["mkdir", args.out+"/Trinity_umap/trinity_out_dir"])
        subprocess.call(["mkdir", args.out+"/Trinity_umap/alignments_trinity"])
        subprocess.call(["mkdir", args.out+"/Trinity_umap/counts_trinity"])

        ''' Extract unmapped reads from each file '''
        for sample in fqFilesList :
            runName = sample[0].split("/")[-1].split(".")[0]
            outDirSam = args.out+"/Trinity_umap/alignments_Umap/"+runName
            outDirFq = args.out+"/Trinity_umap/reads_Umap/"+runName
            inBam = args.out+"/alignments/"+runName+".sorted.bname.bam"

            ### Create Sam of unmapped ###
            outSam=open(outDirSam+".umap.sam","w")
            subprocess.call(["samtools", "view","-@", args.p,"-f4", inBam], stdout = outSam)
            outSam.close()

            ### extract ids of reads, and reads from it ###
            extrIDsCmd = "cut -f1 "+outDirSam+".umap.sam"+" | sort | uniq > "+outDirSam+".umap_id.txt"
            subprocess.call(extrIDsCmd, shell=True)

            # R1 #
            print "PANRNA : EXTRACT R1 READS FROM UNMAPPED"
            outFqR1=open(outDirFq+".R1.fastq","w")
            subprocess.call(["/usr/local/bin/seqtk/seqtk", "subseq",sample[0], outDirSam+".umap_id.txt"], stdout = outFqR1)
            outFqR1.close()
            
            # R2 #
            if not args.single:
                print "PANRNA : EXTRACT R2 READS FROM UNMAPPED"
                outFqR2=open(outDirFq+".R2.fastq","w")
                subprocess.call(["/usr/local/bin/seqtk/seqtk", "subseq",sample[1], outDirSam+".umap_id.txt"], stdout = outFqR2)
                outFqR2.close()

            # compress files... #            
            subprocess.call(["samtools", "view", "-b", "-S","-@", args.p, "-o", outDirSam+".umap.bam",outDirSam+".umap.sam"])

            subprocess.call(["samtools", "sort", "-n","-@", args.p, "-o", outDirSam+".sorted.bname.umap.bam", outDirSam+".umap.bam"])

            # delete .sam and unsorted .bam files #
            subprocess.call(["rm", outDirSam+".umap.sam"])
            subprocess.call(["rm", outDirSam+".umap.bam"])

        
        ''' Run Trinity on concatenated fastq '''

        ### concatenate reads files for each pair ###
        print "PANRNA : CONCATENATE R1 READS FROM UNMAPPED"
        cmdCatR1="cat "+args.out+"/Trinity_umap/reads_Umap/"+"*.R1.fastq"+" > "+args.out+"/Trinity_umap/reads_Umap/"+"allNm.R1.fastq"
        subprocess.call(cmdCatR1,shell=True)      

        if not args.single :
            print "PANRNA : CONCATENATE R2 READS FROM UNMAPPED"
            cmdCatR2="cat "+args.out+"/Trinity_umap/reads_Umap/"+"*.R2.fastq"+" > "+args.out+"/Trinity_umap/reads_Umap/"+"allNm.R2.fastq"
            subprocess.call(cmdCatR2,shell=True)

        ### compress read files ###        
        cmdCompFq = "gzip "+args.out+"/Trinity_umap/reads_Umap/"+"*.fastq"
        subprocess.call(cmdCompFq,shell=True)

        ### RUN TRINITY on concat reads ###
        if args.single :
            print "PANRNA : RUN TRINITY IN SINGLE END READS MODE"
            readTriniR1 = args.out+"/Trinity_umap/reads_Umap/"+"allNm.R1.fastq.gz"
            subprocess.call(["Trinity", "--seqType", "fq", "--max_memory", args.max_memory,
                         "--single", readTriniR1, "--CPU", args.CPU,
                         "--jaccard_clip", "--trimmomatic", "--output", args.out+"/Trinity_umap/trinity_out_dir"])
            
        else :
            print "PANRNA : RUN TRINITY IN PAIRED END READS MODE"
            readTriniR1 = args.out+"/Trinity_umap/reads_Umap/"+"allNm.R1.fastq.gz"
            readTriniR2 = args.out+"/Trinity_umap/reads_Umap/"+"allNm.R2.fastq.gz"

            subprocess.call(["Trinity", "--seqType", "fq", "--max_memory", args.max_memory,
                             "--left", readTriniR1, "--right",  readTriniR2, "--CPU", args.CPU,
                             "--jaccard_clip", "--trimmomatic", "--output", args.out+"/Trinity_umap/trinity_out_dir"])


        ### Filter Trinity output by contig length ###
        outFilteredFasta = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.fasta"
        filterFasta(args.out+"/Trinity_umap/trinity_out_dir/Trinity.fasta", int(args.filt),outFilteredFasta)

        
        ### Find ORFs in the Trinity scaffolds ###
        print "PANRNA : RUN ORFfinder with minimal length : "+ args.ml
        outFilteredCDS = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.CDS.fasta"
        subprocess.call(["ORFfinder", "-in", outFilteredFasta, "-n","true", "-ml", args.ml, "-s", "0","-outfmt", "1",
                         "-out", outFilteredCDS])

        ### Cluster the CDS against the Pangenome with cd-hit-est-2d and create gff###
        outClusteredCDS = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.CDS.clustered.fasta"

        # if no CDS fasta file use ref
        if args.i == "NO":
            subprocess.call(["cd-hit-est-2d", "-d", "0", "-c",args.c, "-s2", args.s2, "-T", args.p,"-i", args.ref,
                             "-i2", outFilteredCDS, "-o", outClusteredCDS])
        else :
            subprocess.call(["cd-hit-est-2d", "-d", "0", "-c",args.c, "-s2", args.s2, "-T", args.p,"-i", args.i,
                             "-i2", outFilteredCDS, "-o", outClusteredCDS])

        outClusteredGFF = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.CDS.clustered.gff"
        createGFF(outClusteredCDS, outClusteredGFF)

        
        ''' Do a mapping and htseq-count for each sample on contigs from Trinity '''

        ### Index the Trinity fasta file ###
        fastaTrinity = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.CDS.clustered.fasta"
        outFilteredGFF = args.out+"/Trinity_umap/trinity_out_dir/Trinity.filtered.CDS.clustered.gff"

        subprocess.call(["bowtie2-build", fastaTrinity, fastaTrinity])
        
        for sample in fqFilesList :
            runName = sample[0].split("/")[-1].split(".")[0]
            outMap = args.out+"/Trinity_umap/alignments_trinity/"+runName
            fastqR1um = args.out+"/Trinity_umap/reads_Umap/"+runName+".R1.fastq.gz"            
            if not args.single :
                
                fastqR2um = args.out+"/Trinity_umap/reads_Umap/"+runName+".R2.fastq.gz"
                print "PANRNA : PAIRED END READS FOR BOWTIE IN TRINITY PART R2 :" + fastqR2um
            

            ### Mapping bowtie2 ###
            if args.single :
                print "PANRNA : RUN BOWTIE ON TRINITY WITH SINGLE END READS"
                subprocess.call(["bowtie2", "-x", fastaTrinity, "-U", fastqR1um, 
                             "-N", "1", "-q","--local", "-p", args.p,
                             "-S", outMap+".sam"])
            else:
                print "PANRNA : RUN BOWTIE ON TRINITY WITH PAIRED END READS"
                subprocess.call(["bowtie2", "-x", fastaTrinity, "-1", fastqR1um, "-2", fastqR2um,
                                 "-N", "1", "-q","--local", "-p", args.p,
                                 "-S", outMap+".sam"])            
            

            subprocess.call(["samtools", "view", "-b", "-S","-@", args.p,
                         "-o", outMap+".bam",
                         outMap+".sam"])

            subprocess.call(["samtools", "sort", "-n","-@", args.p,
                         "-o", outMap+".trinityF.sorted.bname.bam",
                         outMap+".bam"])

            # delete .sam and unsorted .bam files #
            subprocess.call(["rm", outMap+".sam"])
            subprocess.call(["rm", outMap+".bam"])

            

        ### parallell counting with htseq-count###
        htseqCmds = []
        for sample in fqFilesList:
            # run Name is whatever is before the first point and after last / in first read #
            runName = sample[0].split("/")[-1].split(".")[0]
            outMap = args.out+"/Trinity_umap/alignments_trinity/"+runName
            
            cmdHtseq = "samtools view " + outMap +".trinityF.sorted.bname.bam" + \
                           " | htseq-count -r name -a " + args.a + \
                           " --idattr "+ "Name" +" -t "+ "CDS" +" -m "+ args.m +" --stranded="+ args.stranded +" - "+ \
                           outFilteredGFF +" > " + args.out+"/Trinity_umap/counts_trinity/"+runName+".trinityF.counts.table"
            htseqCmds.append(cmdHtseq)

        pool = mp.Pool(processes=int(args.Proc))
        pool.map(runHtSeq, htseqCmds)

    print "PANRNA : ALL ANALYSIS FINISHED"            
        
    

if __name__ == "__main__":
    main()
