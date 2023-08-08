
import matplotlib
matplotlib.use("Agg")
from pylab import *
import re
import pysam
from Bio import  SeqIO
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,BboxConnectorPatch
import matplotlib.font_manager as fm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pandas as pd
import os


# DATA_DIR = "../data/riboseq_patrick/clavuligerus/"
DATA_DIR = "../data/riboseq_patrick/griseus/"
IMG_DIR = "../images/riboseq_patrick/"

mRNA_sequences = os.path.join(DATA_DIR, "AP009493.1.fasta")
gff3_file = os.path.join(DATA_DIR, "AP009493.1.gff3")


def triple_frames(list_profile) :
    zero_f = [0 for coordinate66 in range(len(list_profile))]
    one_f = zero_f[:]
    two_f = zero_f[:]

    for n70 in range(len(list_profile)) :
        if n70 %3 == 0 :
            zero_f[n70] = list_profile[n70]
        elif  n70 %3 == 1 : one_f[n70] = list_profile[n70]
        else :two_f[n70] = list_profile[n70]
    return zero_f, one_f, two_f

def read_bamfile(aligments_A1,chromosome,profiles_dict,riboseq, strand,cdsstart, cdsend):
   
    coding_codordinates2 = [n for n in range( cdsstart, cdsend )]
    coding_codordinates2set = set(coding_codordinates2)
    if strand == "-":
            coding_codordinates2.reverse()

    all_reads        = aligments_A1.fetch(chromosome, cdsstart,cdsend)

    offset = 0
    for read in all_reads :
            readlen = read.qlen
            if strand == "+" :
                if read.is_reverse : continue
                cpositions = read.positions
                cpositions.sort()
                coordinate = cpositions[-1 -11]                
                
                if riboseq :
                    if coordinate in coding_codordinates2set:
                        profiles_dict[coding_codordinates2.index(coordinate)] += 1
                 
                
                else:
                    coordinate = cpositions[0]                
                    if coordinate in coding_codordinates2set:
                        profiles_dict[coding_codordinates2.index(coordinate)] += 1
            else:
                if not read.is_reverse : continue
                cpositions = read.positions
                cpositions.sort()
                coordinate = cpositions[-1-offset]
                coordinate = cpositions[0 +  11]
                if riboseq :
                        if coordinate in coding_codordinates2set:
                            profiles_dict[coding_codordinates2.index(coordinate)] += 1
                        else:
                            coordinate = cpositions[0] 
                        if coordinate in coding_codordinates2set:
                            profiles_dict[coding_codordinates2.index(coordinate)] += 1 
        
    return profiles_dict


    
    

def axesprep(axespanel):
    axespanel.spines['top'].set_visible(False)
    axespanel.spines['right'].set_visible(False)
    axespanel.spines['bottom'].set_visible(False)
    axespanel.xaxis.set_major_locator(MaxNLocator(2))
    axespanel.yaxis.set_major_locator(MaxNLocator(2,prune = "both"))
    axespanel.get_xaxis().tick_bottom()
    axespanel.get_yaxis().tick_left()
    axespanel.set_xticks([],[])


rcParams["xtick.direction"]        = "out"
rcParams["ytick.direction"]        = "out"
rcParams["legend.fontsize"]        = 15
rcParams["ytick.labelsize"]        = 9
rcParams["xtick.labelsize"]        = 9
rcParams['font.size']         =15
rcParams["axes.titlesize"]     = 15
rcParams["legend.frameon"]         = 0
rcParams["axes.axisbelow"]      = False
rcParams["xtick.major.pad"]     = 2.
rcParams["ytick.major.pad"]     = 2
rcParams["xtick.major.size"]    = 2.
rcParams["ytick.major.size"]     = 2
rcParams['axes.linewidth']     = 1.5


rcParams["legend.borderpad"]    = 0.01
rcParams["legend.labelspacing"] = 0.05
rcParams["legend.columnspacing"]= 0.13
rcParams["legend.borderaxespad"]=0.15
rcParams["legend.handlelength"]    = 1


in_seq_handle = open(mRNA_sequences)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()





fig = plt.figure(figsize = (8,8))

gs = gridspec.GridSpec(56, 5)
gs.update(left=0.08, right=0.95, wspace=5.50,top = 0.91,hspace = 0.9,bottom = 0.04)

gs2 = gridspec.GridSpec(56, 5)
gs2.update(left=0.05, right=1.99, wspace=5.60,top = 0.91,hspace = 0.9,bottom = 0.04)






#Using all the alignments to the 3 replicates with unique and reads that mapped up to 3 places
#CONTROL RNA
NormGlu1_rna_file = pysam.Samfile(os.path.join(DATA_DIR, "SRR10212831_SgriRibo_E1_part.bam"), 'rb')
NormGlu1_rna_file2 = pysam.Samfile(os.path.join(DATA_DIR, "SRR10212833_SgriRibo_T1_part.bam"), 'rb')
NormGlu1_rna_file3 = pysam.Samfile(os.path.join(DATA_DIR, "SRR10212835_SgriRibo_L1_part.bam"), 'rb')
NormGlu1_rna_file4 = pysam.Samfile(os.path.join(DATA_DIR, "SRR10212837_SgriRibo_S1_part.bam"), 'rb')


panel1 = 0
panel2 = 0
panel3 = 0
panel1_old = 0
panel2_old = 0
panel4 =9



for index,mode in enumerate(["1"]):
    list_genes = []
    list_transcript_check = []

    
    infileopen = open(gff3_file)
    for line in infileopen:
        if line[0] == "#": continue
        if len(line.split("\t")) != 9 : continue
        if line.split("\t")[2] != "gene": continue
        name = ""
        info = line.split("\t")[8]

        data =  info.split(";")
        for dd in data:
            if "Name=" in dd:
                name = dd
        cdsstart, cdsend =  int(line.split("\t")[3])-1-450,int(line.split("\t")[4])+450
        strand = line.split("\t")[6]
        if name != "Name=adpA":
            continue
    
        if strand == "-":
            sequence = seq_dict["AP009493.1"][cdsstart:cdsend]
            sequence2= str(sequence.reverse_complement().seq)
        elif strand == "+":
            sequence2 = str(seq_dict["AP009493.1"][cdsstart:cdsend].seq)
    
        
        riboseq = 1
        profile_rpf_0a    = [0 for n in range(cdsstart,cdsend )]
        profile_rpf_0a = read_bamfile(NormGlu1_rna_file,"AP009493.1",profile_rpf_0a,riboseq, strand,cdsstart, cdsend)
        if sum(profile_rpf_0a) == 0: continue
        profile_rpf_0a_normalized = profile_rpf_0a
    
    
        transcript = name
        transcript_len = {}
        transcript_len[transcript] = len(profile_rpf_0a)

        profile_rpf_0a = np.array(profile_rpf_0a)/1
        transcript_seq =  str(sequence2)
        seq_line = list(transcript_seq) 
        xlimit = len(profile_rpf_0a)
        
        black_profile = [0 for n in range(len(profile_rpf_0a))]
        len_line = [n+1 for n in range(len(profile_rpf_0a))]
        
        ax1 = plt.subplot(gs[panel1:panel1+8,panel3:panel3+5])
        ax1.plot(profile_rpf_0a_normalized, color = "RoyalBlue", zorder=1)
        ax1.bar(np.arange(len(profile_rpf_0a_normalized)), profile_rpf_0a_normalized, width = 1, color = "RoyalBlue", ec = "RoyalBlue",zorder=2)
        ax1.plot(black_profile,color = "black", zorder=4)
        ax1.set_xlim(0,len(profile_rpf_0a))
        axesprep(ax1)
        xtis = ax1.get_yticks()
        ax1.set_title("Early-exponential phase")
        ylimit = max(profile_rpf_0a_normalized)
        try :
            intiation = 450
            stop =  transcript_len[transcript] - 450
            axvspan(intiation,stop, color = "gold", alpha = .5)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop+1])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'gold', 'alpha': 0.5, 'pad': 5})
            
            intiation, stop = 289,451
            intiation, stop =388,451            

            axvspan(intiation,stop, color = "red", alpha = .3)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*10)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
        except :
            noncoding_reads = sum(profile_rpf_0a_normalized)
            noncoding_reads = str(int(noncoding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, noncoding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})

        
        panel1 += 10
        profile_rpf_0a    = [0 for n in range(cdsstart,cdsend )]
        profile_rpf_0a = read_bamfile(NormGlu1_rna_file2,"AP009493.1",profile_rpf_0a,riboseq, strand,cdsstart, cdsend)
        if sum(profile_rpf_0a) == 0: continue
        profile_rpf_0a_normalized = profile_rpf_0a
    
    
        transcript = name
        transcript_len = {}
        transcript_len[transcript] = len(profile_rpf_0a)

        profile_rpf_0a = np.array(profile_rpf_0a)/1
        transcript_seq =  str(sequence2)
        seq_line = list(transcript_seq) 
        xlimit = len(profile_rpf_0a)
        
        black_profile = [0 for n in range(len(profile_rpf_0a))]
        len_line = [n+1 for n in range(len(profile_rpf_0a))]
        
        ax1 = plt.subplot(gs[panel1:panel1+8,panel3:panel3+5])
        ax1.plot(profile_rpf_0a_normalized, color = "RoyalBlue", zorder=1)
        ax1.bar(np.arange(len(profile_rpf_0a_normalized)), profile_rpf_0a_normalized, width = 1, color = "RoyalBlue", ec = "RoyalBlue",zorder=2)
        ax1.plot(black_profile,color = "black", zorder=4)
        ax1.set_xlim(0,len(profile_rpf_0a))
        axesprep(ax1)
        #ax1.set_ylabel("Rep1_rpf",fontsize = 10,color='white',backgroundcolor="darkred")
        xtis = ax1.get_yticks()
        ax1.set_title("Transition phase")
        ylimit = max(profile_rpf_0a_normalized)
        try :
            intiation = 450
            stop =  transcript_len[transcript] - 450
            axvspan(intiation,stop, color = "gold", alpha = .5)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop+1])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'gold', 'alpha': 0.5, 'pad': 5})
            
            intiation, stop = 289,451
            intiation, stop =388,451            

            axvspan(intiation,stop, color = "red", alpha = .3)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*10)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
        except :
            noncoding_reads = sum(profile_rpf_0a_normalized)
            noncoding_reads = str(int(noncoding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, noncoding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})



        panel1 += 10
        profile_rpf_0a    = [0 for n in range(cdsstart,cdsend )]
        profile_rpf_0a = read_bamfile(NormGlu1_rna_file3,"AP009493.1",profile_rpf_0a,riboseq, strand,cdsstart, cdsend)
        if sum(profile_rpf_0a) == 0: continue
        profile_rpf_0a_normalized = profile_rpf_0a
        
        
        transcript = name
        transcript_len = {}
        transcript_len[transcript] = len(profile_rpf_0a)

        profile_rpf_0a = np.array(profile_rpf_0a)/1
        transcript_seq =  str(sequence2)
        seq_line = list(transcript_seq) 
        xlimit = len(profile_rpf_0a)
        
        black_profile = [0 for n in range(len(profile_rpf_0a))]
        len_line = [n+1 for n in range(len(profile_rpf_0a))]
        
        ax1 = plt.subplot(gs[panel1:panel1+8,panel3:panel3+5])
        ax1.plot(profile_rpf_0a_normalized, color = "RoyalBlue", zorder=1)
        ax1.bar(np.arange(len(profile_rpf_0a_normalized)), profile_rpf_0a_normalized, width = 1, color = "RoyalBlue", ec = "RoyalBlue",zorder=2)
        ax1.plot(black_profile,color = "black", zorder=4)
        ax1.set_xlim(0,len(profile_rpf_0a))
        
        axesprep(ax1)
        #ax1.set_ylabel("Rep1_rpf",fontsize = 10,color='white',backgroundcolor="darkred")
        xtis = ax1.get_yticks()
        ax1.set_title("Late-exponential phase")
        ylimit = max(profile_rpf_0a_normalized)
        try :
            intiation = 450
            stop =  transcript_len[transcript] - 450
            axvspan(intiation,stop, color = "gold", alpha = .5)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop+1])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'gold', 'alpha': 0.5, 'pad': 5})
            intiation, stop = 289,451
            intiation, stop =388,451            

            axvspan(intiation,stop, color = "red", alpha = .3)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*10)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
        except :
            noncoding_reads = sum(profile_rpf_0a_normalized)
            noncoding_reads = str(int(noncoding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, noncoding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})





        panel1 += 10
        profile_rpf_0a    = [0 for n in range(cdsstart,cdsend )]
        profile_rpf_0a = read_bamfile(NormGlu1_rna_file4,"AP009493.1",profile_rpf_0a,riboseq, strand,cdsstart, cdsend)
        if sum(profile_rpf_0a) == 0: continue
        profile_rpf_0a_normalized = profile_rpf_0a
    
    
        fig.suptitle("S. griseus, adpA")
        
        transcript = name
        transcript_len = {}
        transcript_len[transcript] = len(profile_rpf_0a)

        profile_rpf_0a = np.array(profile_rpf_0a)/1
        transcript_seq =  str(sequence2)
        seq_line = list(transcript_seq) 
        xlimit = len(profile_rpf_0a)
        
        black_profile = [0 for n in range(len(profile_rpf_0a))]
        len_line = [n+1 for n in range(len(profile_rpf_0a))]
        
        ax1 = plt.subplot(gs[panel1:panel1+8,panel3:panel3+5])
        ax1.plot(profile_rpf_0a_normalized, color = "RoyalBlue", zorder=1)
        ax1.bar(np.arange(len(profile_rpf_0a_normalized)), profile_rpf_0a_normalized, width = 1, color = "RoyalBlue", ec = "RoyalBlue",zorder=2)
        ax1.plot(black_profile,color = "black", zorder=4)
        ax1.set_xlim(0,len(profile_rpf_0a))
        axesprep(ax1)
        #ax1.set_ylabel("Rep1_rpf",fontsize = 10,color='white',backgroundcolor="darkred")
        xtis = ax1.get_yticks()
        ax1.set_title("Stationary phase")
        ylimit = max(profile_rpf_0a_normalized)
        try :
            intiation = 450
            stop =  transcript_len[transcript] - 450
            axvspan(intiation,stop, color = "gold", alpha = .5)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop+1])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'gold', 'alpha': 0.5, 'pad': 5})
            
            intiation, stop = 289,451
            intiation, stop =388,451

            axvspan(intiation,stop, color = "red", alpha = .3)
            coding_reads = sum(profile_rpf_0a_normalized[intiation:stop])
            coding_reads = str(int(coding_reads))
            ax1.text((xlimit*10)*0.01, (ylimit*80)*0.01, coding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
        except :
            noncoding_reads = sum(profile_rpf_0a_normalized)
            noncoding_reads = str(int(noncoding_reads))
            ax1.text((xlimit*80)*0.01, (ylimit*80)*0.01, noncoding_reads, style='italic',fontsize = 12,bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})
        #ax1.set_ylim(0,400)

        start = [0 for n in range(len(transcript_seq))]
        start2 = [0 for n in range(len(transcript_seq))]
        start3 = [0 for n in range(len(transcript_seq))]
        end    = [0 for n in range(len(transcript_seq))]
        for match in re.finditer(r'(?=(%s))' % re.escape("ATG"), transcript_seq): start[match.start()] += 0.8
        for match in re.finditer(r'(?=(%s))' % re.escape("TTA"), transcript_seq): start2[match.start()] += 0.6
        for match in re.finditer(r'(?=(%s))' % re.escape("GTG"), transcript_seq): start3[match.start()] += 0.6
        for match in re.finditer(r'(?=(%s))' % re.escape("TAG" ),transcript_seq ): end[match.start()] += 1
        for match in re.finditer(r'(?=(%s))' % re.escape("TAA" ),transcript_seq ): end[match.start()] += 1
        for match in re.finditer(r'(?=(%s))' % re.escape("TGA" ),transcript_seq ): end[match.start()] += 1
        zero_frame_s, one_frame_s, two_frame_s = triple_frames(start)
        zero_frame_s2, one_frame_s2, two_frame_s2 = triple_frames(start2)
        zero_frame_s3, one_frame_s3, two_frame_s3 = triple_frames(start3)
        zero_frame_e, one_frame_e, two_frame_e = triple_frames(end)

        panel1 += 10
        #p#anel3 -=5 
        ax1 = plt.subplot(gs[panel1:panel1+2,panel3:panel3+5]) #to get a less high frame 1 panel could put panel1:panel1+1,panel3, but it might not be very legible
        
        index1 = 0
        intiation = 450
        stop =  transcript_len[transcript] - 450        
        axvspan(intiation,stop, color = "gold", alpha = .5)
                
        try :
            intiation = transcript_start[transcript] 
            stop =  transcript_end[transcript]
            if intiation%3 == 0: axvspan(intiation,stop, color = "gold", alpha = .5)
        except :pass
        ax1.spines['bottom'].set_visible(False)
        plt.plot(zero_frame_s, color = "green",alpha = 0.5)
        plt.plot(zero_frame_s2, color = "cyan",alpha = 0.5, lw = 3)
        plt.plot(zero_frame_e, color = "red",alpha = 0.5)
        plt.plot(zero_frame_s3, color = "brown",alpha = 0.5)
        ylim(0,1)
        xlim(0,len(zero_frame_s))
        xticks()
        xticks([],[])
        yticks([],[])

        panel2 += 2
        panel1 += 2
        ax1 = plt.subplot(gs[panel1:panel1+2,panel3:panel3+5]) #to get a less high frame 1 panel could put panel1:panel1+1,panel3, but it might not be very legible
        
        axvspan(388,451, color = "red", alpha = .3)        
        
        try :
            intiation = transcript_start[transcript] 
            stop =  transcript_end[transcript]
            if intiation%3 == 1: axvspan(intiation,stop, color = "gold", alpha = .5)
        except :pass
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        plt.plot(one_frame_s, color = "green",alpha = 0.5)
        plt.plot(one_frame_s2, color = "cyan",alpha = 0.5, lw = 3)
        plt.plot(one_frame_e, color = "red",alpha = 0.5)
        plt.plot(one_frame_s3, color = "brown",alpha = 0.5)
        ylim(0,1)
        xlim(0,len(zero_frame_s))
        xticks([],[])
        yticks([],[])
        panel2 += 2
        panel1 += 2
        
        
        
        ax1 = plt.subplot(gs[panel1:panel1+2,panel3:panel3+5]) #to get a less high frame 1 panel could put panel1:panel1+1,panel3, but it might not be very legible
        ax1.spines['top'].set_visible(False)
        try :
            intiation = transcript_start[transcript] 
            stop =  transcript_end[transcript]
            if intiation%3 == 2: axvspan(intiation,stop, color = "gold", alpha = .5)
        except :pass
        plt.plot(two_frame_s, color = "green",alpha = 0.5)
        plt.plot(two_frame_s2, color = "cyan",alpha = 0.5, lw = 3)
        plt.plot(two_frame_e, color = "red",alpha = 0.5)
        plt.plot(two_frame_s3, color = "brown",alpha = 0.5)
        plt.plot(1, color = "blue")
        ylim(0,1)
        xlim(0,len(zero_frame_s))
        yticks([],[])
        ax1.xaxis.set_major_locator(MaxNLocator(5))
        ax1.get_xaxis().tick_bottom()
        savefig(os.path.join(IMG_DIR, name))
        clf()
        panel1 = 0

