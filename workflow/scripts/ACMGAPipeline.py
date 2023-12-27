#!/home/cactus/cactus_env/bin/python


# coding=utf-8
from cactus.shared.configWrapper import ConfigWrapper
import xml.etree.ElementTree as ET
from cactus.progressive.progressive_decomposition import parse_seqfile
from cactus.progressive.progressive_decomposition import compute_outgroups
from cactus.progressive.progressive_decomposition import get_spanning_subtree
from cactus.progressive.progressive_decomposition import get_subtree
from sonLib.nxnewick import NXNewick
from sonLib.bioio import newickTreeParser
from cactus.paf.paf import get_event_pairs, get_leaves, get_node, get_distances
from cactus.paf.local_alignment import make_ingroup_to_outgroup_alignments_1
import re
import argparse


oldgenome_a=''
oldgenome_b=''
ancestorevent=''
event_tree_str=''
minimap2Paramters=''
proaliParamters=''
genoaliParamters=''
minimap2ForGff=''

fasta=[]
gff=[]
hal={}
treemap={}
sanitized = {}
outgroupsmap={}
input_seq_map={}
newgff = []

fastaVariety=[]
ancestor=[]
replacenamea=""
replacenameb=""
Cactus_threads=""
i=0
def make_paf_alignments(event_tree_string, event_names_to_sequences, ancestor_event_string, params):
    global event_tree_str,i,ancestorevent
    event_tree_str=event_tree_string
    i = i + 1
    print("mkdir " +PATH+"subtree"+str(i))
    event_tree = newickTreeParser(event_tree_string)
    ancestor_event = get_node(event_tree, ancestor_event_string)
    ancestorevent = ancestor_event.iD
    ingroup_events = get_leaves(ancestor_event)  # Get the set of ingroup events
    if ancestor_event.iD in event_names_to_sequences and event_names_to_sequences[ancestor_event.iD]:
        ingroup_events.append(ancestor_event)
    outgroup_events = [event for event in get_leaves(event_tree) if event not in ingroup_events]  # Set of outgroups

    for  ingroup, ingroup2, distance_a_b in get_event_pairs(ancestor_event, ingroup_events):
        print("#align " + event_names_to_sequences[ingroup.iD] + "\tagainst\t" + event_names_to_sequences[ingroup2.iD])
        print("\n" + "#m=============================================================================================================================================================================================================================" + str(i) + "\n")
        # print("line 60 ",event_names_to_sequences[ingroup.iD])
        if event_names_to_sequences[ingroup.iD] not in ancestor:
                for fasta_Variety in fastaVariety:
                    if fasta_Variety in event_names_to_sequences[ingroup.iD]:
                        genome_a_Variety=fasta_Variety
                        # print("line 61",genome_a_Variety)
        else:
            genome_a_Variety=event_names_to_sequences[ingroup.iD]
            # print("line 64",genome_a_Variety)
            gffn=genome_a_Variety+'.ancestor.gff'
            # print("line 69",gffn)
            if gffn not in newgff:
                print("minimap2 "+minimap2ForGff + PATH + genome_a_Variety + ".ref" + "   " +uniquePath+ " > " + PATH + genome_a_Variety + ".cds.sam")
                print("samtools view -O BAM " + PATH + genome_a_Variety + ".cds.sam" + " | " + "samtools sort - > " + PATH + genome_a_Variety + ".cds.bam")
                print("bedtools bamtobed -bed12 -i " + PATH + genome_a_Variety + ".cds.bam" + " | " + " bedtools sort > " + PATH + genome_a_Variety + ".cds.bed")
                print("bedToGenePred " + PATH + genome_a_Variety + ".cds.bed" + " " + PATH + genome_a_Variety + ".cds.genepred")
                print("genePredToGtf"+ " \"file\" " + PATH + genome_a_Variety + ".cds.genepred " + " " + PATH + genome_a_Variety + ".cds.gtf")
                print("gffread -E -M -K -Q -L -O " + PATH + genome_a_Variety + ".cds.gtf" + " -o- | sed 's/;locus=/;Parent=/g' | sed 's/geneID/gene_ID/g' | grep -v \"\slocus\s\">" + PATH + genome_a_Variety + ".gff")
                newgff.append(genome_a_Variety+ ".ancestor.gff")
        if event_names_to_sequences[ingroup2.iD] not in ancestor:
            for fasta_Variety in fastaVariety:
                if fasta_Variety in event_names_to_sequences[ingroup2.iD]:
                    genome_b_Variety = fasta_Variety
        else:
            genome_b_Variety=event_names_to_sequences[ingroup2.iD]
        # print("line 84",genome_a_Variety)
        gffnstr=[n for n in gff if genome_a_Variety in n][0]
        # print("line 87",gffnstr)
        if ((event_names_to_sequences[ingroup.iD] in fasta) and (event_names_to_sequences[ingroup2.iD] in ancestor)):
            print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + event_names_to_sequences[ ingroup.iD] + "  -o " +PATH+ genome_a_Variety + "_cds.fa")
            print("minimap2" +minimap2Paramters + PATH + event_names_to_sequences[ingroup2.iD] + ".ref " + PATH+genome_a_Variety + "_cds.fa > " +PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
            print("minimap2" +minimap2Paramters + PATH + event_names_to_sequences[ingroup.iD] + "  " + PATH+genome_a_Variety + "_cds.fa > " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
            print("/usr/bin/time anchorwave proali -i " + PATH + gffnstr + " -as " + PATH+genome_a_Variety + "_cds.fa  " + " -r " + PATH +event_names_to_sequences[ingroup.iD] + " -a " + PATH+genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH +event_names_to_sequences[ingroup2.iD] +".ref"+ " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + proaliParamters +" > " + genome_a_Variety + "_" + genome_b_Variety + ".log 2>&1")
            print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam")
            print("samtools view -H --reference " + PATH + event_names_to_sequences[ingroup.iD] + "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
            print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" +" >> " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        elif ((event_names_to_sequences[ingroup.iD] in fasta) and (event_names_to_sequences[ingroup2.iD] in fasta)):
            print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + event_names_to_sequences[ingroup.iD] + "  -o " + PATH+genome_a_Variety + "_cds.fa")
            print("minimap2 "+minimap2Paramters + PATH + event_names_to_sequences[ingroup2.iD]+ "  " +PATH+ genome_a_Variety + "_cds.fa > " +PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
            print("minimap2 "+minimap2Paramters + PATH + event_names_to_sequences[ingroup.iD] + "  " +PATH+ genome_a_Variety + "_cds.fa > " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
            print("/usr/bin/time anchorwave genoAli -i " + PATH + gffnstr + " -as " +PATH+ genome_a_Variety + "_cds.fa  " + " -r " + PATH +event_names_to_sequences[ingroup.iD] + " -a " +PATH+  genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH +event_names_to_sequences[ingroup2.iD] + " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + genoaliParamters + " > " + genome_a_Variety + "_" + genome_b_Variety + ".log 2>&1")
            print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam")
            print("samtools view -H --reference " + PATH + event_names_to_sequences[ingroup.iD] + "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
            print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" +" >> " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        elif((event_names_to_sequences[ingroup.iD] not in fasta) and (event_names_to_sequences[ingroup2.iD] in fasta)):
            print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + event_names_to_sequences[ingroup.iD] +".ref"+ "  -o " + PATH+genome_a_Variety + "_cds.fa")
            print("minimap2  "+minimap2Paramters + PATH + event_names_to_sequences[ingroup2.iD] +" " +        PATH+genome_a_Variety + "_cds.fa > " +PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
            print("minimap2  "+minimap2Paramters + PATH + event_names_to_sequences[ingroup.iD] +".ref"+ " " + PATH+genome_a_Variety + "_cds.fa > " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
            print("/usr/bin/time anchorwave proali -i " + PATH + gffnstr + " -as " + PATH+genome_a_Variety + "_cds.fa  " + " -r " + PATH +event_names_to_sequences[ingroup.iD] + ".ref" + " -a " +PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH +
            event_names_to_sequences[ingroup2.iD] + " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + proaliParamters+"  > " + genome_a_Variety + "_" + genome_b_Variety + ".log 2>&1")
            print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam")
            print("samtools view -H --reference " + PATH + event_names_to_sequences[ingroup.iD]+".ref" + "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
            print("grep -v \"^@\" " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" + " >> " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        else:
            print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + event_names_to_sequences[ingroup.iD] +".ref"+ "  -o " + PATH+genome_a_Variety + "_cds.fa")
            print("minimap2  "+minimap2Paramters + PATH + event_names_to_sequences[ingroup2.iD]+".ref"+ " " + PATH+genome_a_Variety + "_cds.fa > " + PATH+genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
            print("minimap2  "+minimap2Paramters + PATH + event_names_to_sequences[ingroup.iD] +".ref"+ " " + PATH+genome_a_Variety + "_cds.fa > " + PATH+genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
            print("/usr/bin/time anchorwave proali -i "+PATH+gffnstr+" -as "+PATH+genome_a_Variety+"_cds.fa  "+" -r "+PATH+event_names_to_sequences[ingroup.iD]+".ref"+" -a "+PATH+genome_b_Variety+"_"+genome_a_Variety+"_cds.sam "+" -ar "+PATH+genome_a_Variety+"_"+genome_a_Variety+"_ref.sam"+" -s "+PATH+event_names_to_sequences[ingroup2.iD]+".ref"+" -n "+PATH + genome_a_Variety + "_" + genome_b_Variety +".anchors" \
                  + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + proaliParamters +" > "+genome_a_Variety+"_"+genome_b_Variety+".log 2>&1")
            print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam")
            print("samtools view -H --reference " + PATH + event_names_to_sequences[ingroup.iD] +".ref"+ "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
            print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.sam" +" >> " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")

        print("k8 /data/workflow/envs/paftools.js" + " sam2paf "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".2"+".anchorwave.sam"+" > "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".tmp.anchorwave.paf")
        for (genome, seq) in input_seq_map.items():
            global replacenamea
            global replacenameb
            if (seq in fasta):
                if genome_a_Variety in genome:
                    replacenamea=genome
                if genome_b_Variety in genome:
                    replacenameb=genome
        if genome_a_Variety in ancestor:
            replacenamea = genome_a_Variety
        if genome_b_Variety in ancestor:
            replacenameb = genome_b_Variety
        if replacenamea and replacenameb:
            print("python /data/workflow/scripts/replace_ref_que.py -i " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".tmp.anchorwave.paf" + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.paf" + " -r " + replacenamea+" -q " +replacenameb)
                # print("python replace_ref_que.py -i " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".tmp.anchorwave.paf" + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.paf" + " -r " + genome_a_Variety +"_protein"+ " -q " + genome_b_Variety+"_protein")
        print("cat "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave.paf"+" > "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave_invert.paf")
        print("paf_invert -i "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave.paf"+" >> "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave_invert.paf")
        print("paf_chain -i "+ PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave_invert.paf"+" --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG "+ " >> "+PATH+'subtree'+str(i)+"/anchorwave_chain.paf")

    distances = get_distances(event_tree)  # Distances between all pairs of nodes
    outgroup_events.sort(key=lambda outgroup: distances[ancestor_event, outgroup])  # Sort from closest to furthest
    if int(params.find("blast").attrib["trimIngroups"]):  # Trim the ingroup sequences
        for  ingroup in ingroup_events:
            if len(outgroup_events) > 0:
                global oldgenome_a
                global oldgenome_b
                oldgenome_a =event_names_to_sequences[ingroup.iD]
                oldgenome_b=event_names_to_sequences[outgroup_events[0] .iD]
                make_ingroup_to_outgroup_alignments_1(ingroup, outgroup_events,dict(event_names_to_sequences), distances, params)
    else:
        for ingroup in ingroup_events:
            for outgroup in outgroup_events:
                make_chunked_alignments(ingroup.iD, event_names_to_sequences[ingroup.iD],outgroup.iD, event_names_to_sequences[outgroup.iD], distances[ingroup, outgroup], params)
                print("m from make_paf_alignments else")

def make_ingroup_to_outgroup_alignments_1(ingroup_event, outgroup_events, event_names_to_sequences,distances, params):
    outgroup = outgroup_events[0]  # The first outgroup
    alignment=make_chunked_alignments(outgroup.iD, event_names_to_sequences[outgroup.iD],ingroup_event.iD, event_names_to_sequences[ingroup_event.iD], distances[ingroup_event, outgroup], params)
    if len(outgroup_events) > 1:
        make_ingroup_to_outgroup_alignments_2(alignment, ingroup_event, outgroup_events[1:],event_names_to_sequences, distances, params)
    else:
        if event_names_to_sequences[outgroup.iD] in fasta:
            for fasta_Variety in fastaVariety:
                if fasta_Variety in event_names_to_sequences[outgroup.iD]:
                    genome_a_Variety = fasta_Variety
        else:
            genome_a_Variety=event_names_to_sequences[outgroup.iD]
        if event_names_to_sequences[ingroup_event.iD] in fasta:
            for fasta_Variety in fastaVariety:
                if fasta_Variety in event_names_to_sequences[ingroup_event.iD]:
                    genome_b_Variety = fasta_Variety
        else:
            genome_b_Variety=event_names_to_sequences[ingroup_event.iD]

        if oldgenome_b in fasta:
            for fasta_Variety in fastaVariety:
                if fasta_Variety in oldgenome_b:
                    oldgenome_a_Variety = fasta_Variety
        else:
            oldgenome_a_Variety = oldgenome_b

        if oldgenome_a in fasta:
            for fasta_Variety in fastaVariety:
                if fasta_Variety in oldgenome_a:
                    oldgenome_b_Variety=fasta_Variety
        else:
            oldgenome_b_Variety = oldgenome_a
        if((genome_a_Variety==oldgenome_a_Variety)and(genome_b_Variety==oldgenome_b_Variety)):
            print("cat " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.paf" +  " > " + PATH + 'cat_' + oldgenome_a_Variety + '_' + oldgenome_b_Variety + '_' + genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave.paf")
        else:
            print("cat " +  PATH + genome_a_Variety + "_" + genome_b_Variety +".anchorwave.paf"+"     "+PATH+oldgenome_a_Variety+"_"+oldgenome_b_Variety+".anchorwave.paf" + " > " +  PATH + 'cat_'+oldgenome_a_Variety+'_'+oldgenome_b_Variety+'_'+genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave.paf")
        print("cat " +PATH + 'cat_'+oldgenome_a_Variety+'_'+oldgenome_b_Variety+'_'+genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave.paf" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave_invert.paf")
        print("paf_invert -i " + PATH + 'cat_'+oldgenome_a_Variety+'_'+oldgenome_b_Variety+'_'+genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave.paf" + " >> " +PATH + genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave_invert.paf")
        print("paf_chain -i " +PATH + genome_a_Variety + "_" + genome_b_Variety + ".2.anchorwave_invert.paf" + " --maxGapLength 1000000 --chainGapOpen 5000 --chainGapExtend 1 --trimFraction 0.02 --logLevel DEBUG " + " >> "+PATH+'subtree'+str(i)+"/anchorwave_chain.paf")

def make_chunked_alignments( event_a, genome_a, event_b, genome_b, distance, params):
    print("\n" + "#m============================================================================================================================================================================================================================="+str(i) + "\n")
    print("#genome_a   "+genome_a+"   genome_b   "+genome_b)
    # print("line 184",genome_a,fasta)
    if (genome_a in fasta):
        for fasta_Variety in fastaVariety:
            if fasta_Variety in genome_a:
                genome_a_Variety = fasta_Variety
                # print("line 198",genome_a_Variety)
    else:
        genome_a_Variety = genome_a
        # print("line 192", genome_a_Variety)
        gffn = genome_a_Variety + '.ancestor.gff'
        if gffn not in newgff:
            print("minimap2 "+minimap2ForGff + PATH + genome_a_Variety + ".ref" + "   " + uniquePath+" > " + PATH + genome_a_Variety + ".cds.sam")
            print("samtools view -O BAM " + PATH + genome_a_Variety + ".cds.sam" + " | " + "samtools sort - > " + PATH + genome_a_Variety + ".cds.bam")
            print( "bedtools bamtobed -bed12 -i " + PATH + genome_a_Variety + ".cds.bam" + " | " + " bedtools sort > " + PATH + genome_a_Variety + ".cds.bed")
            print("bedToGenePred " + PATH + genome_a_Variety + ".cds.bed" + " " + PATH + genome_a_Variety + ".cds.genepred")
            print("genePredToGtf"+ " \"file\" " + PATH + genome_a_Variety + ".cds.genepred " + " " + PATH + genome_a_Variety + ".cds.gtf")
            print("gffread -E -M -K -Q -L -O " + PATH + genome_a_Variety + ".cds.gtf" + " -o- | sed 's/;locus=/;Parent=/g' | sed 's/geneID/gene_ID/g' | grep -v \"\slocus\s\">" + PATH + genome_a_Variety + ".gff")
            newgff.append(genome_a_Variety + ".ancestor.gff")
    if genome_b in fasta:
       for fasta_Variety in fastaVariety:
           if fasta_Variety in genome_b:
               genome_b_Variety = fasta_Variety
    else:
       genome_b_Variety=genome_b
    # print('line 204',genome_a_Variety)
    gffname = [n for n in gff if genome_a_Variety in n]
    # print("line 206",gffname)
    gffnstr = gffname[0]
    if ((genome_a in fasta )and(genome_b not in fasta)):
        print("anchorwave gff2seq -i " + PATH + gffnstr+ "   -r " + PATH +genome_a + "  -o " + PATH+ genome_a_Variety+"_cds.fa")
        print("minimap2  " +minimap2Paramters+ PATH + genome_b + ".ref " + genome_a_Variety + "_cds.fa > " + PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
        print("minimap2  " +minimap2Paramters+ PATH + genome_a + "  " + genome_a_Variety + "_cds.fa > "     +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
        print("/usr/bin/time anchorwave proali -i " + PATH + gffnstr + " -as " + PATH+ genome_a_Variety + "_cds.fa  " + " -r " + PATH + genome_a + " -a " +PATH+  genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " + PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH + genome_b+".ref" + " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + proaliParamters+" > " + genome_a_Variety + "_" + genome_b_Variety + ".log 2>&1")
        print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + "  > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        print("samtools view -H --reference " + PATH + genome_a + "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")
        print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam"+" >> "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")

    elif ((genome_a in fasta) and (genome_b in fasta)):
        print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + genome_a  + "  -o " +PATH+  genome_a_Variety + "_cds.fa")
        print("minimap2 "+minimap2Paramters + PATH + genome_b  + "  " + PATH+ genome_a_Variety + "_cds.fa > " +PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
        print("minimap2 "+minimap2Paramters + PATH + genome_a  + "  " + PATH+ genome_a_Variety + "_cds.fa > " +PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
        print("/usr/bin/time anchorwave genoAli -i " + PATH + gffnstr + " -as " + PATH+ genome_a_Variety + "_cds.fa  " + " -r " + PATH + genome_a  + " -a " + PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " +PATH+  genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH + genome_b + " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + genoaliParamters+" > " + genome_a_Variety + "_" + genome_b_Variety + ".avx2.log 2>&1")
        print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + "  > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        print("samtools view -H --reference " + PATH + genome_a + "  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")
        print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam"+" >> "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")

    elif((genome_a not in fasta)and(genome_b  in fasta)):
        print("anchorwave gff2seq -i " + PATH + gffnstr+ "   -r " + PATH +genome_a+".ref" + "  -o " +PATH+  genome_a_Variety+"_cds.fa")
        print("minimap2 "+minimap2Paramters + PATH + genome_b +"  "+          genome_a_Variety + "_cds.fa > " + PATH+ genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
        print("minimap2 "+minimap2Paramters + PATH + genome_a +".ref"+ "  " + genome_a_Variety + "_cds.fa > " + PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
        print("/usr/bin/time anchorwave proali -i " + PATH + gffnstr + " -as " + PATH+ genome_a_Variety + "_cds.fa  " + " -r " + PATH + genome_a + ".ref" + " -a " +PATH+  genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " + " -ar " +PATH+  genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH + genome_b + " -n " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchors" \
            + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + " -f " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.f.maf" + proaliParamters+"> " + genome_a_Variety + "_" + genome_b_Variety + ".avx2.log 2>&1")
        print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + "  > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        print("samtools view -H --reference " + PATH + genome_a + ".ref  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")
        print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam"+" >> "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")

    else:
        print("anchorwave gff2seq -i " + PATH + gffnstr + "   -r " + PATH + genome_a + ".ref" + "  -o " + PATH+ genome_a_Variety + "_cds.fa")
        print("minimap2 "+minimap2Paramters + PATH + genome_b + ".ref " +       genome_a_Variety + "_cds.fa > " +PATH+  genome_b_Variety + "_" + genome_a_Variety + "_cds.sam")
        print("minimap2 "+minimap2Paramters + PATH + genome_a + ".ref" + "  " + genome_a_Variety + "_cds.fa > " +PATH+  genome_a_Variety + "_" + genome_a_Variety + "_ref.sam")
        print("/usr/bin/time anchorwave proali -i " + PATH + gffnstr  + " -as " +PATH+  genome_a_Variety + "_cds.fa  " + " -r " + PATH + genome_a+".ref" + " -a " +PATH+  genome_b_Variety + "_" + genome_a_Variety + "_cds.sam " +" -ar " + PATH+ genome_a_Variety + "_" + genome_a_Variety + "_ref.sam" + " -s " + PATH + genome_b +".ref"+ " -n " + PATH +genome_a_Variety+"_"+genome_b_Variety + ".anchors" \
              + " -o " + PATH +genome_a_Variety+"_"+genome_b_Variety +".anchorwave.maf" + " -f " + PATH +genome_a_Variety+"_"+genome_b_Variety+ ".anchorwave.f.maf"+proaliParamters+" > " +genome_a_Variety+"_"+genome_b_Variety+".avx2.log 2>&1")
        print("maf-convert sam " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.maf" + "  > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam")
        print("samtools view -H --reference " + PATH + genome_a + ".ref  " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam" + " > " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")
        print("grep -v \"^@\" "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".2" + ".anchorwave.sam"+" >> "+ PATH + genome_a_Variety + "_" + genome_b_Variety + ".3" + ".anchorwave.sam")
    print("k8 /data/workflow/envs/paftools.js  sam2paf " + PATH +genome_a_Variety+"_"+genome_b_Variety +".3"+".anchorwave.sam" + " > "+ PATH +genome_a_Variety+"_"+genome_b_Variety +".tmp.anchorwave.paf")
    for (genome, seq) in input_seq_map.items():
            global replacenamea
            global replacenameb
            if (seq in fasta):
                if genome_a_Variety in genome:
                    replacenamea=genome
                if genome_b_Variety in genome:
                    replacenameb=genome
    if genome_a_Variety in ancestor:
        replacenamea = genome_a_Variety
    if genome_b_Variety in ancestor:
        replacenameb = genome_b_Variety
    if genome_a_Variety and genome_b_Variety:
            print("python /data/workflow/scripts/replace_ref_que.py -i " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".tmp.anchorwave.paf" + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.paf" + " -r " + replacenamea+" -q " +replacenameb)
    # print( "python replace_ref_que.py -i " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".tmp.anchorwave.paf" + " -o " + PATH + genome_a_Variety + "_" + genome_b_Variety + ".anchorwave.paf" + " -r " + genome_a_Variety + "_protein" + " -q " + genome_b_Variety + "_protein")


def make_ingroup_to_outgroup_alignments_2( alignments, ingroup_event, outgroup_events,event_names_to_sequences, distances, params):
    make_ingroup_to_outgroup_alignments_1(ingroup_event, outgroup_events,event_names_to_sequences, distances, params)

def cactus_cons_with_resources(spanning_tree, event, config_node, subtree_eventmap,og_map):
    resultPath = "result/"
    tmpdict={}
    for event, seq_id in subtree_eventmap.items():
        if(event in sanitized):
            tmpdict[event]=sanitized[event]
    outgroups = og_map[ancestorevent] if ancestorevent in og_map else []
    pairs = [[genome, faPath] for genome, faPath in list(tmpdict.items())]
    print("paf_tile -i "   +PATH+'subtree'+str(i)+"/anchorwave_chain.paf" + " --logLevel DEBUG " + " > " +  PATH+'subtree'+str(i)+"/anchorwave_tile.paf")
    print("grep -v tp:A:S "+PATH+'subtree'+str(i)+"/anchorwave_tile.paf"+ " > "+ PATH + str(ancestorevent) + "_primary.paf || true")
    print("grep  tp:A:S "  +PATH+'subtree'+str(i)+"/anchorwave_tile.paf"+ "> " + PATH  +str(ancestorevent) + "_secondary.paf || true")
    if(outgroups==[]):
        print("cactus_consolidated --sequences " +"\" "+" ".join([item for tmpdict in pairs for item in tmpdict])+ " \" "+" --speciesTree " +" \" "+ event_tree_str+ "\""+" --logLevel DEBUG --alignments  " + PATH + str(ancestorevent) + "_primary.paf "+ "  --params "+ xmlPath+ " --outputFile " +PATH+ ancestorevent + ".c2h " + " --outputHalFastaFile " +PATH+ ancestorevent + ".c2h.fa"+ " --outputReferenceFile " +PATH+ ancestorevent+ ".ref"+" ".join(outgroups)+ " --referenceEvent " + ancestorevent+" "+Cactus_threads+ " --secondaryAlignments " +PATH + str(ancestorevent) + "_secondary.paf")
    else:
        print("cactus_consolidated --sequences " +"\" "+" ".join([item for tmpdict in pairs for item in tmpdict])+ " \" "+" --speciesTree " +" \" "+ event_tree_str+ "\""+" --logLevel DEBUG --alignments  " + PATH + str(ancestorevent) + "_primary.paf "+ "  --params "+ xmlPath+ " --outputFile " +PATH+ ancestorevent + ".c2h " + " --outputHalFastaFile " +PATH+ ancestorevent + ".c2h.fa"+ " --outputReferenceFile " +PATH+ ancestorevent+ ".ref" + " --outgroupEvents " +" ".join(outgroups)+ " --referenceEvent " + ancestorevent+ " "+Cactus_threads+ " --secondaryAlignments " +PATH + str(ancestorevent) + "_secondary.paf")
    pattern=f"\(([^(]*){ancestorevent}"
    haltree=re.search(pattern,event_tree_str)
    if(outgroups==[]):
        hal[ancestorevent]="halAppendCactusSubtree "+PATH+ancestorevent+".c2h"+"  "+PATH+ancestorevent+".c2h.fa"+" "+" \""+haltree.group(0)+ ";"+"\" "+"  "+" "+resultPath+"evolverPlants.hal"
    else:
        hal[ancestorevent]="halAppendCactusSubtree "+PATH+ancestorevent+".c2h"+"  "+PATH+ancestorevent+".c2h.fa"+" "+" \""+haltree.group(0)+ ";"+"\" "+"  "+" "+resultPath+"evolverPlants.hal"+" --outgroups "+",".join(outgroups)
    treemap[ancestorevent]=event_tree_str
    outgroupsmap[ancestorevent]=' '.join(outgroups)

def progressive_step(config_node, seq_id_map, tree, og_map, event):

    subtree = get_subtree(tree, event, ConfigWrapper(config_node), og_map)

    # get our events
    subtree_eventmap = {}
    for leaf in subtree.getLeaves():
        subtree_eventmap[subtree.getName(leaf)] = seq_id_map[subtree.getName(leaf)]

    # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
    if event in seq_id_map and seq_id_map[event]:
        subtree_eventmap[event] = seq_id_map[event]

    # get the spanning tree (which is what consolidated wants)
    spanning_tree = get_spanning_subtree(tree, event, ConfigWrapper(config_node), og_map)

    make_paf_alignments(NXNewick().writeString(spanning_tree), subtree_eventmap, event, config_node)
    if int(config_node.find("blast").attrib["trimOutgroups"]):  # Trim the outgroup sequences
        outgroups = og_map[event] if event in og_map else []
        progressive_step_2(config_node, subtree_eventmap,spanning_tree, og_map, event)

    else:  # Without outgroup trimming
        cactus_cons_with_resources(spanning_tree, event, config_node, subtree_eventmap,og_map)

def halAppendCactusSubtree():
    items=list(hal.items())
    for key,value in  reversed(items):
        print(value)


def progressive_step_2(config_node, subtree_eventmap,spanning_tree, og_map, event):
        cactus_cons_with_resources(spanning_tree, event, config_node, subtree_eventmap,og_map)

def main(xml,evolverMammals,snakemakePath,snakemakeFasta,snakemakeGff,uniqueCds,sVariety,sminimap2ForGff,sminimap2Paramters,sproaliParamters,sgenoaliParamters,scactus_threads):
    global xmlPath
    global PATH
    global fasta
    global gff
    global uniquePath
    global fastaVariety
    global ancestor
    global minimap2ForGff
    global minimap2Paramters
    global proaliParamters
    global genoaliParamters
    global input_seq_map
    global Cactus_threads
    xmlPath="/data/"+xml
    PATH="/data/"+snakemakePath
    uniquePath="/data/"+uniqueCds
    fasta=snakemakeFasta.split(',')
    gff=snakemakeGff.split(',')
    fastaVariety=sVariety.split(',')

    minimap2ForGff=sminimap2ForGff
    minimap2Paramters=sminimap2Paramters
    proaliParamters=sproaliParamters
    genoaliParamters=sgenoaliParamters
    Cactus_threads=scactus_threads
    config_node = ET.parse(xml).getroot()
    config_wrapper = ConfigWrapper(config_node)
    config_wrapper.substituteAllPredefinedConstantsWithLiterals()

    mc_tree, input_seq_map, og_candidates = parse_seqfile(evolverMammals, config_wrapper)
    og_map = compute_outgroups(mc_tree, config_wrapper, set(og_candidates), None)


    root_event = mc_tree.getRootName()
    dep_table = {}
    subtree_event_set = set(mc_tree.getSubtreeRootNames())
    ancestor=list(subtree_event_set)
    for i in range(1000):
        ancestor.append("Anc"+str(i))
        gff.append("Anc"+str(i)+".gff")
    # print("line 394",ancestor,"\n",gff)
    iteration_0_events = []
    for (genome, seq) in input_seq_map.items():
        if(seq in fasta):
            print("cactus_sanitizeFastaHeaders  "+ PATH+seq+"  "+genome +" > " +PATH+genome+".sanitized.fa")

    for node in mc_tree.postOrderTraversal(mc_tree.getNodeId(root_event)):
        event = mc_tree.getName(node)
        global sanitized
        if event not in subtree_event_set:
            sanitizedPATH=PATH+event+'.sanitized.fa'
            sanitized[event]=sanitizedPATH
        if event in subtree_event_set:
            sanitizedPATH=PATH+event+'.ref'
            sanitized[event]=sanitizedPATH
            dep_table[event] = []
            subtree = get_subtree(mc_tree, event, config_wrapper, og_map)
            input_count = 0
            for leaf_id in subtree.getLeaves():
                leaf_name = subtree.getName(leaf_id)
                dep_table[event].append(leaf_name)
                if leaf_name in input_seq_map:
                    input_count += 1
            if input_count == len(dep_table[event]):
                iteration_0_events.append(event)
    assert len(dep_table) > 0

    # make the jobs.  jobs must be made after their dependencies.  we just brute force it
    job_table = {}
    fasta_results = input_seq_map


    while len(job_table) != len(dep_table):
        events = [k for k in dep_table.keys()]
        for event in events:
            if event in job_table:
                continue
            unscheduled_deps = False
            for dep in dep_table[event]:
                if dep not in fasta_results:
                    unscheduled_deps = True
                    break
            if not unscheduled_deps:
                # found an event that a) hasn't been made into a job yet an b) has an existing job for every dep
                event_id_map = {}
                for dep in dep_table[event]:
                    event_id_map[dep] = fasta_results[dep]
                # to be consistent with pre-refactor (and work with updating tests), we include the root when its id is input
                if event in input_seq_map and input_seq_map[event]:
                    event_id_map[event] = input_seq_map[event]

                progressive_step(config_node, event_id_map, mc_tree, og_map, event)
                job_table[event] = "event_job"
                fasta_results[event] = event
    halAppendCactusSubtree()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ACMGA required parameters')
    parser.add_argument('--xml',type=str,help='Cactus runs the required xml')
    parser.add_argument('--tree',type=str,help='A guide tree tuns a progressive algorithm')
    parser.add_argument('--path',type=str,help='Path of FASTA and GFF files, as well as non-duplicate CDS data')
    parser.add_argument('--fasta',help='Collection of input FASTA file names')
    parser.add_argument('--gff',help='Collection of input gff and ancestor gff file names')
    parser.add_argument('--nonDuplicateCDS',type=str,help='Collection of non-duplicate CDS for input genomes')
    parser.add_argument('--species',type=str,help='Species names corresponding to FASTA and GFF files')
    parser.add_argument('--minimap2ForGff', type=str,help='Parameter settings for minimap2 in ACMGA to liftover ancestor.fa')
    parser.add_argument('--minimap2Paramters',type=str,help='Parameter settings for minimap2 in ACMGA using AnchorWave')
    parser.add_argument('--proaliParamters',type=str,help='Parameter settings for proali using AnchorWave in ACMGA')
    parser.add_argument('--genoaliParamters',type=str,help='Parameter settings for genoali using AnchorWave in ACMGA')
    parser.add_argument('--cactus_threads',type=str,help='Number of threads used in the cactus')
    args = parser.parse_args()

    main(args.xml,args.tree,args.path,args.fasta,args.gff,args.nonDuplicateCDS,args.species,args.minimap2ForGff,args.minimap2Paramters,args.proaliParamters,args.genoaliParamters,args.cactus_threads)
