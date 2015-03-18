#!/usr/bin/python

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This script is not part of this project... it's accessory code from before that i am using as a template
# for a new script.  it will be removed later, so if you are reading this piece of code, it not being
# used in the project.  Sorry to waste your time.

# I do not know how many of these i need, i will edit later
import os
import sys
import time
import argparse
import shutil
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Application
from Bio.Application import _Option
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
import subprocess

# This folder and its contents will be removed after this program concludes its run. i just need a place to store a bunch of intermediate files.
tmp_directory = "./msa_tmp/"
try:
    os.mkdir(tmp_directory)
except:
    pass


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is to build a newick format phylogenetic tree from a list of genomes and a marker gene.')
    
    parser.add_argument("-i", "--infolder", dest="infolder", metavar="DIRECTORY", default='./genomes/',
                help="Directory containing all genbank files for use by the program.")
    
    parser.add_argument("-o", "--outfile", dest="outfolder", metavar="FILE", default='./out_tree.nwk',
                help="Newick format tree containing the result of the phylogenetic tree built with a marker gene.")
                
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='./phylo_order.txt',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    if parsed_args.filter == 'NONE' or os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
    
    outfile = parsed_args.outfolder
    
    return infolder, outfile, filter_file
    
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

def return_file_list(infolder, filter_file):
    if filter_file == '':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]
        

def create_distmat(fname, method = 1):
    cline = ClustalwCommandline('clustalw', infile=fname)
    print "cline", cline
    #base = fname.split('/')[len(fname.split('/')) - 1].split('.')[0]
    #better = fname.split('.')[0]
    #print 'better: ', better
    #base = fname.replace('.', '_tst.aln')

    return_code = subprocess.call(str(cline),stdout = open(os.devnull), stderr = open(os.devnull),shell=(sys.platform!="win32"))
    #print "return_code", return_code
    
    
    distmat_infile = '.'.join(fname.split('.')[:2]) + '.aln'
    #distmat_outfile = '.'.join(fname.split('.')[:2]) + '.distmat'
    distmat_outfile = '.'.join(fname.split('.')[:2])
    print "fname", fname, "distmat_infile", distmat_infile, "distmat_outfile", distmat_outfile
    #distmat_line = "distmat %s.aln -outfile %s.distmat  -protmethod %i" % (better, better, method) 
    distmat_line = "distmat %s -outfile %s.distmat -protmethod %i" % (distmat_infile, distmat_outfile, method) 
    print 'distmat_line: ', distmat_line
    return_code = subprocess.call(distmat_line, stdout = open(os.devnull), stderr = open(os.devnull), shell=(sys.platform!="win32"))
    print "return_code", return_code
    
    
    # make the newick tree
    cline = ClustalwCommandline('clustalw', infile=fname, tree=1)
    print "cline", cline
    return_code = subprocess.call(str(cline),stdout = open(os.devnull), stderr = open(os.devnull),shell=(sys.platform!="win32"))
    


        
# this function will take a marker gene and a list of genbank files
# and construct a fasta file for their sequences.  For the time being at least
# this function will only do protein marker genes. Later i will expand this function
# to take RNA genes like 16s as well.
# TODO : expand the functionality of this function to make use of RNA genes and give the user the ability to list more than one marker that they are interested in
def make_target_fasta(marker, infolder, filter_file):
    org_paths = return_file_list(infolder, filter_file)
    result = [] # this is not great
    common_to_accession_dict = {}
    orgs_with_marker = []
    marker = marker.lower()
    for org in org_paths:
        genes_found = []
        seq_record = SeqIO.parse(open(org), "genbank").next()
        accession = seq_record.annotations['accessions'][0]
        organism_tmp = seq_record.annotations['organism'].replace(' ', '_')
        # Put code here to determine the format of the organisms' english name. currently i am using genus species, but strain can also be used
        organism = '_'.join(organism_tmp.split('_')[:2])
        
        # Here we store the {organism:accession} information so that we can build a new list that is needed for the visualization pipeline
        common_to_accession_dict.update({organism:accession})
        
        found = False
        for fnum, feature in enumerate(seq_record.features):            
            if feature.type == 'CDS':
                #start = feature.location._start.position
                start = feature.location.start
                stop = feature.location.end
                try: 
                    gene = feature.qualifiers['gene'][0]
                    gene = gene.lower()
                except:
                    gene = 'unknown'
                if gene == marker:
                    genes_found.append(gene)
                    seq = feature.qualifiers['translation'] # this returns the protein product, not suitable for RNA products like 16s
                    #result.append(">%s|%s|%s" % (accession, organism, gene))
                    result.append(">%s" % organism)
                    result.append(''.join(seq))
                    orgs_with_marker.append(accession)
                    break
            
    outfile = tmp_directory + "distmat_marker.fa"
    print "outfile", outfile
    handle = open(outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()
    print "orgs_with_marker", len(orgs_with_marker)
    create_distmat(outfile)
    
    print "common_to_accession_dict", common_to_accession_dict.keys()
    return "./msa_tmp/distmat_marker.ph", common_to_accession_dict


def return_tree_order_list(newick_tree_file, common_to_accession_dict):
    result = []
    tree = Phylo.read("out_tree.nwk", "newick")
    for clade in tree.find_clades():
        if clade.name != None:
            #accession, common_tmp, marker = clade.name.split('|')
            #common = '_'.join(common_tmp.split('_')[:2])
            common = clade.name
            accession = common_to_accession_dict[common]
            result.append(','.join([accession, common]))
            clade = common
    print '\n'.join(result)
    
    Phylo.write(tree, "./test_tree.nwk", "newick")


def main():
    
    start = time.time()    

    parsed_args = parser_code()
    
    infolder, outfile, filter_file = check_options(parsed_args)
    
    marker_gene = "rpob"
    
    print infolder, outfile, filter_file
    
    tree_file, common_to_accession_dict = make_target_fasta(marker_gene, infolder, filter_file)
    shutil.copy(tree_file, outfile)
    
    return_tree_order_list(tree_file, common_to_accession_dict)
    
    # get rid of the temp files
    shutil.rmtree(tmp_directory)
    
    print time.time() - start
    
if __name__ == '__main__':
    main()    

