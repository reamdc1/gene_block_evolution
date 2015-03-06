#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from homolog4 import *
from collections import defaultdict

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="Parse the results of a BLAST -m8 search and organize the results by specific operons. The program will save the results in a folder designated by the user or the default './blast_parse/'.")
                
    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_result/', metavar="DIRECTORY",
                help="A file that contains the path to every organism database that you are interested in.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./blast_parse/',
                help="Folder where the BLAST results will be stored. Default is the folder './blast_result/'.")
    
    parser.add_argument("-q", "--operon_query", dest="operon_query", default='./regulonDB/operon_names_and_genes.txt', metavar="FILE",
                help="A file that contains the names and genes comprising the operons that are under investigation.")

    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the accession numbers of the organisms that are under investigation.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
        
    if os.path.exists(parsed_args.operon_query):
        operon_query = parsed_args.operon_query
    else:
        print "The file %s does not exist." % parsed_args.operon_query
        sys.exit()
    
    
    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
        

    return infolder, outfolder, operon_query, filter_file, num_proc


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# The filter file here i think should be either a vacant value (such as '') or a user defined
# value.  I do not think by default it should be given, like I have made the default behavior.    
def parallel_blast_parse_dict(in_folder, out_folder, num_proc, filter_file, operon_dict):
    result = {}
    operon_out_folder = out_folder
    if not os.path.isdir(operon_out_folder):
        os.makedirs(operon_out_folder)
    if filter_file != '':
        tmp = returnRecursiveDirFiles(in_folder)
        nc_list = [i.strip() for i in open(filter_file).readlines()]
        fname_list = [i for i in tmp if i.split('/')[-1].split('.')[0] in nc_list]
    else:
        fname_list = returnRecursiveDirFiles(in_folder)
    
    #print fname_list
    for fname in fname_list:
        for line in [i.strip() for i in open(fname).readlines()]:
            try:
                hlog = Homolog.from_blast(line)
            except:
                print "ERROR", line
            #hlog.Print()
            
            # this might have to be changed.... 
            #accession = hlog.accession()
            try:
                accession = hlog.accession()
            except:
                print line
            
            
            predicted_gene = hlog.blast_annatation()
            
            '''
            try: # faster implementation than "if predicted_gene in operon_dict.keys():"
                operon = operon_dict[predicted_gene]
                # Debugging the missing operon casABCDE12... no idea right now.
                if operon == 'casABCDE12':
                    print 'AFDFAFDSF'
                if operon in result.keys():    # check if the operon is in the result dict already, if not make a new entry in the else clause
                    if accession in result[operon].keys(): # Check if the organism has been added to the operon
                        result[operon][accession].append(hlog.ret_str())
                    else: # if the organims is not part of the operon, add it
                        result[operon].update({accession:[hlog.ret_str()]})
                else: # add the operon to the result
                    result.update({operon: {accession: [hlog.ret_str()]}})
            except:
                pass
                
            '''
            try: # faster implementation than "if predicted_gene in operon_dict.keys():"
                operon = operon_dict[predicted_gene]
                # Debugging the missing operon casABCDE12... no idea right now.
                # ok, it's there, so omitting this, leaving the comment in for the time being though
                #if operon == 'casABCDE12':
                #    print 'AFDFAFDSF'
                if operon in result.keys():    # check if the operon is in the result dict already, if not make a new entry in the else clause
                    if accession in result[operon].keys(): # Check if the organism has been added to the operon
                        result[operon][accession].append(hlog.to_file())
                    else: # if the organims is not part of the operon, add it
                        result[operon].update({accession:[hlog.to_file()]})
                else: # add the operon to the result
                    result.update({operon: {accession: [hlog.to_file()]}})
            except:
                pass
            
    print sorted(result.keys()), len(result)
    
    # For the time being, i am going to cause each intermediate step in this pipeline to save in a folder called intermediate_for_debug
    intermediate_folder = './intermediate_for_debug/'
    if not os.path.isdir(intermediate_folder):
        os.makedirs(intermediate_folder)
    
    # For this step i will save the result in 'unfiltered_operon/'
    
    unfilter_folder = 'unfiltered_operon/'
    if not os.path.isdir(intermediate_folder + unfilter_folder):
        os.makedirs(intermediate_folder + unfilter_folder)
        
    for operon in result.keys():
        
        # this code is omitted because it used to debugging purposes, and is currently unneeded
        '''
    	outfile = intermediate_folder + unfilter_folder + operon + '.txt'
        #print "outfile", outfile
        handle = open(outfile, 'w')
        
        for accession in result[operon].keys():
        	handle.write('\n'.join(result[operon][accession]) + '\n')
        handle.close()
        '''
        
        # save results where i actually want them to go:
        #print "outfile", out_folder + operon + '.txt'
        handle = open(out_folder + operon + '.txt', 'w')
        for accession in result[operon].keys():
        	handle.write('\n'.join(result[operon][accession]) + '\n')
        handle.close()

# I have to figure out a better name for this function. The jist of what I am doing here is as follows:
# First, I will provide a file name that contains all of the hits for every organism that we are interested in.
# Then it sorts this homolog list first by organism, then by locus. (By the required input, the files already have 
# been screened for both eval cutoff and operon membership. The function will then return a dictionary for that
# operon. The intention is that this data structure will then be used to find the best hit for the locus out of the
# many redundant hits, however this functionality will be handled another function that i have yet to write/test.
def return_operon_list(fname):
    operon = fname.split('/')[-1].split('.')[0]
    hlog_list = [Homolog.from_file(i.strip()) for i in open(fname).readlines()]
    result_dict = {}
    for hlog in hlog_list:
        accession = hlog.accession()
        locus = hlog.locus()
        if accession not in result_dict.keys():
            result_dict.update({accession: {}})
        if locus not in result_dict[accession].keys():
            result_dict[accession].update({locus: [hlog]})
        else:
            result_dict[accession][locus].append(hlog)
    #print result_dict[accession]
    
    return result_dict
        

# might  not use this: will see
def parallel_return_operon_list(infolder, outfolder, num_proc):
    pool = Pool(processes = num_proc)
    organism_dict_for_recovery = dict(pool.map(parallel_operon_fasta, genome_of_interest_list))


# This function will take the organism-locus dict (per operon file) and determine the best homolog.
# 
def best_homolog_list(operon_dict, outfile):
    result = []
    
    for org in sorted(operon_dict.keys()):
        for locus in sorted(operon_dict[org].keys()):
            hlog_list = operon_dict[org][locus][1:]
            best_hit = operon_dict[org][locus][0]
            gene_count = defaultdict(int) # i am goign to use this, to see if a locus has more than one predicted gene, and the count ratio
            gene_count[best_hit.predicted_gene()] +=1
            for hlog in hlog_list:
                gene_count[hlog.predicted_gene()] +=1
                if best_hit.e_val() > hlog.e_val():
                    best_hit = hlog
            print gene_count.keys()
            result.append(best_hit)
    handle = open(outfile, 'w')
    handle.write('\n'.join([i.ret_str() for i in result]))
    handle.close()


def return_gene_to_operon_dict(fname):
    operon_dict = {}
    for line in [i.strip().split('\t') for i in open(fname).readlines()]:
        operon = line[0]
        for gene in line[1:]:
            operon_dict.update({gene: operon})
    return operon_dict
 
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, operon_query, filter_file, num_proc = check_options(parsed_args)
    
    print infolder, outfolder, operon_query, filter_file, num_proc
    
    # This code makes a dictionary mapping gene annotation to the operon that it belong to
    operon_dict = return_gene_to_operon_dict(operon_query)
    #print operon_dict
    
    #parallel_blast_parse_dict('./blast_parse/organism_raw_info/', './blast_parse/filtered_homologs/', num_proc, './genbank_pathway_lists/nc_filter_file.txt', operon_dict)
    parallel_blast_parse_dict(infolder, outfolder, num_proc, filter_file, operon_dict)
    
    
    #operon_dict = return_operon_list('./blast_parse/filtered_homologs/atpIBEFHAGDC.txt')
    
    #best_homolog_list(operon_dict, './blast_parse/processed_operon_files/atpIBEFHAGDC.txt')
    
    print time.time() - start

    # ./blast_parse.py -f phylo_order.txt
if __name__ == '__main__':
    main()
