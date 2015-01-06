#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
import math
from homolog4 import *
from collections import defaultdict
import itertools

# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will calculate the number of events that an organismal pair do not have in common.  The return is a pickled dictionary.")
    
    parser.add_argument("-i", "--infolder", dest="infolder", default='./optimized_operon/', metavar="FOLDER",
                help="A folder that contains the final operon results. The files will have operons grouped by organism, arranged by start, and have spurious BLAST results removed.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./operon_distance_matrices/',
                help="Folder where the filtered results will be stored. Default is the folder './operon_distance_matrices/'.")

    parser.add_argument("-F", "--operon_filter", dest="operon_filter", default='NONE', metavar="FILE",
                help="A file that contains the operons that are under investigation.  All others will be omitted from analysis an results.")            
    
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
     
    if os.path.exists(parsed_args.operon_filter):
        operon_filter_file = parsed_args.operon_filter
    elif parsed_args.operon_filter == 'NONE':
        operon_filter_file = parsed_args.operon_filter
    else:
        print "The file %s does not exist." % parsed_args.operon_filter
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    return infolder, outfolder, operon_filter_file, num_proc

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
    if filter_file == 'NONE':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]

    
# This function will take the folder of the operon result files, and it will return a dictionary that has as a primary key
# the operon name, the secondary key will be nc. The values of the attributes will then follow.
def parse_operon_result_files(in_folder, dest_folder, operon_filter_file):
    file_list = return_file_list(in_folder, operon_filter_file)
    
    print "file_list", file_list
    
    result = {}
    #distmat_dict = {}
    
    iddo_result = {}
    
    #gene_dict = make_operon_gene_string_dict(operon_file)
    #print len(gene_dict)
    
    print "in_folder", in_folder

    #for i in open(distmat_file).readlines():
    #    nc = i.strip().split('\t')[0]
    #    val = i.strip().split('\t')[1].strip() # fixes a wonky error in distmat.... ugh
    #    distmat_dict.update({nc: val})
    for fname in file_list:
        operon = fname.split('/')[-1].split('.')[0]
        #fname = in_folder + f
        print "fname", fname, operon
        result.update({operon: {}})
        
        iddo_result.update({operon: {}})
        
        file_output = []
        summary_info = [i.strip().split('\t')[1:] for i in open(fname).readlines() if i[:2] == '##']
        homolog_entries = []
        
        print "summary_info", summary_info
        # This section of code has a single purpose. ad-hoc the goodness of rearrangements and distance
        # between grouped genes.
        
        summmary_info = []
        tmp_hlog_list_for_grouping = []
        ignore_list = ['#', '$', '@', '+']
        #for i in [i.strip() for i in open(fname).readlines() if len(i) > 1]: # ugh, corrects for some garbage in the files (occurs once it seems)
        for i in [i.strip() for i in open(fname).readlines() if i.split('\t')[0] == '##']:
        #for i in [i.strip() for i in open(fname).readlines() if len(i) > 1 and i[0] == '#']:
            #if i[:2] == '##':
            if len(i) < 2:
                print "err", i, fname
            if i[0] == '#':
                comprehensive_list, group_list = group_homologs(tmp_hlog_list_for_grouping, INTERGENIC_MAX_LENGTH)
                for group in group_list:
                    #print gene_dict[operon]['reference_string']
                    rearrangements = return_operon_string_distance(gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict))
                    print gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict), return_operon_string_distance(gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict))
                    print "vals recorded", gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict), return_operon_string_distance(gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict)), rearrangements
                    print "i line", i, len(i.split('\t'))
                    try:
                        a,nc,c,d,e,f,g,h,i,j,k = i.split('\t') # only interested in a few of the fields here: 
                        distmat_val = distmat_dict[nc]
                        print c
                        common = '_'.join(c.split('_')[:2])
                        print "common", common
                    except:
                        print "Error in line", i, fname
                    
                    
                #print len(tmp_hlog_list_for_grouping), tmp_hlog_list_for_grouping, len(group_list)
                tmp_hlog_list_for_grouping = []
            elif i[0] not in ignore_list:
                tmp_hlog_list_for_grouping.append(return_homolog(i))
                pass
                #print i.split('\t')
        
        
        
        ignore_list = ['#', '$', '@', '+']
        
        #for line in [i.strip() for i in open(fname).readlines() if i[0] not in ignore_list ]:
        #    print line 
        #print homolog_entries
        #print summary_info
        print "summary_info2", summary_info
        for item in summary_info:
            nc = item[0]
            common = ' '.join(item[1].split('_')[:2])
            if nc in distmat_dict:
                distmat_val = distmat_dict[nc]
                #print nc, common
                line = [nc, common, distmat_val] + [i.strip() for i in item[2:]]
                #print line
                file_output.append(line)
        header = 'NC_Number,Common,Distance,Splitting,Migration,Duplicates,Deletions,Splits,Longest_Group,Total_Fusions,Group_Fusions\n'
        handle = open(dest_folder + operon + '.csv', 'w')
        handle.write(header)
        handle.write('\n'.join([','.join(i) for i in file_output]))
        handle.close()
    #print result.keys()    
    
    
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, operon_filter_file, num_proc = check_options(parsed_args)
    
    print infolder, outfolder, operon_filter_file, num_proc
    
    #file_list = return_file_list(infolder, opeon_filter_file)
    #print "file_list", file_list
        
    #parallel_list_param = [(i, outfolder, max_gap, e_val) for i in file_list]
    
    parse_operon_result_files(infolder, outfolder, operon_filter_file)

    print time.time() - start


if __name__ == '__main__':
    main()   
