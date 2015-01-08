#!/usr/bin/python

import time
import os
import sys
import argparse

#TODO:  Rename this script, it's horrible!

# Copyright(C) 2014 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

#########################################################################################################################################
# I am putting some globals here, they are command line arguments for some of the scripts that we are using.  They are not				#
# important enough, at least at this time, to justify making them command line arguments for them.  This can be revised					#
# later, or changed by someone who cares too much about these trivial things. after we get everything running to our satisfaction.  	#
# Most likely all/almost all will be removed because it may tempt someone to ruin what already seems to be working well.				#
#########################################################################################################################################

# removed from regulondb_dl_parse.py as a command line param for this master script
regulon_url = 'http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt'
regulon_outfolder = './regulonDB/'
# the followoing two require additional code to work, fix later
regulon_download = 'True'
regulon_experimental_only = 'True'

# removed from format_db.py as a command line param for this master script
format_protein = 'True'
BLAST_database_folder = './db/'

# removed from make_operon_query.py as a command line param for this script
refrence_organism = 'NC_000913'
operon_query_outfile = './operon_query.fa'


# removed from blast_script.py as a command line param for this script
blast_outfolder = './blast_result/'

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():
    
    parser = argparse.ArgumentParser(description='The purpose of this script is to run the full software suite that we have developed to study operons using as few inputs as possible.  This will facilitate the ease of use as much as possible.')

    parser.add_argument("-i", "--infile", dest="infile", metavar="FILE", default='./regulonDB/operon_names_and_genes.txt',
                help="Input file for the operon query step of the pipeline.")
    
    parser.add_argument("-I", "--infolder", dest="infolder", metavar="FOLDER", default='./genomes/',
                help="Folder containing all genbank files for use by the program.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./regulonDB/',
                help="Folder where results will be stored.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='NONE',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-m", "--min_genes", dest="min_genes", metavar="INT", default = 5, type=int,
                help="Minum number of genes that an operon must contain before it can be considered for further analysis. The default is 5 because that is what we are currently using in the study.")
                       
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
    if parsed_args.outfolder[-1] != '/':
       outfolder = parsed_args.outfolder + '/'
    else:
       outfolder = parsed_args.outfolder
    
    if parsed_args.filter == 'NONE' or os.path.exists(parsed_args.filter):
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
    
    if parsed_args.min_genes <= 0:
        min_genes = 1
    else:
        min_genes = parsed_args.min_genes

    return infolder, outfolder, filter_file, num_proc, min_genes

    
    
def main():
    
    start = time.time()

    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc, min_genes = check_options(parsed_args)
    
    print infolder, outfolder, filter_file, num_proc, regulon_download, regulon_url, regulon_experimental_only, min_genes
    
    '''
    # Stage 1: Get operon set and parse into something that we can use
    cmd1 = "./regulondb_dl_parse.py -f %s -i %s -o %s -n %i -u %s -m %i" % (filter_file, infolder, regulon_outfolder, num_proc, regulon_url, min_genes)
    # print "cmd1", cmd1
    os.system(cmd1)

    #Stage 2: Create BLAST searchable databases. (I am limiting this to protein databases right now since that is what we do in the paper)
    cmd2 = "./format_db.py -f %s -i %s -o %s -n %i" % (filter_file, infolder, BLAST_database_folder,  num_proc)
    
    # Set the database formatting option[Protein or DNA], even though we don't use it
    if format_protein == 'True':
        pass
    else:
        cmd2 = cmd2 + ' -d'
    print cmd2
    os.system(cmd2)
    
    #Stage 3: make the operon query fasta file(s)
    # This is the default location for where the parsed regulonDB file will be stored, and is a necessary input.
    # TODO: change the regulonDB script to allow the explicit naming of the desired output file. I can make this non-selectable
    # as well since no one should care, like ever about something so trivial.
    operon_file = regulon_outfolder + 'operon_names_and_genes.txt'
    
    
    cmd3 = "./make_operon_query.py -i %s -o %s -p %s -n %i -r %s" % (infolder, operon_query_outfile, operon_file, num_proc, refrence_organism)
    print cmd3
    os.system(cmd3)

    #Stage 4: run BLAST with the query that we made in stage 3, using the databases that we used in stage 2.
    # TODO: add eval filtering here, going with default since i'm low on time.  i will fix in the nex few days
    cmd4 = "./blast_script.py -d %s -o %s -f %s -n %i -q %s" % (BLAST_database_folder, blast_outfolder, filter_file, num_proc, operon_query_outfile)
    print cmd4
    os.system(cmd4)
    
    # Stage 5: Parse the BLAST result and sort it by gene block
    
    # i'm just trying to get this out the door, everything works how it should, but i am saving time to get this out the door. the final
    # version will implement some ability to control this program's behavior.
    cmd5 = "./blast_parse.py -f %s -n %i" % (filter_file, num_proc)
    print cmd5
    os.system(cmd5)
    
    '''
    
    # Stage 6: filter out spurious results and report the gene blocks that best represent the origional.
    cmd6 = "./filter_operon_blast_results.py -n %i" % (num_proc)
    print cmd6
    os.system(cmd6)
    
    print time.time() - start
    
if __name__ == '__main__':
    main()
