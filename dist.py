#!/usr/bin/python

# I do not know how many of these i need, i will edit later
import Bio
import re
import os
import sys
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Application
from Bio.Application import _Option
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
import math
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import math
from scipy import stats
from Levenshtein import distance
import cPickle as pickle
import itertools

from homolog import *  #allows me to use "Homolog" without extra refrencing



# Globals - I will likely have more platform-specific tweaks as time progresses
if os.name == "posix":
    gl_new_line = '\n'
else:
    gl_new_line = '\r\n'
    
NUM_PROCESSORS = os.sysconf("SC_NPROCESSORS_CONF") - 1
INTERGENIC_MAX_LENGTH = 200

# not really necessary for this program!!!!!! This should be something that is known before this piece of code is run.
# What this functions does is to search the genbank_source_directory for files with a specified criteria that will
# identify chromosomal genbank files and exclude plasmids or any other extra-chromosomal genetic elements. Or at least
# reasonable likely to only be a chromosomal file... turns out there is no standard verbage for chromosome files, ugh.
def find_complete_genome_genbank_files(genbank_source_directory, result_file):
    file_list = []
    #grep_command = 'grep -l -i -e "complete *genome" `find %s -name "*.gbk"` > %s' % (genbank_source_directory, result_file)
    grep_command = 'grep -l -i -e "complete" `find %s -name "*.gbk"` > %s' % (genbank_source_directory, result_file)
    print "grep_command ", grep_command
    os.system(grep_command)

    #handle = open(fname, 'r')
    for rec in [i.strip() for i in open(result_file).readlines()]:
        temp = rec.split('/')
        # the file will be returned as a tuple, the first value will be the folder within the genbank 
        # source directory.  The second will be the genbank file name.  
        #file_name = temp[len(temp)-1].split('\n')[0]
        file_list.append((temp[len(temp)-2], temp[len(temp)-1].strip()))
        #file_list = 

    return file_list


# Find the marker genes of interest (located in the marker file) in a reference organism file (genbank format),
# and create a fasta file to use as a query for blast search.
def make_marker_fasta(reference_organism_file, marker_list, result_file):
    err = False # returns error status to calling function
    genes_found = []
    result = []
    seq_record = SeqIO.parse(open(reference_organism_file), "genbank").next()
    accession = seq_record.annotations['accessions'][0]
    organism = seq_record.annotations['organism'].replace(' ', '_')
    for fnum, feature in enumerate(seq_record.features):            
        if feature.type == 'CDS':
	    start = feature.location._start.position
            stop = feature.location._end.position
            try: 
                gene = feature.qualifiers['gene'][0]
                gene = gene.lower()
            except:
                gene = 'unknown'
            if gene in marker_list:
                genes_found.append(gene)
                seq = feature.qualifiers['translation'] # this returns the protein product, not suitable for RNA products like 16s
                result.append(">%s|%s|%s" % (organism, accession, gene))
                result.append(''.join(seq))
    handle = open(result_file, 'w')
    handle.write('\n'.join(result) + '\n')
    handle.close()
    if set(marker_list).issubset(set(genes_found)):
        err = False
    else:
        print "Reference organism is missing marker genes: %s" % ', '.join(list(set(marker_list) - set(genes_found)))
        err = True
    return err

# This function takes a file containing the list of organisms that are of interest, and a file of the full path
# to organism files on the local machine, and returns a list of the match organism file to the full path.
def return_full_path(org_file, path_file):
    path_dict = {}
    org_list = [i.strip() for i in open(org_file).readlines()]
    path_list = [i.strip() for i in open(path_file).readlines()]
    for item in path_list:
        key = item.split('/')[len(item.split('/')) -1].split('.')[0]
        path_dict.update({key:item})
    result = []
    for org in org_list:
        if org in path_dict:
            result.append(path_dict[org])
        else:
            print 'Missing: ', org
    return result


def create_distmat(fname, method = 1):
    cline = ClustalwCommandline('clustalw', infile=fname)
    #print cline
    base = fname.split('/')[len(fname.split('/')) - 1].split('.')[0]
    better = fname.split('.')[0]
    #print 'better: ', better
    return_code = subprocess.call(str(cline), stdout = open(os.devnull), stderr = open(os.devnull), shell=(sys.platform!="win32"))
    #print "return_code", return_code
    distmat_line = "distmat %s.aln -outfile %s.distmat  -protmethod %i" % (better, better, method) 
    #print 'distmat_line: ', distmat_line
    return_code = subprocess.call(distmat_line, stdout = open(os.devnull), stderr = open(os.devnull), shell=(sys.platform!="win32"))


# this function will take a list of marker genes, as well as a file containing a list of the NC numbers of organisms that we
# are interest in, and a file that contains the full paths to all genbank files we have.
# The return will be marker gene files... and i need to prohibit multiple copies of the genes.
def make_target_fasta(marker_list, org_file, path_file):
    #print "Got here make_target_fasta"
    org_paths = return_full_path(org_file, path_file)
    result = [] # this is not great
    for org in org_paths:
        genes_found = []
        seq_record = SeqIO.parse(open(org), "genbank").next()
        accession = seq_record.annotations['accessions'][0]
        organism = seq_record.annotations['organism'].replace(' ', '_')
        for fnum, feature in enumerate(seq_record.features):            
            if feature.type == 'CDS':
                start = feature.location._start.position
                stop = feature.location._end.position
                try: 
                    gene = feature.qualifiers['gene'][0]
                    gene = gene.lower()
                except:
                    gene = 'unknown'
                if gene in marker_list:
                    genes_found.append(gene)
                    seq = feature.qualifiers['translation'] # this returns the protein product, not suitable for RNA products like 16s
                    result.append(">%s|%s|%s" % (accession, organism, gene))
                    result.append(''.join(seq))
        if not set(marker_list).issubset(set(genes_found)):
            pass
            #print "%s, %s is missing marker genes: %s" % ( organism, accession, ', '.join(list(set(marker_list) - set(genes_found))))
    # to complete this piece of code quickly, i am shortcutting this step. ideally we will store (into a dict) the genes from each organism keyed on the marker gene.
    # When we miss a gene in an organism, we will then BLAST for it, and if we hit something, then add it. I am skipping this. I am not making this function take more
    # than one marker gene at a time either. I just want something fast for the lab meeting.
    # this whole thing should be iteratable, etc..... ugh time constraints
    
    handle = open("%s_for_msa.fa" % marker_list[0], 'w')
    handle.write('\n'.join(result))
    handle.close()
    create_distmat("%s_for_msa.fa" % marker_list[0])

# returns a datastructure that has constant time access to distmat data   
# the keying on this will be the NC number  
def ReadDistmat(fname, delim = '\t'):
    result = {}
    cnt = 0
    organism = []
    temp = []
    handle = open(fname, 'r')
    for rec in handle:
        #print rec
        if cnt > 6:
            organism.append(rec.split(delim)[len(rec.split(delim)) - 1 ].split('|')[0])
            junk = rec.split('\n')[0].split(delim)[1:]
            hold = []
            for i in junk:
                hold.append(i.strip())
            #temp.append(rec.split('\n')[0].split(delim)[1:])
            temp.append(hold)
            #print organism
        cnt = cnt + 1
    #print 'length: ', len(organism) , ' ' ,organism
    #print 'temp: ', temp
    #print temp[2]
    # so outer and inner are with respect to the keys that are eithr in the outer or inner loops
    for i in range(0,len(temp)):
        outer = organism[i]
        result.update({outer:{}})
        for j in range(i, len(temp)):
            inner = organism[j]
            #print outer, ' ', inner
            value = temp[i][j]
            #print value
            t = result[outer]
            t.update({inner:value})
            result[outer].update(t)
            if inner not in result:
                result.update({inner:{}})
            t = result[inner]
            t.update({outer:value})
    return result

def MakeSingleOrganismDistanceFile(fname, nc_number, save_file, delim = '\t'):
    dist_dict = ReadDistmat(fname)
    handle = open(save_file, 'w')
    #print 'Length of the organism distance dict: ', len(dist_dict[nc_number])
    for key in dist_dict[nc_number]:
        #print dist_dict[nc_number][key]
        write_str = key + delim + dist_dict[nc_number][key]
        handle.write(write_str + '\n')
        
def read_distmat(distmat_file, org_file, out_file):
    file_as_list = [i.strip() for i in open(distmat_file).readlines()][7:]
    ref_line = file_as_list[0].split('\t')
    nc_to_index_dict = {}
    for item in file_as_list:
        tmp = item.split('\t')[len(item.split('\t')) -1]
        index = tmp.split(' ')[1]
        nc = tmp.split('|')[0]
        nc_to_index_dict.update({nc:int(index) - 1})
    #print nc_to_index_dict['NC_003198']
    result = []
    # this is here to preserve the order of the organism file only, can be omitted if order is irrelevant
    org_list = [i.strip() for i in open(org_file)]
    for org in org_list:
        if org in nc_to_index_dict:
            index = nc_to_index_dict[org]
            val = ref_line[index]
            result.append("%s\t%s" % (org, val))
    handle = open(out_file, 'w')
    handle.write('\n'.join(result))
    handle.close()

def return_group_str(group_homologs, operon, gene_dict):
    result = ''
    for h_log in group_homologs:
        gene = h_log.predicted_gene()
        result = result + gene_dict[operon][gene]
    return result# , gene_dict[operon]['refrence_string']
    #return result


def return_operon_string_distance(operon_string, gene_string):
    len_operon = len(operon_string)
    len_gene_group = len(gene_string)
    length_difference = len_operon - len_gene_group
    
    reverse_gene_string = gene_string[::-1]
    
    d1 = distance(operon_string, gene_string) - length_difference
    d2 = distance(operon_string, reverse_gene_string) - length_difference
    
    return min(d1, d2)



# This function will take the folder of the operon files, and it will return a dictionary that has as a primary key
# the operon name, the secondary key will be nc. The values of the attributes will then follow. I will post order
# after i get this up and running.
def parse_operon_result_files(in_folder, distmat_file, dest_folder, operon_file):
    file_list = os.listdir(in_folder)
    result = {}
    distmat_dict = {}
    
    iddo_result = {}
    
    gene_dict = make_operon_gene_string_dict(operon_file)
    print len(gene_dict)
    
    print "in_folder", in_folder

    for i in open(distmat_file).readlines():
        nc = i.strip().split('\t')[0]
        val = i.strip().split('\t')[1].strip() # fixes a wonky error in distmat.... ugh
        distmat_dict.update({nc: val})
    for f in file_list:
        operon = f.split('.')[0]
        fname = in_folder + f
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
                print "errr", i, fname
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

# This function assumes that Distance will be at index number 2
def return_operon_line(fname, significance = .05, min_number_operon_organisms = 3):
    #operon = fname.split('/')[len(fname.split('/'))-1].split('.')[0]
    operon = fname.split('/')[-1].split('.')[0]

    #print fname
    success = True
    result = {}
    res_list = []
    index_dict = {}
    res = [i.strip().split(',') for i in open(fname).readlines()]
    #print len(res)
    if len(res) <= (min_number_operon_organisms):
        #print operon, 'crap', len(res)
        success = False
        return (False, '')
        
    else:
        print operon, 'good', len(res)
        cnt = 0
        for item in res[0]:
            index_dict.update({cnt:item})
            result.update({item:[]})
            cnt +=1
        for item in res[1:]:
            cnt = 0
            #for i in [int(math.fabs(j)) for j in item]:
            for i in item:
                #print 
                index = index_dict[cnt]
                result[index].append(i)
                cnt +=1
        #print result
        binary_vector = []
        for index in range(3, len(res[0])):
            name = index_dict[index]
            #print "name", name
            #print result[index_dict[index]]
            if sum([float(i) for i in result[index_dict[index]]]) == 0:
                #print 'No statistics available, no observation set for', name
                p_val = 'None'
                spear_r = 'None'
                sig = 'None'
                binary_vector.append('N')
                res_list.append('\t'.join([str(i) for i in ['##', operon, name, spear_r, p_val, sig]]))
            else:
                #print name, result['Distance'], result[index_dict[index]]
                #tmp = [str(math.fabs(float(i))) for i in result[index_dict[index]]]
                #tmp = [type(i) for i in result[index_dict[index]]]
                #tmp = [str(math.fabs(float(i))) for i in result['Distance']]
                #print tmp
                spear_r, p_val = stats.spearmanr(result['Distance'], result[index_dict[index]])
                #print stats.spearmanr(result['Distance'], result[index_dict[index]])   
                if p_val < significance:
                    sig = 'Yes'
                    binary_vector.append(1)
                else:
                    sig = 'No'
                    binary_vector.append(0)
                res_list.append('\t'.join([str(i) for i in ['##', operon, name, spear_r, p_val, sig]]))
            #spear_r, p_val = stats.spearmanr(result['Distance'], result[index_dict[index]])
            #print operon, name, spear_r, p_val, sig
            #res_list.append('\t'.join([str(i) for i in ['##', operon, name, spear_r, p_val, sig]]))
        #res_list.append("$$\t%s\t" % operon + ','.join([str(i) for i in binary_vector]))
        tmp = ['$$', operon, len(result['Distance']), ','.join([str(i) for i in binary_vector])]
        res_list.append('\t'.join([str(i) for i in tmp]))
    #print res_list
    #print stats.spearmanr(result['Distance'], result['Splitting'])
    #print result['Distance']
    return (True, '\n'.join(res_list))
    
        
# This function will calculate the spearman rank coefficient

def calculate_spearman(folder, outfile, significance):
    file_list = sorted(os.listdir(folder))
    result = []
    header = 'Tag   Operon  Attribute   r-val   p-val   Significant'

    for fname in [folder + i for i in file_list]:
        success, line = return_operon_line(fname, significance)
        if success:
            result.append(line)
        #result.append(return_operon_line(fname, significance))
    handle = open(outfile, 'w')
    handle.write(header + '\n')
    handle.write('\n'.join(result))
    handle.close()


# this program will create a summary table of all operon data by organism
# the 

def create_organism_summary(folder, outfile = 'tally_table.csv', columns_of_interest = [5,6,7,9]):
    #print "got here"
    file_list = sorted(os.listdir(folder))
    #print file_list
    first_time = True
    header = ''
    # dictionary used to keep tallies on totals. keyed on the NC number
    # will contain nc_number, common_name, and distance in the first field (a list) that is not modified
    # the second field will contain a list of summed individual markers. a running total will be contained
    # in the final field (single element list) because it is easy later
    organism_dict = {}
    for fname in file_list:
        # making the headers for the summary file
        #print fname
        if first_time:
            # the first 3 colums have identification data and the distance, which are necessary, but not tallied.
            header_list = [0,1,2] + columns_of_interest
            line = [i.strip() for i in open(folder+fname, 'r').readlines()][0].split(',')
            #print "folder+fname", folder+fname
            
            header_tmp = []
            for i in header_list:
                #print i
                print "line", line, "fname", fname
                header_tmp.append(line[i])
            header_tmp.append('Total')
            header = ','.join(header_tmp)
            #print  'header:', header
            first_time = False
        if not first_time: # no more setup for the summary file needed, do real work
            for item in [i.strip().split(',') for i in open(folder+fname, 'r').readlines()[1:]]:
                #print item
                nc_number = item[0]
                common_name = item[1]
                distance = item[2]
                indentification = [nc_number, common_name, distance]
                vals_to_sum = []
                for index in columns_of_interest:
                    vals_to_sum.append(int(item[index]))
                if nc_number not in organism_dict.keys(): # add new dictionary element
                    organism_dict.update({nc_number:[indentification, vals_to_sum, [sum(vals_to_sum)]]})
                else:
                    ident, vals, total = organism_dict[nc_number]
                    tmp = []
                    for i in range(0, len(vals)):
                        tmp.append(vals[i] + vals_to_sum[i])
                    total = [sum(tmp)]
                    organism_dict[nc_number] = [ident, tmp, total]
                    
                #print vals_to_sum
        #print organism_dict
    res = []
    total_line = []
    for n in range(0,len(vals) + 1):
        total_line.append(0)
    #print "total_line", total_line
    for org in sorted(organism_dict.keys()):
        ident, vals, total = organism_dict[org]
        #print ident + [str(i) for i in vals + total]
        res.append(ident + [str(i) for i in vals + total])
        total_line = [(a + b) for a, b in zip(vals + total, total_line)]
    handle = open(outfile, 'w')
    handle.write(header + '\n')
    handle.write('\n'.join([','.join(i) for i in res]) + '\n')
    handle.write(',,Totals,' + ','.join([str(i) for i in total_line]))
    handle.close()
         


def calculate_cost(folder, processed_operon_file, atribute_list = ['Splits', 'Duplicates', 'Deletions', 'Total_Fusions']):
    # this section of code creates a dictionary to directly access the length of an operon based on its name
    #print folder
    print processed_operon_file
    operon_length_dict = {}
    for line in [i.strip().split('\t') for i in open(processed_operon_file).readlines()]:
        operon = line[0]
        gene_list = line[1:]
        number_genes = len(gene_list)
        operon_length_dict.update({operon:{'number_genes':number_genes, 'gene_list':gene_list}})

    # this code block will ultimately build the data structure that allows the calculation Po/Pe

    # determine all files in the directory where observations are summarized
    file_list = sorted(os.listdir(folder))
    #print file_list

    # dict keeping track of observations, and the possible number of observations for each attribute
    result_dict = {'total_specific':{}, 'grand_total':{'observed':0, 'possible':0}}
    #print atribute_list
    for att in atribute_list:
        result_dict['total_specific'].update({att:{'observed':0, 'possible': 0}})
    #print 'result', result_dict
    
    #print operon_length_dict, len(operon_length_dict)
    
    operons_investigated = operon_length_dict.keys()
    #print operons_investigated
    #for file_name in file_list:
    for file_name in [i for i in file_list if i.split('.')[0] in operons_investigated]:
        operon = file_name.split('.')[0]
        #print operon, file_name, operon_length_dict[operon]
        fname = folder + file_name
        
        # all the lines in the operon-specific csv file, will be used for quick calculations for cost analysis
        lines = [i.strip().split('\t') for i in open(fname).readlines()]
        number_of_orgs = len(lines) - 1
        #print "lines", lines
        #print 'number_of_orgs', number_of_orgs
        header = lines[0][0]
        #print "header", header
        column_dict = {}
        cnt = 0
        #print 'Ahh', [i.split(',') for i in header]
        for column_name in header.split(','):
            column_dict.update({column_name:cnt})
            cnt += 1
        #print column_dict
        for line in [i.strip().split('\t') for i in open(fname).readlines()[1:]]:
            tmp = line[0].split(',')
            for attribute in atribute_list:
                att_index = column_dict[attribute]
                att_value = int(tmp[att_index])
                #if attribute == 'Splits' or attribute == 'Total_Fusions': # you can split things in a group one fewer time than there are total genes
                if attribute == 'Splits' or attribute == 'Deletions': # you can split things in a group one fewer time than there are total genes
                    #possible_observations = int(operon_length_dict[operon]['number_genes'] - 1)
                    possible_observations = int(operon_length_dict[operon]['number_genes'] - 2)
                elif attribute == 'Total_Fusions': # you can fuse genes together 
                    possible_observations = int(operon_length_dict[operon]['number_genes'] - 1)
                elif attribute == 'Inversions':
                    possible_observations = int((operon_length_dict[operon]['number_genes'] - 1)/2)
                    #print "Inversions", operon_length_dict[operon]['number_genes'],  possible_observations
                else:
                    possible_observations = int(operon_length_dict[operon]['number_genes'])
              
                # now it is time to take the values we extract, and places them into the correct place in the result dict
                if operon not in result_dict.keys():
                    result_dict.update({operon:{'total_observed': 0, 'total_possible': 0, 'num_orgs':number_of_orgs}})
                    #print result_dict
                if attribute not in result_dict[operon].keys():
                    tmp_tot = result_dict[operon]['total_observed'] + att_value
                    tmp_possible = result_dict[operon]['total_possible'] + possible_observations
                    result_dict[operon].update({attribute:{'operon_specific_observations': att_value, 'possible_operon_observations': possible_observations}, 'total_observed': tmp_tot, 'total_possible': tmp_possible})
                else:
                    # this is the verbose, but easy to understand implementation. add the total to the new val, then replace the entry in the result dict
                    temp_specific = result_dict[operon][attribute]['operon_specific_observations'] + att_value
                    temp_possible = result_dict[operon][attribute]['possible_operon_observations'] + possible_observations
                    temp_total = result_dict[operon]['total_observed'] + att_value
                    temp_total_possible = result_dict[operon]['total_possible'] + possible_observations
                    result_dict[operon].update({attribute:{'operon_specific_observations': temp_specific, 'possible_operon_observations': temp_possible}, 'total_observed': temp_total, 'total_possible': temp_total_possible})
                
   
                # update the running total
                temp_total_observations = result_dict['grand_total']['observed'] + att_value
                temp_total_possible = result_dict['grand_total']['possible'] + possible_observations
                result_dict['grand_total'].update({'observed': temp_total_observations, 'possible': temp_total_possible})
                # update running totals for each specific attribute
                temp_observed_operon_total = result_dict['total_specific'][attribute]['observed'] + att_value
                temp_possible_operon_total = result_dict['total_specific'][attribute]['possible'] + possible_observations
                result_dict['total_specific'][attribute].update({'observed':temp_observed_operon_total, 'possible': temp_possible_operon_total})
        #print result_dict
        
        
    # ok now we are done making the dictionary that will solve everything! (in a limited context)
    # time to leverage the datastructure's information and report a nice summary
    # note that we are solving the equation: Po/Pe where:
    # Pe = probability of expected observation: total attribute observation/ number of possible observations
    # Po = probability of observed: observations in an operon/total possible observations possible in this operon
    
    total_observed_attributes = result_dict['grand_total']['observed']
    total_possible_observations = result_dict['grand_total']['possible']
    Pe_all_observations = float(total_observed_attributes)/float(total_possible_observations)
    
    # Header list so that i can add additional information without greatly modifying existing code
    h_list = ['Operon', 'Num_Orgs', "Length", 'Total_ratio']
    #print "result_dict.keys()", result_dict.keys()
    garbage = sorted(result_dict.keys())[0]
    ignore_list = ['total_possible', 'total_observed', 'num_orgs']
    
    for att in [i for i in result_dict[garbage] if i not in ignore_list]:
        h_list.append(att+'_ratio')
    header = ','.join(h_list)
    #print header
    
    result = [header]
    
    for operon in sorted([i for i in result_dict.keys() if i not in ['total_specific', 'grand_total']]):
        total_operon_observations = result_dict[operon]['total_observed']
        total_operon_observations_possible = result_dict[operon]['total_possible']
        Po_all_operon_observations = float(total_operon_observations)/float(total_operon_observations_possible)
        ratio_whole_operon = Po_all_operon_observations/Pe_all_observations
        #print operon, total_observed_attributes, total_possible_observations, Pe_all_observations, total_operon_observations, total_operon_observations_possible, Po_all_operon_observations, ratio_whole_operon
        
        number_of_genes = str(operon_length_dict[operon]['number_genes'])

        number_of_organisms = str(result_dict[operon]['num_orgs'])
        
        # This list will be used to create each line of the summary
        line = [operon, number_of_organisms, number_of_genes, "%5.3f" % math.log(ratio_whole_operon)]
        for att in [i for i in result_dict[operon].keys() if i not in ignore_list]:
        
            att_total_observed = result_dict['total_specific'][att]['observed']
            att_total_possible = result_dict['total_specific'][att]['possible']
            # Pe
            att_total_ratio = float(att_total_observed)/float(att_total_possible)
            
            att_observed = result_dict[operon][att]['operon_specific_observations']
            att_possible = result_dict[operon][att]['possible_operon_observations']
            # Po
            att_ratio = float(att_observed)/float(att_possible)
            
            # Po/Pe
            if att_ratio == 0:
                prob_ratio = 0
            else:
                prob_ratio = math.log(att_ratio/att_total_ratio)
            
            line.append("%5.3f" % prob_ratio)
        #print line
        #result.append(','.join([str(i) for i in line]))
        result.append(','.join(line))
    
    handle = open('ratio_summary.csv', 'w')
    handle.write('\n'.join(result))
    handle.close()
            
            
                
# returns a Homolog object from a line in a file
# i think that i should make this a function in the homolog class, but for now it stays as a function here
def return_homolog(line):
    if len(line.strip().split('\t')) == 15:
        a,b,c,d,e,f,g,h,i,j,k,l,m,n,o = line.strip().split('\t')
        f_list = f.split(':')
        return Homolog(a,b,c,d,e,f_list,float(g),float(h),float(i),float(j),int(k),int(l),int(m),n,o)
    else:
        a,b,c,d,e,f,g,h,i,j,k,l,m,n = line.strip().split('\t')
        f_list = f.split(':')
        return Homolog(a,b,c,d,e,f_list,float(g),float(h),float(i),float(j),int(k),int(l),int(m),n)

def parse_operon_name(operon):
    result = []
    group = operon.split('-')
    for i in group:
       prefix = i[:3]
       genes = i[3:]

       if len(genes) > 0:
           for gene in genes:
               result.append(prefix + gene)
       else:
           result.append(prefix)
    return result

    
# The point of the function is to create a unique string of characters from a list of operon genes.
# This unique string will then be used to calculate the Levinstein edit distance for groups of homologs
# with respect the operon in E.coli. 
def make_operon_gene_string_dict(operon_file = './operon_name_and_genes.txt'):
    result = {}
    
    print operon_file
    for line in [i.strip().split('\t') for i in open(operon_file).readlines()]:
        operon = line[0]
        print operon
        result.update({operon:{'reference_string': ''}})
        
        # Returns the genes in order and a corresponding index. This index will be used to generate the
        # unicode integer of capital letters (as they are lower numerically than lower case).
        
        operon_genes_in_order = parse_operon_name(operon)
        #for gene, index in zip(line[1:], range(0,len(line[1:]))):
        for gene, index in zip(operon_genes_in_order, range(0,len(operon_genes_in_order))):
            print gene
            result[operon].update({gene:chr(65+index)})
            result[operon].update({'reference_string': result[operon]['reference_string'] +  chr(65+index)})
        print operon, result[operon], len(result)
    return result
        
# this function will take a list of homologs and return a list of lists
# each element will be a list of grouped homologs. the criteria is that
# each item will be no further than intergenic_max_length away from its closest neighbor
# strand is completely ignored for this function. (Strand needs to be optimized before it can 
# be used as an optimality constraint)
def group_homologs(list_h, intergenic_max_length):
    list_homologs = sorted(list_h, key = lambda a: a.start())
    comprehensive_list = []
    local_group = []
    ungrouped = True
    for i in range(0, len(list_homologs)-1):
        #look at current
        start = list_homologs[i].start()
        stop = list_homologs[i].stop()
        # look at next
        start_n = list_homologs[i+1].start()
        stop_n = list_homologs[i+1].stop()
        
        #if math.fabs(start - stop_n) < intergenic_max_length or math.fabs(stop - start_n) < intergenic_max_length:
        #print math.fabs(stop - start_n)
        
        # Check for gene fusion
        if start == start_n:
            if ungrouped:
                local_group.append(list_homologs[i])
                local_group.append(list_homologs[i+1])
                ungrouped = False
            else:
                local_group.append(list_homologs[i+1])
        elif math.fabs(stop - start_n) < intergenic_max_length:
            if ungrouped:
                local_group.append(list_homologs[i])
                local_group.append(list_homologs[i+1])
                ungrouped = False
            else:
                local_group.append(list_homologs[i+1])
        else: # get ready for next possible set of matches (in a segmented OTU)
            ungrouped = True
            if len(local_group) > 0:
                comprehensive_list.append(local_group)
                local_group = []
            else:
                comprehensive_list.append([list_homologs[i]])
    if ungrouped:
        comprehensive_list.append([list_homologs[len(list_homologs)-1]])
    else:
        comprehensive_list.append(local_group)
    group_list = []
    
    # more readable version determining if there are groups of homologs in the list of results
    for i in comprehensive_list:
        if len(i) > 1:
            group_list.append(i)
    #print 'group_list', group_list
    #return comprehensive_list, [[j for j in i] for i in comprehensive_list if len(i) > 1] # this last one is only grouped homologs
    return comprehensive_list, group_list

def return_longest_group(lst):
    if len(lst) == 0:
        return 0
    else:
        return max([len(i) for i in lst])


'''
# this function is inappropriate for optimization of inversions. optimizing this param may be quite important
# though.  long story, inversions should not constantly occur when we reconstruct, but we canot know how to minimize
# them without trying to reconstruct the ancestry in the first place. 
def return_inversions_from_group_line(group_str):
    tmp = [i for i in group_str.split(' ') if len(i) > 0]
    print "tmp", tmp
    pos_strand = 0
    neg_strand = 0
    for group in tmp:
        if group[0] == '<':
            neg_strand += 1
        else:
            pos_strand += 1
    print "tmp", tmp, str(min(pos_strand, neg_strand))
    return str(min(pos_strand, neg_strand))
'''

# this function is inappropriate for optimization of inversions. optimizing this param may be quite important
# though.  long story, inversions should not constantly occur when we reconstruct, but we canot know how to minimize
# them without trying to reconstruct the ancestry in the first place. 
def return_inversions_from_group_line(group_str):
    tmp = [i for i in group_str.split(' ') if len(i) > 0]
    print "tmp", tmp
    pos_strand = 0
    neg_strand = 0
    
    #for group in tmp:
    #    if group[0] == '<':
    #        neg_strand += 1
    #    else:
    #        pos_strand += 1
    
    for group in tmp:
        strand = 'not_eval'
        for i in group:
            if i == '<': # start minus strand
                neg_strand +=1
                strand = 'minus'
            elif i == '>': # end neg strand, get ready to count something if it exists
                strand = 'not_eval'
            elif strand == 'not_eval':
                strand = 'positive'
                pos_strand += 1
            else:
                pass
    #print "pos_strand", pos_strand, "neg_strand", neg_strand
    print "tmp", tmp, str(min(pos_strand, neg_strand))
    return str(min(pos_strand, neg_strand))


# This function will take the folder of the operon files, and it will return a dictionary that has as a primary key
# the operon name, the secondary key will be nc. The values of the attributes will then follow. I will post order
# after i get this up and running.
def parse_operon_result_files2(in_folder, distmat_file, dest_folder, operon_file):
    file_list = os.listdir(in_folder)

    result = {}
    distmat_dict = {}
    
    iddo_result = {}
    
    common_result = {}
    
    gene_dict = make_operon_gene_string_dict(operon_file)
    print len(gene_dict)
    
    print "in_folder", in_folder

    for i in open(distmat_file).readlines():
        nc = i.strip().split('\t')[0]
        val = i.strip().split('\t')[1].strip() # fixes a wonky error in distmat.... ugh
        #print "val",val
        distmat_dict.update({nc: val})
    for f in file_list:
        operon = f.split('.')[0]
        fname = in_folder + f
        print "fname", fname, operon
        result.update({operon: {}})
        
        iddo_result.update({operon: {}})
        common_result.update({operon: {}})
        
        file_output = []
        file_output2 = []
        summary_info = [i.strip().split('\t')[1:] for i in open(fname).readlines() if i[:2] == '##']
        homolog_entries = []
        
        #print "summary_info", summary_info
        # This section of code has a single purpose. ad-hoc the goodness of rearrangements and distance
        # between grouped genes.
        
        summmary_info = []
        tmp_hlog_list_for_grouping = []
        ignore_list = ['$', '@', '+']
        print "fname", fname
        #for i in [i.strip() for i in open(fname).readlines() if len(i) > 1]: # ugh, corrects for some garbage in the files (occurs once it seems)
        for i in [i.strip() for i in open(fname).readlines()]:
            # ok back to stripping more information. we will do inversions.  this next piece of code has no other purpose
            # due to oversights on my part from a year ago, i have to do the following:
            # Parse the group data line, which is denoted with the "++" marker. Then determine the number of inversions. This will be accomplished by
            # a function that takes in the grouping string. The function will return the minimum value of inversions with respect to both strands.
            # this value will be saved in a variable, and when results are stored, will be appended at the end.
            if len(i) <= 2:
                print "input line has no data, this is a problems with file", fname
            elif i[0] == '+':
                g_line = i.split('\t')[1]
                inversions = return_inversions_from_group_line(g_line)
            elif i[0] == '#':
                comprehensive_list, group_list = group_homologs(tmp_hlog_list_for_grouping, INTERGENIC_MAX_LENGTH)
                rearrangement_list = []
                for group in group_list:
                    #print gene_dict[operon]['reference_string']
                    rearrangements = return_operon_string_distance(gene_dict[operon]['reference_string'], return_group_str(group, operon, gene_dict))
                    rearrangement_list.append(rearrangements)
                total_rearrangements = str(sum(rearrangement_list))
                try:
                    tmp = i.split('\t')[3:]
                    a,nc,c,d,e,f,g,h,z,j,k = i.split('\t') # only interested in a few of the fields here: 
                    common = '_'.join(c.split('_')[:2])
                    
                    # we don't have a distmat value for vibrio_cholerae, so in teh event we have the organism, in the dataset, ignore hte data. not a good long term solution mind
                    try:
                        distmat_val = distmat_dict[nc]
                    except:
                        pass
                except:
                    print "error in the house"
                    
                to_store = [nc, common, distmat_val] + tmp + [total_rearrangements] + [inversions]
                duplicates = to_store[5]
                deletions = to_store[6]
                splits = to_store[7]
                total_rearrangements = to_store[11]
                inversions = to_store[12]
                print "to_store", to_store
                file_output2.append([nc, common, distmat_val] + tmp + [total_rearrangements] + [inversions])  
                #iddo_result[operon].update({common:[duplicates, deletions, splits, total_rearrangements]})  
                iddo_result[operon].update({nc:[duplicates, deletions, splits, total_rearrangements, inversions]}) 
                common_result[operon].update({common:[duplicates, deletions, splits, total_rearrangements, inversions]}) 
                    
                #print len(tmp_hlog_list_for_grouping), tmp_hlog_list_for_grouping, len(group_list)
                tmp_hlog_list_for_grouping = []
            elif i[0] not in ignore_list:
                tmp_hlog_list_for_grouping.append(return_homolog(i))
        
        
        
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
        #header = 'NC_Number,Common,Distance,Splitting,Migration,Duplicates,Deletions,Splits,Longest_Group,Total_Fusions,Group_Fusions\n'
        header = 'NC_Number,Common,Distance,Splitting,Migration,Duplicates,Deletions,Splits,Longest_Group,Total_Fusions,Group_Fusions,Inversions\n'
        handle = open(dest_folder + operon + '.csv', 'w')
        handle.write(header)
        handle.write('\n'.join([','.join(i) for i in file_output]))
        handle.close()
        
        header2 = 'NC_Number,Common,Distance,Splitting,Migration,Duplicates,Deletions,Splits,Longest_Group,Total_Fusions,Group_Fusions,Total_Rearrangements,Inversions\n'
        handle = open('./iddo2/' + operon + '.csv', 'w')
        handle.write(header2)
        handle.write('\n'.join([','.join(i) for i in file_output2]))
        handle.close()
    print "iddo_result", iddo_result
    #print result.keys()
    
    
    # current format of the dict that is stored here:
    #iddo_result[operon].update({nc:[duplicates, deletions, splits, total_rearrangements, inversions]}) 
    pickle.dump(iddo_result, open("event_dict.p", "wb"))
    pickle.dump(common_result, open("common_event_dict.p", "wb"))


# this function will take a pickled dict of all the information that we are looking at for the project, adn compute per operon
# an all vs. all distmat of observed events.  (it is (total org1 + total org2) / 2)  so really shitty measure but will look to refine
# in a later incarnation
def all_vs_all_distmat(fname, outfolder):
    # the format of this dict currently is: {nc:[duplicates, deletions, splits, total_rearrangements, inversions]}
    operon_result_dict = pickle.load(open(fname, "rb"))
    #print operon_result_dict
    
    
    for operon in sorted(operon_result_dict.keys()):
        print 'operon', operon
        org_list = sorted(operon_result_dict[operon].keys())
        print "org_list", org_list
        #print len(org_list), org_list
        header = 'Organism,' + ','.join(org_list)
        print "header", header
        lines = [header]
        result = {}
        #print "Length_combinations", itertools.combinations_with_replacement(org_list, 2)
        for pair in itertools.combinations_with_replacement(org_list, 2):
            print "Pair", pair
            org1, org2 = pair
            if org1 == org2:
                print "Same", org1, org2
                result.update({org1:{org2: 0}})
            else:
                org1_sum = sum([int(i) for i in operon_result_dict[operon][org1]])
                org2_sum = sum([int(i) for i in operon_result_dict[operon][org2]])
                cost = (float(float(org1_sum) + float(org2_sum))/2) # yes, too much float, deal with it
                result.update({org1:{org2: cost}})
                print "Stored information", org1, org2, cost, result[org1][org2]
                #result.update({org2:{org1: cost}})
                #print "Different", org1, org2
                #print "cost", cost, result[org1][org2], result[org2][org1]
        print 'result', len(result), result   
        for org1 in org_list:
            line = [org1]
            for org2 in sorted(result[org1].keys()):
                print "rorg2
                if org1 == org2:
                    line.append('0')
                else:
                    try:
                        line.append(str(result[org1][org2]))
                    except:
                        print 'org1', org1, 'org2', org2
            print "line", line
            lines.append(line)
            print 'lines', lines
        outfile = outfolder + operon + '.csv'
        handle = open(outfile, 'w')
        handle.write('\n'.join(lines))
        handle.close()
                
        
            
    


                
def main():
    ##find_complete_genome_genbank_files('/home/dave/Desktop/all_genbank/', './genbank_genome_path.txt')
    #marker_list = [i.strip().lower() for i in open('marker_genes.txt')]
    #error = make_marker_fasta('NC_000913.gbk', marker_list, 'distance/marker_fasta.fa')
    #if error:
    #    print "Error!"
        
    ##org_paths = return_full_path('proteobacteria_list.txt', 'genbank_genome_path.txt')
    ##org_paths = return_full_path('base_filename_file_possibly_short.txt', 'genbank_genome_path.txt')
    
    
    ##make_target_fasta(marker_list, 'proteobacteria_list.txt', 'genbank_genome_path.txt')
    #create_distmat('rpob_for_msa.fa')
    #print "got past create_distmat"
    
    distmat_file = 'Dist_test.txt'
    
    ##read_distmat('rpob_for_msa.distmat', 'proteobacteria_list.txt', distmat_file)
    
    
    # restore if function below fails
    result_folder = 'results_3/'
    #parse_operon_result_files('optimized_results_proteobacteria/', distmat_file, result_folder, './operon_name_and_genes.txt')
    #print "got past parse_operon_result_files"
    
    
    #parse_operon_result_files2('optimized_results_proteobacteria/', distmat_file, './iddo/', './operon_name_and_genes.txt')
    parse_operon_result_files2('optimized_results_proteobacteria/', distmat_file, result_folder, './operon_name_and_genes.txt')
    
    
    calculate_spearman(result_folder, 'summary_01.txt', .01)
    
    print "got past calculate_spearman"
    
    #create_organism_summary(result_folder)
    create_organism_summary('./iddo2/', 'tally2.csv', [5,6,7,9,12])
    
    # operon_name_and_genes.txt
    #calculate_cost(result_folder, 'operon_name_and_genes.txt')
    
    #calculate_cost(result_folder, 'operon_filtered_results.txt')
    calculate_cost('./iddo2/', 'operon_filtered_results.txt', ['Splits', 'Duplicates', 'Deletions', 'Total_Fusions', 'Inversions'])
        
    all_vs_all_distmat("common_event_dict.p", './operon_dist/')


if __name__ == '__main__':
    main()
