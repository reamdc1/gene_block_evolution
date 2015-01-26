#!/usr/bin/python

from collections import defaultdict
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse


# Copyright(C) 2014 Iddo Friedberg and David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment


# all_evop == a dicitonary with all events and operons


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will compute a z-score matrix for a pairwise event distance matrix.")
    
    parser.add_argument("-i", "--infile", dest="infile", default='NONE', metavar="FILE", required = True,
                help="A pickled dictionary containing pairwaise event compairisons.")
                
    return parser.parse_args()

def check_options(parsed_args):

    if os.path.exists(parsed_args.infile):
        infile = parsed_args.infile
    else:
        print "The file %s does not exist." % parsed_args.infile
        sys.exit()  

    return infile

event_types = ['deletions','splits','duplications']
def read_events(inpath):
    all_evop = {}
    # This returns the "all_evop" variable used later
    return cPickle.load(open(inpath))

def reduce_event_matrix(inpath):
    all_evop = {}
    s = 0.
    # This returns the "all_evop" variable used later
    ae = cPickle.load(open(inpath))
    for operon in ae:
        all_evop[operon] = {}
        for species1 in ae[operon]:
            for species2 in ae[operon][species1]:
                if species1[:5] == 'Total' or species2[:5] == 'Total':
                    continue
                elif species1 == species2:
                    continue
                elif species1 < species2:
                    all_evop[operon][(species1, species2)] = ae[operon][species1][species2]
                else:
                    all_evop[operon][(species2, species1)] = ae[operon][species1][species2]
    return all_evop
            
"""            
def all_loop_operons(all_evop):
    bar = open("all_loop_operons.csv","w")
    for event_type in event_types:
        event_sum, s = loop_operons(all_evop, event_type)
        bar.write("%s\t%.d\t%d\t%.3f\n" % (event_type,
                  num, event_sum/s ))
    bar.close()

def loop_operons(all_evop, event_type):
    #
    # Loop over all operons with a specific event type.
    # event_sum: the sum of the count of that event type between all species
    event_sum = 0.
    s = 0.
    for operon in all_evop:
        for species_i, species_j in all_evop[operon]:
            for event_type in event_types:
                print operon, species_i,species_j, event_type
                event_sum += all_evop[operon][(species_i,species_j)][event_type]
                s += 1
    return event_sum, s
"""
import numpy as np

#def all_loop_events(all_evop):
#    expected_events = dict.fromkeys(event_types) # number of expected events per pair
#    foo = open("all_loop_events.csv","w")
#    foo.write("operon\tduplications\tsplits\tdeletions\tnpairs\t")
#    foo.write("duplications\tsplits\tdeletions\n")
#    for operon in all_evop:
#        event_vecs = np.array(loop_events(all_evop, operon))
##        event_sum, npairs = loop_events(all_evop, operon)
#        foo.write("%s\t%d\t%d\t%d\t%d\t" % (operon,  
#                  event_vec['duplications'].sum(),
#                  event_vec['splits'].sum(),
#                  event_vec['deletions'].sum(),
#                  npairs
#                  ))
#        foo.write("%.2f\t%.2f\t%.2f\t" % (  
#                  event_vec['duplications'].std(),
#                  event_vec['splits'].std(),
#                  event_vec['deletions'].std(),
#                  ))
#        foo.write("%.2f\t%.2f\t%.2f\t" % (  
#                  event_sum['duplications'].std(),
#                  event_sum['splits'].std(),
#                  event_sum['deletions'].std(),
#                  ))
##        expected_events[operon] = {'duplications':0, 'splits':0, 'deletions':0}
#        expected_events[operon]['duplications'] = event_sum['duplications']/npairs,
#        expected_events[operon]['splits'] = event_sum['splits']/npairs,
#        expected_events[operon]['deletions'] = event_sum['deletions']/npairs,
#    foo.close()
#    return expected_events

def all_loop_events(all_evop):
    event_zscore = {}
    all_event_zscore = {}
    for operon in all_evop:
        event_zscore = loop_events(all_evop, operon)
        all_event_zscore[operon] = event_zscore
    return all_event_zscore
        

def loop_events(all_evop, operon):
    # Get all events for a single operon
    event_vec = defaultdict(list)
    event_zscore = defaultdict(dict)
    event_mean = {}
    event_std = {}
    s = 0.
    npairs = 0.
    for species_i, species_j in all_evop[operon]:
        npairs += 1
        for event_type in all_evop[operon][(species_i,species_j)]:
            event_vec[event_type].append(all_evop[operon][(species_i,species_j)][event_type])
    for event_type in event_vec:
        event_mean[event_type] = np.mean(event_vec[event_type])
        event_std[event_type] = np.std(event_vec[event_type])
    for species_i, species_j in all_evop[operon]:
        event_zscore[(species_i,species_j)] = dict.fromkeys(event_types) 
        for event_type in all_evop[operon][(species_i,species_j)]:
#            print event_type
#            print "mean, std", event_mean[event_type], event_std[event_type]
            event_zscore[(species_i,species_j)][event_type] = ((all_evop[operon][(species_i,species_j)][event_type] - event_mean[event_type]) 
            / float(event_std[event_type]))

    return event_zscore

            

if __name__ == '__main__':

    parsed_args = parser_code()
    
    infile = check_options(parsed_args)

    #all_evop = reduce_event_matrix(sys.argv[1])
    all_evop = reduce_event_matrix(infile)
    
#    all_loop_operons(all_evop)
    all_event_zscore = all_loop_events(all_evop)
    #cPickle.dump(all_event_zscore,open("%s.events.pik" % sys.argv[1],"w"))
    cPickle.dump(all_event_zscore,open("%s.events.pik" % infile,"w"))

    
