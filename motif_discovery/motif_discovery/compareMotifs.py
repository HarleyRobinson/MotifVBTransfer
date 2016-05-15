#! /usr/bin/env python2.4
#
# Compare motifs in tamo format
#

from   TAMO              import MotifTools
from   TAMO.MotifMetrics import ProbeSet
from   TAMO.Clustering   import MotifCompare
from   TAMO.Clustering   import Kmedoids
import sys
import pickle
import pprint


file_unknown = sys.argv[1]# Unknown
file_tfbs = sys.argv[2]# TF db
motifs_unknown = MotifTools.load(file_unknown) 
motifs_tfbs = MotifTools.load(file_tfbs) 

match_dict = {}
for unknown in motifs_unknown:
  tf_list = []
  for tfbs in motifs_tfbs:
    #print 
    #print "Comparing motifs:"
    #print "    %s  vs  %s" % (unknown.source, tfbs.source)
    #print "    Unknown motif ( %s ) vs TFBS ( %s ) " % (unknown, tfbs)
    #print
    joined_motifs = []
    joined_motifs.append(unknown)
    joined_motifs.append(tfbs)
    print joined_motifs
    Dmat = MotifCompare.computeDmat(joined_motifs)
    cluster = Kmedoids.bestaveKMedoids_cluster(Dmat,kmax=2,min_dist=0.2)
    for i in cluster[1].keys():
      cluster_list = cluster[1][i]
      if len(cluster_list) > 1:
        dist = Dmat[0][1]
        tf_list.append( (dist, tfbs.source, tfbs.oneletter) )
        #print 
        #print "*** Motif match found! ***" 
        #tfbs.giflogo("%s" % tfbs.oneletter)
        #for i in cluster_list:
          #print joined_motifs[i]
  novel_id = unknown.source
  tf_list.sort()
  match_dict[novel_id] = {'query_motif': unknown.oneletter,
                          'subject': tf_list }
      #print "//"
  #print "////"
#print all_motifs
#print cluster
#pprint.pprint(match_dict)


fileout = open("COMPARE_MOTIF.txt", 'w')
fileout.write("NovelID\tNovel_motif\tTFBS\tTFBS_motif\tSimilarity\n")

for key in match_dict.keys():
  if len(match_dict[key]['subject']) > 0:
    for i in range(len(match_dict[key]['subject'])):
      novel_id = key
      novel_motif = match_dict[key]['query_motif']
      tfbs_id = match_dict[key]['subject'][i][1]
      tfbs_motif = match_dict[key]['subject'][i][2]
      similarity = 1.0 - match_dict[key]['subject'][i][0]
      fileout.write("%s\t%s\t%s\t%s\t%s\n" % (
        novel_id,novel_motif,tfbs_id,tfbs_motif,similarity))

fileout.close()

