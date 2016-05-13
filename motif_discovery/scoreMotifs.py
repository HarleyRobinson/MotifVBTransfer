#! /usr/bin/env python2.4

# Author: Alexandre S. Cristino
# email:  alexscristino@gmail.com
# Date:   16 june 2008
# Script: Scoring motifs in upstream regions of specific databases


import os,sys,string
from   TAMO              import MotifTools
from   TAMO.seq          import Fasta
from   TAMO.MotifMetrics import ProbeSet

promoters = ProbeSet(sys.argv[1])
geneset_ids = open(sys.argv[2]).read().split('\n')[:-1]
match_ids = []
prom_ids = promoters.probes.keys()
for id in geneset_ids:
  if id in prom_ids:
    match_ids.append(id)

motifs = MotifTools.load(sys.argv[3])
church = 0.05
rocauc = 0.1
pvalue = 0.05

print "Name\tMotif\tChurch\tRoc-auc\tP-value"
for m in motifs:
  m.church   = promoters.church  (m, match_ids)
#  m.ROC_auc  = promoters.ROC_AUC (m, match_ids)
  m.pvalue   = promoters.p_value (m, match_ids)
  if m.church <= church and m.pvalue <= pvalue:
    print "%s\t%s\t%s\t%s" %\
    (m.source, m, m.church, m.pvalue) 


