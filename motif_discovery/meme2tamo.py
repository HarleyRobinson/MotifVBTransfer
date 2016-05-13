# Author: Alexandre S. Cristino
# email:  alexscristino@gmail.com
# Date:   july 2015
# Script: convert meme motifs output file (txt) to TAMO format 

import os
import sys
import pickle
import re
from TAMO import MotifTools

filename = sys.argv[1]
motif_list = open(filename).read().split('\nMOTIF')[1:]
tamo_list = []
motif_counter = 1
nsites_pat = re.compile("(w= [0-9]+)")

for motif in motif_list:
  m_info1, m_info2 = motif.split('letter-probability matrix: ')
  m_mat = m_info2.split('--------------------------------------------------------------------------------', 1)[0]
  m_mat_header, m_prob_mat = m_mat.split('\n', 1) 
  nsites = int(nsites_pat.findall(m_mat_header)[0].split('= ')[1])
  count_pos = m_prob_mat.split('\n')[:-1]
  count_mat = []
  site_list = []
  for count in count_pos:
    sites = [float(i) for i in count.split()]
    site_list.append(sites)
    count_dict = {'A': int(sites[0] * nsites),
                  'C': int(sites[1] * nsites),
                  'G': int(sites[2] * nsites),
                  'T': int(sites[3] * nsites)}
    count_mat.append(count_dict)
  m = MotifTools.Motif_from_counts(count_mat)
  m.source = "Motif%s | %s" % (motif_counter, m_mat_header)
  tamo_list.append(m)
  motif_counter += 1  
  
MotifTools.save_motifs(tamo_list, "MEME_motifs_%s.tamo" % filename.split('.')[0])

