# run meme
/usr/local/meme/bin/meme inputfile -dna -maxsites 20 -mod anr -nostatus -text -minw 4 -maxw 10 > output 

# convert meme to tamo
python meme2tamo.py output.meme

# score motifs
python scoreMotifs.py background_seqs target_ids motif_file

# map sites
python /usr/local/lib/python2.7/dist-packages/TAMO/Sitemap.py -f target_seqs -m motifs_tamo -t 0.7
