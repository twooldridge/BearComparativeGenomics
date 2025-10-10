import sys
import os
from Bio import SeqIO

# usage: python combined_fastas.py [input directory with fastas] [outdir]

indir = sys.argv[1]
outprefix = sys.argv[2]

collection = {}
file_number = 1

for file in os.listdir(indir):
    if file.endswith(".fa"):
        for record in SeqIO.parse(file, "fasta"):
            print('Parsing file ' + str(file_number))
            if file_number == 1:
                collect_id = record.id.split('_')[1]
                species = record.id.split('_')[0]
                record.id = species
                record.description = ''
                collection[collect_id] = [record]
            elif file_number > 1:
                collect_id = record.id.split('_')[1]
                species = record.id.split('_')[0]
                record.id = species
                record.description = ''
                if collect_id in collection:
                    collection[collect_id].append(record)
        file_number += 1

print('Writing matching sequences from ' + str(file_number-1) +' input fastas in ' + outprefix)

for item in collection:
    filename = outprefix + '/' + item + '.unaligned.8spec.fa'
    with open(filename,'w') as out:
        SeqIO.write(collection[item],out,'fasta')
print ('Finished!')