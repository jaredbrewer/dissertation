#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os, sys

script_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_path)

def parser(ref, esx_utr):

  fasta = open(str(ref))
  fastas = {}

  for f in SimpleFastaParser(fasta):
    name = ">" + f[0].split("|")[1]
    seq = f[1]
    fastas[name] = seq

# This code block is good for writing out the intermediate, but is not strictly necessary. Overwrites the original, but probably okay.

# os.remove(ref)
# with open(ref, 'w') as cdna:
#   for gene, seq in fastas.items():
#     cdna.write(gene + "\n" + seq + "\n")

### Everything above is very general in parsing these files for any use, the subsequent lines are to find and replace the esxMJKPW sequences with their corresponding 5' UTR sequences for the purposes of the subsequent alignment. ###

  esx = open(str(esx_utr))
  esxs = {}

  for f in SimpleFastaParser(esx):
    name = ">" + f[0]
    seq = f[1]
    esxs[name] = seq

  esx_genes = ["esxM", "esxJ", "esxP", "esxW", "esxK"]

  [fastas.pop(">" + key) for key in esx_genes]
  new = fastas | esxs

  with open('Mtb_H37Rv.fa', 'w') as cdna:
    for gene, seq in new.items():
      cdna.write(gene + "\n" + seq + "\n")

if __name__ == '__main__':
  globals()[sys.argv[1]](sys.argv[2], sys.argv[3])

