#!/usr/bin/python3

import argparse
import os
import random

"""
    Generate variants from a template DNA sequence
    This starting sequence will be flanked by separate sequences on the 5' and 3' sides
    
    The sequences generated can be used for trial datapoints COMET 
"""


parser = argparse.ArgumentParser()
parser.add_argument(
	'-n', '--num_variants', default=2000, type=int,
	help='Number of generated variants'
)
parser.add_argument(
	'-me', '--mut_exp', default=60, type=int,
	help='Percent chance of mutating a codon in experimental set'
)
parser.add_argument(
	'-mb', '--mut_bg', default=100, type=int,
	help='Percent chance of mutating a codon in background set'
)
args = parser.parse_args()


# Original sequences
seqDNA = [
    'CCGGAGCATAAACAGTACCGCGTT', 'GCAATCGTTGGGAATGTGAAGTTA',
    'AGCTGCACCGTGAACCTCGAGCCG', 'CAGGGCGCCTTTATTCGCCATGAA',
    'TGGAACCGCTGCATGAAGGCTTAC', 'CACGCTCGCACCTCACAGTTCATG',
    'TTCCACAAGTACCGCGCTCCCCGC', 'GTTGACCAGCACTGGTGGGTTTGG',
    'ATGTGGATGCACCAGATTTGGAAC', 'ATTATGTGGCACGCTGAAATGCGT'
] # Starting protein sequence
seq5Prime = 'AAAGGCAGT' # 5' flanking sequence
seq3Prime = 'GGTGGAAGT' # 3' flanking sequence

print(f'Generating {args.num_variants:,} variants\n\n'
      # f'Starting sequence: {seqDNA}\n'
      f'Full sequence: {seq5Prime}-{seqDNA[0]}-{seq3Prime}\n'
      f'Flanking sequences:\n'
      f'* 5\' {seq5Prime}\n'
      f'* 3\' {seq3Prime}\n\n'
      
      f'Mutation Odds:\n'
      f'* Experimental: {args.mut_exp} %\n'
      f'* Background: {args.mut_bg} %\n')

dir = 'dset/test/'
if not os.path.exists(dir):
    os.makedirs(dir, exist_ok=True)
pathExp = os.path.join(dir, 'variantsExp.fastq')
pathBg = os.path.join(dir, 'variantsBg.fasta')


def generateVariants(sequences, mutationOdds, numVariants, path):
    bases = ['A', 'T', 'C', 'G']
    stop = ['TAG', 'TAA', 'TGA']
    variants = [('variant_0', f'{seq5Prime}{sequences[0]}{seq3Prime}')]
    mutationOdds = mutationOdds

    while len(variants) < numVariants:
        sequence = sequences[random.randint(0, len(sequences) - 1)]
        var = list(sequence)
        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]
            if  random.randint(0, 100) <= mutationOdds:
                bp = codon[random.randint(0, 2)]
                bpNew = random.choice([b for b in bases if b != bp])
                var[index] = bpNew
        var = ''.join(var)
        if var != sequence:
            keepSub = True
            for i in range(0, len(var), 3):
                codon = var[i:i + 3]
                if codon in stop:
                    keepSub = False
                    break
            if keepSub: # Dont save seq with stop codons
                name = f"variant_{len(variants)}"
                subCassette = seq5Prime + var + seq3Prime
                variants.append((name, subCassette))

    # Save to FASTA and FASTQ
    saveSeqs(variants,  path)


def saveSeqs(variants, fileName):
    if '.fastq' in fileName:
        with open(fileName, 'w') as fastq:
            for name, seq in variants:
                quality = ''.join(chr(random.randint(53, 73)) for _ in seq)
                fastq.write(f"@{name}\n{seq}\n+\n{quality}\n")
    else:
        with open(fileName, 'w') as fasta:
            for name, seq in variants:
                fasta.write(f">{name}\n{seq}\n")


# Generate variants
generateVariants(
    seqDNA, mutationOdds=args.mut_exp, numVariants=args.num_variants, path=pathExp
)
generateVariants(
    seqDNA, mutationOdds=args.mut_bg, numVariants=args.num_variants, path=pathBg
)

print(f'The variants were saved at:\n'
      f'     {pathExp}\n'
      f'     {pathBg}')
