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
	'-n', '--num_variants', default=1000, type=int,
	help='Number of generated variants'
)
parser.add_argument(
	'-me', '--mut_exp', default=70, type=int,
	help='Percent chance of mutating a codon in experimental set'
)
parser.add_argument(
	'-mb', '--mut_bg', default=10, type=int,
	help='Percent chance of mutating a codon in background set'
)
args = parser.parse_args()

# Original sequences
seqDNA = 'GTTATTTTACAACTTGAACGTGTT' # Starting protein sequence
seq5Prime = 'AAAGGCAGT' # 5' flanking sequence
seq3Prime = 'GGTGGAAGT' # 3' flanking sequence

print(f'Generating {args.num_variants} variants\n'
      f'Starting sequence: {seqDNA}\n'
      f'Full sequence: {seq5Prime}-{seqDNA}-{seq3Prime}\n'
      f'Mutation Odds:\n'
      f'* Experimental: {args.mut_exp} %\n'
      f'* Background: {args.mut_bg} %\n')

setInit = False
dir = 'data'
if not os.path.exists(dir):
    os.makedirs(dir, exist_ok=True)
if setInit:
    pathExp = os.path.join(dir, 'variantsExp.fastq')
    pathBg = os.path.join(dir, 'variantsBg.fasta')
else:
    seqDNA = 'CCTTATATTCAGATTGATAATGCG'
    pathExp = os.path.join(dir, 'variantsExp2.fastq')
    pathBg = os.path.join(dir, 'variantsBg2.fasta')


def generateVariants(sequence, mutationOdds=4, numVariants=50):
    bases = ['A', 'T', 'C', 'G']
    variants = [('variant_0', f'{seq5Prime}{sequence}{seq3Prime}')]
    numVariants = 100 - numVariants

    while len(variants) < numVariants:
        var = list(sequence)
        for index in range(0, len(sequence), 3):
            codon = sequence[index:index + 3]
            if  random.randint(0, 100) <= mutationOdds:
                bp = codon[random.randint(0, 2)]
                bpNew = random.choice([b for b in bases if b != bp])
                var[index] = bpNew
        if var != sequence:
            name = f"variant_{len(variants)}"
            subCassette = seq5Prime + "".join(var) + seq3Prime
            variants.append((name, subCassette))

    return variants


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
variantsExp = generateVariants(
    seqDNA, mutationOdds=args.mut_exp, numVariants=args.num_variants
)
variantsBg = generateVariants(
    seqDNA, mutationOdds=args.mut_bg, numVariants=args.num_variants
)

# Save to FASTA and FASTQ
saveSeqs(variantsExp,  pathExp)
saveSeqs(variantsBg,  pathBg)

print(f'Saved {args.num_variants} variants at:\n'
      f'     {pathExp}\n'
      f'     {pathBg}')
