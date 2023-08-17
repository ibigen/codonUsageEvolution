""" Has all values that are constants """
from collections import OrderedDict

class Constants(object):
	ordered_codons = OrderedDict()
	codons_per_aminoacid = {'AAA': 'Lys', 'AAG': 'Lys', 'AAC': 'Asn', 'AAU': 'Asn',
							'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACU': 'Thr',
							'AGA': 'Arg', 'AGG': 'Arg', 'CGA': 'Arg', 'CGC': 'Arg',
							'CGG': 'Arg', 'CGU': 'Arg', 'AGC': 'Ser', 'AGU': 'Ser',
							'UCA': 'Ser', 'UCC': 'Ser', 'UCG': 'Ser', 'UCU': 'Ser',
							'AUA': 'Ile', 'AUC': 'Ile', 'AUU': 'Ile', 'AUG': 'Met',
							'UUC': 'Phe', 'UUU': 'Phe', 'UAC': 'Tyr', 'UAU': 'Tyr',
							'UGC': 'Cys', 'UGU': 'Cys', 'UGG': 'Trp', 'CUA': 'Leu',
							'CUC': 'Leu', 'CUG': 'Leu', 'CUU': 'Leu', 'UUA': 'Leu',
							'UUG': 'Leu', 'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro',
							'CCU': 'Pro', 'CAC': 'His', 'CAU': 'His', 'CAA': 'Gln',
							'CAG': 'Gln', 'GUA': 'Val', 'GUC': 'Val', 'GUG': 'Val',
							'GUU': 'Val', 'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala',
							'GCU': 'Ala', 'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly',
							'GGU': 'Gly', 'GAC': 'Asp', 'GAU': 'Asp', 'GAA': 'Glu',
							'GAG': 'Glu', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'}
	# to order them by amino acid

	TOTAL_CODONS = [str(key).upper().replace('U', 'T') for key in codons_per_aminoacid]

	for key, value in codons_per_aminoacid.items():
		if value not in ordered_codons:
			ordered_codons[value] = [str(key).upper().replace('U', 'T')]
		else:
			ordered_codons[value].append(str(key).upper().replace('U', 'T'))



	GENOME_KEY = "genome"
	GENE_CAI = "GENE CAI"
	GENOME_CAI = "GENOME CAI"



	def __init__(self):
		pass
