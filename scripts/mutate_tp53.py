from Bio import SeqIO
from Bio.Seq import Seq
record = SeqIO.read("TP(53).gb", "genbank")
for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.extract(record.seq)

dna = cds_seq
pos=40
dna[:pos] + dna[pos+1:]
print(dna[:pos] + dna[pos+1:])

healthy_protein = Seq(dna).translate()
print(healthy_protein)

mutated_protein = Seq(mutated).translate()
print(mutated_protein)

print("*" in mutated_protein)

mutated_protein.find("*")

print(len(healthy_protein))
print(len(mutated_protein))

stop_aa_index=mutated_protein.find("*")
stop_codon_start=stop_aa_index*3

print(dna[stop_codon_start:stop_codon_start+3])
print(mutated[stop_codon_start:stop_codon_start+3])

window_start = stop_codon_start-15
window_end = stop_codon_start+12

print("Healthy window", dna[window_start:window_end])
print("Mutated window", mutated[window_start:window_end])
print(stop_codon_start)
stop=stop_aa_index

BRCT_START=1646
BRCT_END=1859
print("STOP at", stop)

if stop < BRCT_START:
  print("STOP IS BEFORE BRCT-> HARMFUL! DOMAIN WON'T REFORM")
elif stop <= BRCT_START:
  print("STOP IS INSIDE THE BRCT-> HARMFUL!")
else: stop <= BRCT_END
print("STOP IS AFTER BRCT-> LIKELY MILD.")

mutation_type='nonsense'

if mutation_type == "frameshift":
    print("Severity: harmful (frameshift)")

elif mutation_type == "nonsense":
    if stop < BRCT_START:
        print("Severity: harmful (early STOP)")

    elif stop <= BRCT_END:
        print("Severity: harmful (STOP inside domain)")

    else:
        print("Severity: mild/uncertain (late STOP)")

elif mutation_type == "missense":

    if BRCT_START <= stop <= BRCT_END:
        print("Severity: harmful (missense inside domain)")

    else:
        print("Severity: mild/uncertain (missense outside domain)")

elif mutation_type == "silent":
    print("Severity: mild (silent mutation)")

