"""
Comparison between the formation energies of structures within:
interval pairs that were generated around basepair of the global structure;
and intervals pairs that were generated around randomly selected positions;
"""
import RNA
import generate_intervals as gi
import subprocess
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def run_RNAduplex(intervals, seq):
    energies = []
    structures = []

    for interval in intervals:
        seq1 = seq[interval[0][0]:interval[0][1]+1]
        seq2 = seq[interval[1][0]:interval[1][1]+1]
        RNAduplex_result = subprocess.run(["RNAduplex"], input=f'{seq1}\n{seq2}', encoding='ascii', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = RNAduplex_result.stdout.strip()

        structure, interval1, _, interval2, energy = output.replace("( ", "(").split()
        struct1, struct2 = structure.split('&')
        pos1_str1, pos2_str1 = map(int, interval1.split(','))
        pos1_str2, pos2_str2 = map(int, interval2.split(','))

        # Adjust struct1 to match the length of seq1
        struct1 = '.' * (pos1_str1 - 1) + struct1
        if len(struct1) < len(seq1):
            struct1 += '.' * (len(seq1) - len(struct1))

        # Adjust struct2 to match the length of seq2
        struct2 = '.' * (pos1_str2 - 1) + struct2
        if len(struct2) < len(seq2):
            struct2 += '.' * (len(seq2) - len(struct2))

        energy = float(energy.strip('()'))

        structures.append((struct1, struct2))
        energies.append(energy)

    return energies, structures


def run_RNAcofold(intervals, seq):
    energies = []
    structures = []

    for interval in intervals:
        seq1 = seq[interval[0][0]:interval[0][1]+1]
        seq2 = seq[interval[1][0]:interval[1][1]+1]
        RNAduplex_result = subprocess.run(["RNAcofold"], input=f'{seq1}&{seq2}', encoding='ascii', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = RNAduplex_result.stdout.strip()

        _, structure, energy = output.replace("( ", "(").replace("( ", "(").split()
        struct1, struct2 = structure.split('&')
        energy = float(energy.strip('()'))
        structures.append((struct1, struct2))
        energies.append(energy)

    return energies, structures


seq_len = 3000
interval_lengths = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]
nr_intervals = 300
max_attempts = 8000
seq = RNA.random_string(seq_len, 'ACGU')

(RNAfold_structure, _) = RNA.fold(seq)

base_pairs = gi.parse_dot_bracket(RNAfold_structure)
all_possible_pairs = [(i, j) for i in range(1, seq_len-5) for j in range(i+4, seq_len-1)]

energies_cofold = []
energies_duplex = []
rand_duplex = []
rand_cofold = []

structures_bp_cofold = []
structures_rand_cofold = []
structures_bp_duplex = []
structures_rand_duplex = []
intervals = []
for interval_length in interval_lengths:

    _, selected_bp_intervals, attempt_count = gi.pick_N_basepairs_wo_overlapping_intervals(base_pairs, interval_length, interval_length, interval_length, interval_length, seq_len, nr_intervals, max_attempts)
    nr_intervals = len(selected_bp_intervals)

    _, selected_rand_intervals, attempt_count = gi.pick_N_basepairs_wo_overlapping_intervals(all_possible_pairs, interval_length, interval_length, interval_length, interval_length, seq_len, nr_intervals, max_attempts)

    energy_bp_cofold, structure_bp_cofold = run_RNAcofold(selected_bp_intervals, seq)
    energy_rand_cofold, structure_rand_cofold = run_RNAcofold(selected_rand_intervals, seq)

    energy_bp_duplex, structure_bp_duplex = run_RNAduplex(selected_bp_intervals, seq)
    energy_rand_duplex, structure_rand_duplex = run_RNAduplex(selected_rand_intervals, seq)

    energies_cofold += (energy_bp_cofold)
    rand_cofold += ["around basepairs"]*nr_intervals

    energies_cofold += (energy_rand_cofold)
    rand_cofold += ["random positions"]*nr_intervals

    energies_duplex += (energy_bp_duplex)
    rand_duplex += ["around basepairs"]*nr_intervals

    energies_duplex += (energy_rand_duplex)
    rand_duplex += ["random positions"]*nr_intervals

    intervals += [interval_length*2+1]*nr_intervals*2

    structures_bp_cofold += (structure_bp_cofold)
    structures_rand_cofold += (structure_rand_cofold)
    structures_bp_duplex += (structure_bp_duplex)
    structures_rand_duplex += (structure_rand_duplex)

    print(f"Interval size {interval_length*2+1} done")


data = pd.DataFrame({
    'Interval Length': intervals,
    'RNACofold Free Energy [kcal/mol]': energies_cofold,
    'RNADuplex Free Energy [kcal/mol]': energies_duplex,
    'Interval Creation RNACofold': rand_cofold,
    'Interval Creation RNADuplex': rand_duplex,

})

data.to_csv("baseline.csv")
data['RNADuplex Free Energy [kcal/mol]'] = data['RNADuplex Free Energy [kcal/mol]'].replace(100000, 0)

sns.set(font_scale=2.3)
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(11.7, 8.27))  # initializes figure and plots

# p = sns.lineplot(data, errorbar="ci", x='Interval Length', y='RNACofold Free Energy [kcal/mol]', hue='Interval Creation RNACofold', legend="brief")
p = sns.boxplot(data, x='Interval Length', y='RNACofold Free Energy [kcal/mol]', hue='Interval Creation RNACofold', legend="brief")
p.legend_.set_title(None)
plt.savefig("RNAcofold_baseline.pdf")
plt.show()

sns.set(font_scale=2.3)
sns.set_style("whitegrid")
fig, ax1 = plt.subplots(figsize=(11.7, 8.27))  # initializes figure and plots

# p = sns.lineplot(data, errorbar="ci", x='Interval Length', y='RNADuplex Free Energy [kcal/mol]', hue='Interval Creation RNADuplex', legend="brief")
p = sns.boxplot(data, x='Interval Length', y='RNADuplex Free Energy [kcal/mol]', hue='Interval Creation RNADuplex', legend="brief")
p.legend_.set_title(None)
plt.savefig("RNAduplex_baseline.pdf")
plt.show()
