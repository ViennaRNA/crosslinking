"""
Script used to compare
the structures predicted locally using RNAcofold within interval pairs generated around basepairs in the global structure
to the part of the structure predicted globally on the complete RNA sequence using RNAfold that falls within the same interval pairs.
"""
import csv
import generate_intervals as GeIn
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import RNA
import seaborn as sns
from subprocess import run, PIPE
import time

# ____________________________Sequence_generation____________________________#


def random_rna(length):
    """
    Generates a random RNA sequence of a given length.

    Args:
        length (int): The length of the RNA sequence to generate.

    Returns:
        str: A random RNA sequence composed of 'A', 'C', 'G', and 'U'.
    """
    return ''.join(random.choice('ACGU') for _ in range(length))

    # ____________________________Ptable_intervals____________________________#


def extract_ptable_intervals(ptable, base1_interval, base2_interval):
    """
    Extracts and returns the elements from the ptable based on the given intervals.

    Args:
        ptable (list): List representing the ptable.
        base1_interval (tuple): Interval for the first set of positions.
        base2_interval (tuple): Interval for the second set of positions.

    Returns:
        tuple: Lists of elements for base1_interval and base2_interval.
    """

    interval1_elements = ptable[max(0, base1_interval[0]+1):min(len(ptable), base1_interval[1] + 2)]
    interval2_elements = ptable[max(0, base2_interval[0]+1):min(len(ptable), base2_interval[1] + 2)]

    return interval1_elements, interval2_elements

    # _____________________________Loose_upaired____________________________#


def remove_unpaired(structure, elements, interval1, interval2):
    """
    Modifies the structure by replacing positions with dots if their paired base is outside the target interval.

    Args:
        structure (str): Original dot-bracket notation structure.
        elements (list): List of paired positions from the ptable (1-based index).
        interval (tuple): Interval for the paired positions to keep.

    Returns:
        str: Modified structure.
    """
    interval2_start, interval2_end = interval2
    interval1_start, interval1_end = interval1
    modified_structure = []

    # Ignore the first element in 'elements' since it only indicates length
    for i, elem in enumerate(elements):
        # Check if paired position exists and is OUTSIDE the interval
        if elem > 0 and not ((interval2_start+1 <= elem) and (elem <= interval2_end+1)) and not ((interval1_start+1 <= elem) and (elem <= interval1_end+1)):
            modified_structure.append('.')
        else:
            modified_structure.append(structure[i])

    return ''.join(modified_structure)

# ____________________________RNA_duplex____________________________#


def run_RNAcofold(seq1, seq2):
    """
    Runs the RNAcofold command to predict RNA-RNA interactions and parses the output.

    Args:
        seq1, seq2 (str): RNA sequences to analyze.

    Returns:
        tuple: Parsed output components (struct1, struct2, pos1_str1, pos2_str1, pos1_str2, pos2_str2, energy).
    """
    RNAcofold_result = run(["RNAcofold"], input=f'{seq1}&{seq2}', encoding='ascii', stdout=PIPE, stderr=PIPE)
    output = RNAcofold_result.stdout.strip()

    # parts = output.split()
    # if len(parts) != 5:
    #     raise ValueError(f"Unexpected RNAduplex output format: {output}")

    # Parse the output
    _, structure, energy = output.replace("( ", "(").replace("( ", "(").split()
    struct1, struct2 = structure.split('&')
    energy = float(energy.strip('()'))

    return struct1, struct2, energy
    # ____________________________Save_results____________________________#


def run_RNAcofold_C(seq1, const1, seq2, const2):
    """
    Runs the RNAduplex command to predict RNA-RNA interactions and parses the output.

    Args:
        seq1, seq2 (str): RNA sequences to analyze.

    Returns:
        tuple: Parsed output components (struct1, struct2, pos1_str1, pos2_str1, pos1_str2, pos2_str2, energy).
    """
    RNAduplex_result = run(["RNAcofold", "-C"], input=f'{seq1}&{seq2}\n{const1}&{const2}', encoding='ascii', stdout=PIPE, stderr=PIPE)
    output = RNAduplex_result.stdout.strip()

    _, structure, energy = output.replace("( ", "(").replace("( ", "(").split()
    struct1, struct2 = structure.split('&')
    energy = float(energy.strip('()'))

    return struct1, struct2, energy


def save_results_to_csv(filename, sorted_data, RNAsequence, original_structure, ptabletest, distances, energies):
    """
    Saves the results to a CSV file.

    Args:
        filename (str): Name of the output CSV file.
        sorted_data (list): Sorted base pairs and intervals.
         (str): RNA .
        original_structure (str): Original dot-bracket notation structure.
        ptabletest (list): Parsed table from RNA structure.
    """
    file_exists = os.path.isfile(filename) and os.path.getsize(filename) > 0

    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        if not file_exists:
            writer.writerow(["Base Pair (u, v)", "Interval around u", "Interval around v", "first_base", "second_base", "full_sequence", "RNAcofold_1st_structure", "RNAcofold_2nd_structure", "RNAcofold", "Distance", "Energy"])

        for i, (selected_bp, intervals) in enumerate(sorted_data):
            u, v = selected_bp
            base1_interval, base2_interval = intervals
            first_base = RNAsequence[base1_interval[0]+1:base1_interval[1]+2]
            second_base = RNAsequence[base2_interval[0]+1:base2_interval[1]+2]
            RNAfold_1st_structure = original_structure[base1_interval[0]+1:base1_interval[1]+2]
            RNAfold_2nd_structure = original_structure[base2_interval[0]+1:base2_interval[1]+2]

            full_sequence = RNAsequence[u+1:v+2]
            RNAcofold_result = run_RNAcofold(first_base, second_base)
            RNAcofold_1st_structure, RNAcofold_2nd_structure, _ = RNAcofold_result
            interval1_elements, interval2_elements = extract_ptable_intervals(ptabletest, base1_interval, base2_interval)

            writer.writerow([
                f"({selected_bp})",
                f"({base1_interval[0]}, {base1_interval[1]}, {list(interval1_elements)})",
                f"({base2_interval[0]}, {base2_interval[1]}, {list(interval2_elements)})",
                f"{first_base}_{second_base}", f"{RNAfold_1st_structure}_{RNAfold_2nd_structure}",
                full_sequence, RNAcofold_1st_structure, RNAcofold_2nd_structure, RNAcofold_result,
                distances[i], energies[i]
            ])

    # _____________________________DATA_PLOTTING_______________________________#


def data_plotting(data):
    sns.set(font_scale=2.3)
    sns.set_style("whitegrid")
    fig, ax1 = plt.subplots(figsize=(11.7, 8.27))
    ax2 = ax1.twinx()

    sns.lineplot(data=data, x='d', y='Mean Distance', label='Mean Dist.', marker="o", ax=ax1, color="blue")
    ax1.fill_between(data['d'],
                     [max(0, m+s) for m, s in zip(data['Mean Distance'], data['StdDev'])],
                     [max(0, m-s) for m, s in zip(data['Mean Distance'], data['StdDev'])],
                     alpha=0.3, color='blue', label='± Std Dev')

    sns.lineplot(data=data, x='d', y='Median Distance', label='Median Dist.', color='red', marker='X', ax=ax1)

    sns.lineplot(data=data, x='d', y='Mean Ener', label=r'$\Delta E$', color='green', marker='s', ax=ax2)
    ax2.fill_between(data['d'],
                     [m+s for m, s in zip(data['Mean Ener'], data['StdDev Ener'])],
                     [m-s for m, s in zip(data['Mean Ener'], data['StdDev Ener'])],
                     alpha=0.3, color='green', label='± Std Dev')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="upper center", fontsize="20")
    ax1.legend(bbox_to_anchor=(100, 100))

    for spine in ax1.spines.values():
        spine.set_edgecolor('black')
    for spine in ax2.spines.values():
        spine.set_edgecolor('black')

    ax1.xaxis.label.set_color('black')
    ax1.yaxis.label.set_color('black')
    ax2.yaxis.label.set_color('black')

    ax1.grid(True, which='both', axis='both', color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
    ax2.grid(False, which='both', axis='x')

    ax1.tick_params(axis='both', which='major', labelsize=20, colors='black')
    ax2.tick_params(axis='both', which='major', labelsize=20, colors='black')

    ax1.set_xlabel('Interval Length [nt]')
    ax1.set_ylabel('Base Pair Distance []')
    ax2.set_ylabel(r'$\Delta$ Energy [kcal/mol]')

    plt.savefig("RNAcofold.pdf")
    plt.show()

    # ____________________________MAIN_FUNCTION_______________________________#


def main(seq_len, d1, d2, d3, d4, N, max_attempts, output_file, output2_file):
    """
    Main function to generate a random RNA sequence, fold it, select base pairs, and save the results.

    Args:
        length (int): Length of the RNA sequence.
        d1, d2, d3, d4 (int): Interval parameters.
        N (int): Number of base pairs to select.
        max_attempts (int): Maximum number of attempts to select base pairs.
        output_file (str): Name of the output CSV file.
    """

    RNAsequence = random_rna(length)                        # sequence comes from the function random_rna in predefined lenght
    RNAfold_structure, mfe = RNA.fold(RNAsequence)          # using the generated sequence RNAfold returns structure and energy for this structure
    ptabletest = RNA.ptable(RNAfold_structure)              # in order to find out "whoparis with who" we use ptable for the structure - BP relations
    base_pairs = GeIn.parse_dot_bracket(RNAfold_structure)  # list of base pair positions - BP localisation

    selected_bps, selected_intervals, attempt_count = GeIn.pick_N_basepairs_wo_overlapping_intervals(base_pairs, d1, d2, d3, d4, seq_len, N, max_attempts)
    sorted_data = sorted(zip(selected_bps, selected_intervals), key=lambda x: x[0][0])

    distances = []
    energies = []

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(["Base Pair (u, v)", "Interval around u", "Interval around v", "Sequence around u", "Sequence around v", "Ptable u", "Ptable v", "Score", "Base Pair Distance", "d1", "d2", "d3", "d4", "Length"])

    if sorted_data:
        for bp, intervals in sorted_data:
            u, v = bp
            base1_interval, base2_interval = intervals
            first_base = RNAsequence[base1_interval[0]:base1_interval[1]+1]
            second_base = RNAsequence[base2_interval[0]:base2_interval[1]+1]
            structure1 = RNAfold_structure[base1_interval[0]:base1_interval[1]+1]
            structure2 = RNAfold_structure[base2_interval[0]:base2_interval[1]+1]

            RNAcofold_result = run_RNAcofold(first_base, second_base)
            RNAcofold_1st_structure, RNAcofold_2nd_structure, RNAcofold_energy = RNAcofold_result

            interval1_elements, interval2_elements = extract_ptable_intervals(ptabletest, base1_interval, base2_interval)

            trimmed_structure1 = remove_unpaired(structure1, interval1_elements, base1_interval, base2_interval)
            trimmed_structure2 = remove_unpaired(structure2, interval2_elements, base2_interval, base1_interval)

            _, _, RNAfold_energy = run_RNAcofold_C(first_base, trimmed_structure1.replace('.', 'x'), second_base, trimmed_structure2.replace('.', 'x'))

            RNAcofold_structure_interval = RNAcofold_1st_structure + RNAcofold_2nd_structure
            RNAfold_structure_interval = trimmed_structure1 + trimmed_structure2
            bp_distance = RNA.bp_distance(RNAcofold_structure_interval, RNAfold_structure_interval)
            distances.append(bp_distance)
            energies.append(RNAcofold_energy-RNAfold_energy)

        save_results_to_csv(output_file, sorted_data, RNAsequence, RNAfold_structure, ptabletest, distances, energies)
    else:
        print("No base pairs found where the intervals do not overlap.")

    with open(output2_file, mode='a', newline='') as file2:
        writer2 = csv.writer(file2, delimiter='|')
        writer2.writerow(["random RNA Sequence", "RNAfold Structure", "Energy", "Time", "Length", "Max Attempts", "N", "d1", "d2", "d3", "d4"])
        runtime = time.time() - start_time
        writer2.writerow([RNAsequence, RNAfold_structure, mfe, runtime, length, max_attempts, N, d1, d2, d3, d4])

    return distances, energies


if __name__ == "__main__":

    length = int(input("Enter the desired length for the random RNA sequence generation: ").strip())
    start_time = time.time()
    default_settings = input("Do you want to use default settings for interval generation? (yes/no): ").strip().lower()

    mean_ener_diff = []
    stddev_ener_diff = []
    distances = []
    energies = []
    means = []
    medians = []
    stddev = []
    num_dist_zero = []

    ds = []
    interval_lengths = []
    N = 0
    max_attempts = 0

    if default_settings == 'yes':
        ds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
        N = 3000
        max_attempts = 8000
        intverval_lengths = [d*2+1 for d in ds]

    else:
        N = int(input("Enter the maximal number of non-overlapping base pairs: ").strip())
        max_attempts = int(input("Enter the maximal number of attempts: ").strip())
        d1 = int(input("Enter parameter d: ").strip())
        step = int(input("Enter the step size: ").strip())
        repetitions = int(input("Enter how many intervals should be generated: ").strip())

        ds = [d1 + step * i for i in range(repetitions)]
        intverval_lengths = [d*2+1 for d in ds]

    for d in ds:
        dist, ener = main(
            seq_len=length,
            d1=d,
            d2=d,
            d3=d,
            d4=d,
            N=N,
            max_attempts=max_attempts,
            output_file="output_RNAcofold.csv",
            output2_file="output2_RANcofold.csv"
        )
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Run time: {elapsed_time:.2f} s")
        mean_ener_diff.append(np.mean(np.array(ener)))
        stddev_ener_diff.append(np.std(np.array(ener)))
        distances.append(dist)
        energies.append(ener)
        means.append(np.mean(np.array(distances[-1])))
        medians.append(np.median(np.array(distances[-1])))
        stddev.append(np.std(np.array(distances[-1]), ddof=1))
        num_dist_zero.append(100*dist.count(0)/len(dist))

    data = pd.DataFrame({
        'd': intverval_lengths,
        'Mean Distance': means,
        'Median Distance': medians,
        'Zero Distance': num_dist_zero,
        'StdDev': stddev,
        'Mean Ener': mean_ener_diff,
        'StdDev Ener': stddev_ener_diff
    })
    data_plotting(data)
