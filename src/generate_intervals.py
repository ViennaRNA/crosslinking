import random
import argparse
import csv

def parse_dot_bracket(structure):
    """
    parses the dot-bracket structure into single basepairs and prints unmatched parenthesis
    returns: list of bps in the ss
    """
    stack = []
    base_pairs = []
    
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                u = stack.pop()
                v = i
                base_pairs.append((u, v))
            else:
                print(f"Unmatched closing parenthesis at position {i}")
    if stack:
        print("Unmatched opening parentheses at positions:", stack)
    
    return base_pairs

def define_intervals(u, v, d1, d2, d3, d4, seq_len):
    start_u = max(0, u - d1)
    end_u = min(seq_len - 1, u + d2)
    start_v = max(0, v - d3)
    end_v = min(seq_len - 1, v + d4)
    
    return (start_u, end_u), (start_v, end_v)

def intervals_overlap(interval1, interval2):
    """
    checks if intervals overlap
    returns: boolean
    """
    start1, end1 = interval1
    start2, end2 = interval2
    return max(start1, start2) <= min(end1, end2)

def pick_N_basepairs_wo_overlapping_intervals(base_pairs, d1, d2, d3, d4, seq_len, N, max_attempts):
    """
    picks N amount of unique random bps out of the ss (u, v)
    interval length defined by (u-d1, u+d2), (v-d3, v+d4)
    intervals are not allowed to overlap!
    u-d1 truncated to 0 and v+d4 truncated to len(ss) if over limit -> intervals not always the same length!
    only has max_attemps to find bps, to mimick seq depth??
    returns: 
    list of selected bps [(u,v),(u1,v2),...]
    list of selected intervals [((u-d1,u+d2),(v-d3,v+d4)),...]
    number of tries to find the bp, dependent on different unique bp positions and d1-4 interval length
    """
    selected_bps = []
    selected_intervals = []
    tried_basepairs = set()
    attempt_count = 0

    while len(selected_bps) < N and attempt_count < max_attempts:
        selected_bp = random.choice(base_pairs)
        u, v = selected_bp

        attempt_count += 1

        if selected_bp in tried_basepairs:
            continue
        tried_basepairs.add(selected_bp)

        interval_u, interval_v = define_intervals(u, v, d1, d2, d3, d4, seq_len)

        if intervals_overlap(interval_u, interval_v):
            continue

        selected_bps.append(selected_bp)
        selected_intervals.append((interval_u, interval_v))

    return selected_bps, selected_intervals, attempt_count

def main(dot_bracket_structure, d1, d2, d3, d4, N, max_attempts, output_file):

    # Step 1: Parse the dot-bracket notation to get base pairs
    base_pairs = parse_dot_bracket(dot_bracket_structure)
    seq_len = len(dot_bracket_structure)

    # Step 2: Pick N bps without overlapping intervals
    selected_bps, selected_intervals, attempt_count = pick_N_basepairs_wo_overlapping_intervals(
        base_pairs, d1, d2, d3, d4, seq_len, N, max_attempts)

    # Step 3: Sort the selected base pairs and intervals by u in ascending order
    sorted_data = sorted(zip(selected_bps, selected_intervals), key=lambda x: x[0][0])

    # Step 4: Output the selected bps and intervals
    if sorted_data:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(["Base Pair (u, v)", "Interval around u", "Interval around v"])
            
            for bp, intervals in sorted_data:
                u, v = bp
                interval_u, interval_v = intervals
                writer.writerow([f"({u}, {v})", f"({interval_u[0]}, {interval_u[1]})", f"({interval_v[0]}, {interval_v[1]})"])
                
        print(f"Selected {len(sorted_data)} base pairs after {attempt_count} attempts.")
        print(f"Results written to {output_file}")
    else:
        print("No base pairs found where the intervals do not overlap.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
            This script parses a dot-bracket structure, selects base pairs without overlapping intervals around each
            nucleotide, and saves the selected pairs and intervals to a TSV file.
        """
    )

    parser.add_argument(
        'dot_bracket_structure',
        type=str,
        help='The dot-bracket notation of the secondary structure (e.g., ".(((...)))..")'
    )
    parser.add_argument(
        '-N', '--num_base_pairs', 
        type=int, 
        default=10, 
        help='Number of base pairs to pick without overlapping intervals'
    )
    parser.add_argument(
        '-d1', type=int, default=5, 
        help='Interval distance to the left of nucleotide u'
    )
    parser.add_argument(
        '-d2', type=int, default=5, 
        help='Interval distance to the right of nucleotide u'
    )
    parser.add_argument(
        '-d3', type=int, default=5, 
        help='Interval distance to the left of nucleotide v'
    )
    parser.add_argument(
        '-d4', type=int, default=5, 
        help='Interval distance to the right of nucleotide v'
    )
    parser.add_argument(
        '-m', '--max_attempts', 
        type=int, 
        default=1000, 
        help='Maximum number of attempts to find non-overlapping base pairs (proxy for sequencing depth)'
    )
    parser.add_argument(
        '-o', '--output_file', 
        type=str, 
        default="output.tsv", 
        help='Output TSV file for selected base pairs and intervals'
    )

    args = parser.parse_args()

    main(
        dot_bracket_structure=args.dot_bracket_structure,
        d1=args.d1,
        d2=args.d2,
        d3=args.d3,
        d4=args.d4,
        N=args.num_base_pairs,
        max_attempts=args.max_attempts,
        output_file=args.output_file
    )
