# Welcome to RNA Crosslinking
`RNA folding algorithms` `Pseudo-energy` `Dynamic programming`  `RNA crosslinking`

*This is the source code to generate the data for the pre-print "[Integrating High-Throughput RNA-RNA Interaction Data into RNA Secondary Structure Prediction](https://hal.science/hal-05032055)".*

<img width="100%" src="https://github.com/user-attachments/assets/eb91f281-15b4-4b14-b310-d3dc83d9dd25" />


## Abstract
*In recent years several methods for detecting RNA-RNA in- teractions have become available that use a combination of crosslinking, ligation, and sequencing of the resulting chimeric reads. In principle, such data also convey information on intramolecular helices. They are, however, not accurate enough to identify basepairs directly. Instead only regions of direct contacts can be inferred. Here we show that such data can be incorporated as pseudo-energies into RNA secondary structure prediction algorithms by assigning a bonus term to all potential pairs between crosslinked intervals. Using simulated data we show that given sufficient coverage such data can push the accuracy of the predicted structure to a basepair-wise MCC of above 90%. Moreover, we observe that the beneficial effect of such interval-wise pseudo-energies is quite robust w.r.t. the length of the interval and the value of the bonus term, but depends strongly on the fraction of the sequence that is covered by significant interaction data.*


## Recquirements
This project requires Python 3.8 or higher. 
Whole project is based on the functions from the `ViennaRNA` package.

In addition, the following external Python packages must be installed `numpy`, `pandas`, `matplotlib` and `seaborn`.

## SRC
### RNAcofold_intervals.py and RNAduplex_intervals.py

|  RNAcofold (`RNAcofold_intervals.py`) |  RNAduplex (`RNAduplex_intervals.py`) |
|--------------------------------------|---------------------------------------|
| Compares the **locally predicted structures** using RNAcofold within interval pairs generated around base pairs in the global structure to the **globally predicted part of the structure** on the entire RNA sequence using RNAfold that falls within the same interval pairs. | Runs the RNAduplex function and generates a graph showing duplexed interval pairs. |
| <img width="500" alt="RNA_cofold" src="https://github.com/user-attachments/assets/7cc4fff9-7584-42dc-999e-37fd9bb70533" /> | <img width="500" alt="RNA_duplex" src="https://github.com/user-attachments/assets/78487732-79c4-4251-9b18-f6ed80c2e9a5" /> |




### comparison_random_intervals.py
Ensures comparison between the formation energies of structures within interval pairs, generated around base pairs of the global structure, and those generated around randomly selected positions.


|          |  Default values                                                                |
|------------------|------------------------------------------------------------------------|
| `seq_len`        | `3000` – Length of the RNA sequence                                    |
| `interval_lengths` | `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]` – List of interval half-lengths to test |
| `nr_intervals`   | `300` – Number of intervals to select                                  |
| `max_attempts`   | `8000` – Maximum number of attempts to find non-overlapping intervals  |
| `seq`            | `RNA.random_string(seq_len, 'ACGU')` – Randomly generated RNA sequence of length 3000 |


### generate_data.py
Generates synthetic data and writes:
- sequence,
- MFE (of the sequence),
- perturbed structures (of the sequence)
- energy parameters used for the pertubation of file

A minimum Basepair distance to the unperturbed MFE as well as a minmum energy difference to the MFE can be specified. A new, perturbed structure is generated until these requirements are met. After 50 unsuccessful attempts, a structure that does not meet the requirements is returned.

Use our Visualisation File [mcc_visualisation.ipynb](visualization/mcc_visualization.ipynb) `.ipynb`  in order to view the script and graphs used in the preprint/paper.


### generate_intervals.py
This script parses a dot-bracket RNA secondary structure, selects base pairs with non-overlapping intervals around each nucleotide, and saves the results to a TSV file.
