import argparse
import re
import random
import RNA


def read_params(path):
    inFile = open(path, 'r')
    paramsList = []

    for line in inFile:
        lineList = line.split()
        for i, token in enumerate(lineList):
            if re.match(r'^[+-]?[0-9]+$', token):
                lineList[i] = int(token)
        paramsList.append(lineList)

    return paramsList


def append_int_to_paramStr(paramsString, number):
    number_str = str(number)
    for _ in range(6-len(number_str)):
        paramsString += ' '
    paramsString += number_str

    return paramsString


def append_INF_to_paramStr(paramsString):
    for _ in range(3):
        paramsString += ' '
    paramsString += 'INF'

    return paramsString


def append_comment_str_to_paramStr(paramsString, string):
    for _ in range(4):
        paramsString += ' '
    paramsString += string

    return paramsString


def make_parameter_string(paramsList):
    # Up until now, the parameters are stored as a list.
    # But the ViennaRNA package needs them as a specially formatted string.
    # This function generates such a string out of a parameter List.
    paramsString = ''
    for line in paramsList:
        for i, token in enumerate(line):
            if isinstance(token, int):
                paramsString = append_int_to_paramStr(paramsString, token)
            elif token == 'INF':
                paramsString = append_INF_to_paramStr(paramsString)
            elif token[0] == '/':
                paramsString = append_comment_str_to_paramStr(paramsString, token)
            elif i == 0:
                paramsString += token
            else:
                paramsString += f' {token}'
        paramsString += '\n'

    return paramsString


def get_temp(maxTempDiff):
    return random.uniform(37-maxTempDiff, 37+maxTempDiff)


def get_salt(maxSaltDiff):
    return random.uniform(1.021-maxSaltDiff, 1.021+maxSaltDiff)


def get_random_pertubation_factor(maxParamFactor):
    return random.uniform(1-maxParamFactor, 1+maxParamFactor)


def perturbe_params(inParams, maxParamFactor):
    newParams = []

    for line in inParams:
        newParams.append([])
        for token in line:
            if (isinstance(token, int)):
                newParams[-1].append(round(token * get_random_pertubation_factor(maxParamFactor)))
            else:
                newParams[-1].append(token)

    return newParams


def make_2d_symetric(paramsList, start_i, end_i):
    # Assumes the entries in paramsList[start_i:end_i+1] [:] correspond to a quadratic 2d matrix
    # Makes this 2d matrix symetric s.t. a[i][j] == a[j][i]

    row_i = -1
    col_i = -1

    for line in paramsList[start_i:end_i+1]:
        row_i += 1
        col_i = -1
        for token in line:
            if not isinstance(token, int):
                continue
            col_i += 1
            if col_i >= row_i:
                continue
            paramsList[start_i+row_i][col_i] = paramsList[start_i+col_i][row_i]

    return paramsList


def make_4d_symetric(paramsList, start_i, end_i):
    # Assumes the entries in paramsList[start_i:end_i+1] [:] correspond to a 4d matrix of size (7,7,5,5)
    # Makes this 4d matrix symetric s.t. a[i][j][k][l] == a[j][i][l][k]

    index_i = 0
    index_j = 0
    index_k = 0
    index_l = -1

    for list_row, line in enumerate(paramsList[start_i:end_i+1]):
        for list_col, token in enumerate(line):
            if not isinstance(token, int):
                continue

            index_l += 1
            if index_l == 5:
                index_l = 0
                index_k += 1
                if index_k == 5:
                    index_k = 0
                    index_j += 1
                    if index_j == 7:
                        index_j = 0
                        index_i += 1
                        if index_i == 7:
                            print("ERROR")
                            exit()

            """
            fourth_j = list_col
            third_j = list_row % 5
            second_j = math.floor((list_row/5) % 7)
            first_j = math.floor((list_row/35) % 35)
            """

            if index_i > index_j or index_k > index_l:
                paramsList[start_i+35*index_i+5*index_j+index_k][index_l] = paramsList[start_i+35*index_j+5*index_i+index_l][index_k]
                # paramsList[start_i+35*i+5*j+k][l] = paramsList[start_i+35*j+5*i+l][k]

    return paramsList


def make_6d_symetric(paramsList, start_i, end_i):
    # Assumes the entries in paramsList[start_i:end_i+1] [:] correspond to a 6d matrix of size (6,6,4,4,4,4)
    # Makes this 6d matrix symetric s.t. a[i][j][k][l][m][n] == a[j][i][m][n][k][l]

    index_i = 0
    index_j = 0
    index_k = 0
    index_l = 0
    index_m = 0
    index_n = -1

    for list_row, line in enumerate(paramsList[start_i:end_i+1]):
        for list_col, token in enumerate(line):
            if not isinstance(token, int):
                continue
            index_n += 1
            if index_n == 4:
                index_n = 0
                index_m += 1
                if index_m == 4:
                    index_m = 0
                    index_l += 1
                    if index_l == 4:
                        index_l = 0
                        index_k += 1
                        if index_k == 4:
                            index_k = 0
                            index_j += 1
                            if index_j == 6:
                                index_j = 0
                                index_i += 1
                                if index_i == 6:
                                    print("ERROR")
                                    exit()

            if index_i > index_j or index_k > index_m or index_l > index_n:
                paramsList[start_i+384*index_i+64*index_j+16*index_k+4*index_l+index_m][index_n] = paramsList[start_i+384*index_j+64*index_i+16*index_m+4*index_n+index_k][index_l]

    return paramsList


def get_unperturbed_MFE_struc(sequence):

    RNA.params_load_RNA_Turner2004()
    fc = RNA.fold_compound(sequence)
    (structure, mfe) = fc.mfe()

    return (structure, mfe)


def get_ener_of_structure(sequence, structure):

    RNA.params_load_RNA_Turner2004()
    fc = RNA.fold_compound(sequence)
    energy = fc.eval_structure(structure)

    return energy


def structure_fits_requirements(sequence, perturbedStructure, minimalBasePairDistance, minimalEnergyDifference):

    unperturbedStructure, unperturbedMFE = get_unperturbed_MFE_struc(sequence)

    if RNA.bp_distance(perturbedStructure, unperturbedStructure) < minimalBasePairDistance:
        return False

    perturbedStructureEnergy = get_ener_of_structure(sequence, perturbedStructure)

    if abs(perturbedStructureEnergy - unperturbedMFE) < minimalEnergyDifference:
        return False

    return True


def generate_perturbed_structure(maxTempDiff, maxSaltDiff, maxParamFactor, sequence, minBPD, minEnerDist):
    # Perturbes the default parameters and then uses them to generate a perturbed structure. Then checks if
    # the structure fits the requirements:
    #   The basepair distance between the perturbes structure and the mfe structure should be larger than minBPD
    #   The difference between the energy of the perturbed structure and mfe (both using default parameters) larger than minEnerDist
    # If those requirements are not met, calculate a different structure with differently perturbed parameters. Try this 50 times
    # and return the latest structure if the requirements are never fulfilled.

    defaultParams = read_params('params/rna_turner2004.par')
    maxTries = 50

    for try_i in range(maxTries):

        perturbedParams = perturbe_params(defaultParams, maxParamFactor)
        perturbedTemp = get_temp(maxTempDiff)
        perturbedSalt = get_salt(maxSaltDiff)

        # The parameters were perturbed individually, but there is
        # some symetry within the parameters.
        # This symetry is now restored
        perturbedParams = make_2d_symetric(perturbedParams, 3, 9)
        perturbedParams = make_2d_symetric(perturbedParams, 12, 18)
        perturbedParams = make_4d_symetric(perturbedParams, 505, 749)
        perturbedParams = make_4d_symetric(perturbedParams, 752, 996)
        perturbedParams = make_6d_symetric(perturbedParams, 3453, 5757)
        perturbedParams = make_6d_symetric(perturbedParams, 5759, 8062)

        perturbedParamsString = make_parameter_string(perturbedParams)
        RNA.params_load_from_string(perturbedParamsString)

        md = RNA.md(temperature=perturbedTemp, salt=perturbedSalt)
        fc = RNA.fold_compound(sequence, md)
        (perturbedStructure, mfe) = fc.mfe()

        if structure_fits_requirements(sequence, perturbedStructure, minBPD, minEnerDist):
            break

    return (perturbedStructure, (perturbedTemp, perturbedSalt, perturbedParamsString))


def generate_random_seq(lenght):
    return RNA.random_string(lenght, 'CGAU')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
            Generate Synthetic Data and writes a sequence, the MFE of that sequence, perturbed structurs of that sequence and the energy parameters used for the pertubation of file.
            A minimum Basepair distance to the unperturbed MFE as well as a minmum energy difference to the MFE can be specified. A new, perturbed structure is generated until those
            requirement are met. After 50 unsuccessful tries, a structure that doesn't satisfy the requirements is returned.
            Needs two folders:
                params: contains the default energy parameters that will be perturbed (params/runa_turner2004.par). The perturbed parameter files will also be written into this folder.
                data: will be used to write the output to. Each file contains the random sequence and the calcualted strucutes. The structure in the first line is calculated using the default parameters.
            """
    )

    parser.add_argument('-n', '--nrSeq', type=int, default=1, help='Nr of sequences to generate structures for')
    parser.add_argument('-s', '--nrStruc', type=int, default=10, help='Nr of structures to be generated for each sequence')
    parser.add_argument('-t', '--maxTempDiff', type=float, default=10, help='How much the temperature is allowed to be changed from the default setting')
    parser.add_argument('-c', '--maxSaltDiff', type=float, default=0.2, help='How much the salt concentration is allowed to be changed from the default setting')
    parser.add_argument('-p', '--maxParamFactor', type=float, default=0.25, help='The energy parameters will be perturbed by multiplication with a random number in the range [1-p,1+p],')
    parser.add_argument('-l', '--length', type=int, default=100, help='Length of the generated sequences/lengths')
    parser.add_argument('-i', '--sequence', type=str, default='', help='Use only the provided sequence. If no sequence is provided, generates a random one. Overrides the length given by -n')
    parser.add_argument('-o', '--out', type=str, default='test', help='Prefix for the output files')
    parser.add_argument('-b', '--minBPD', type=int, default=10, help='Gives a minimum Basepair distance between the perturbed structure and the MFE. Structures that dont satisfy this requirement will be discarded.')
    parser.add_argument('-e', '--minEnerDist', type=float, default=2, help='Gives a minimum Energy difference between the perturbed structure and the MFE. Structures that dont satisfy this requirement will be discarded.')

    args = parser.parse_args()

    if args.nrSeq < 0:
        print('nrSeq needs to be positive')
        exit()
    elif args.nrStruc < 0:
        print('nrStruc needs to be positive')
        exit()
    elif args.maxTempDiff > 37 or args.maxTempDiff < 0:
        print('maxTempDiff should be in [0,37]')
        exit()
    elif args.maxSaltDiff > 1 or args.maxSaltDiff < 0:
        print('maxTempDiff should be in [0,1]')
        exit()
    elif args.maxParamFactor >= 1 or args.maxParamFactor < 0:
        print('maxTempDiff should be in [0,1[')
        exit()
    elif args.length <= 0:
        print('length should be larger than 0')
        exit()

    outPrefix = args.out

    for seqIndex in range(args.nrSeq):

        if args.sequence == '':
            currSeq = generate_random_seq(args.length)
        else:
            currSeq = args.sequence

        currOutFile = open(f'data/{outPrefix}Seq{seqIndex}.csv', 'w')
        currOutFile.write("seq,stuc,temp,salt,ParamPath\n")

        for strucIndex in range(args.nrStruc):

            pertubedStructure, (perturbedTemp, perturbedSalt, perturbedParamsStr) = generate_perturbed_structure(args.maxTempDiff, args.maxSaltDiff, args.maxParamFactor, currSeq, args.minBPD, args.minEnerDist)

            currParamPath = f'params/{outPrefix}Seq{seqIndex}_struc{strucIndex}.par'
            currParamFile = open(currParamPath, 'w')
            currParamFile.write(perturbedParamsStr)
            currParamFile.close()

            currOutFile.write(f"{currSeq},{pertubedStructure},{perturbedTemp},{perturbedSalt},{currParamPath}\n")

        currOutFile.close()
