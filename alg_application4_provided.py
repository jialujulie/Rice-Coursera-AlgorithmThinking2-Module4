"""
Provide code and solution for Application 4
"""

DESKTOP = True

import math
import random
import urllib2

if DESKTOP:
    import matplotlib.pyplot as plt
    import project4 as student
else:
    import simpleplot
    import userXX_XXXXXXX as student

# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"


###############################################
# provided code

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urllib2.urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict


def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urllib2.urlopen(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urllib2.urlopen(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print "Loaded a dictionary with", len(word_list), "words"
    return word_list


# question 1
scoring_matrix = read_scoring_matrix(PAM50_URL)
seq_x = read_protein(HUMAN_EYELESS_URL)
seq_y = read_protein(FRUITFLY_EYELESS_URL)
consensusseq = read_protein(CONSENSUS_PAX_URL)

alignment_matrix = student.compute_alignment_matrix(seq_x,seq_y,scoring_matrix,False)
score, string_Hu, string_Fr = student.compute_local_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix)
print string_Hu

newstring_Hu = ""
for elem in string_Hu:
    if elem != '-':
        newstring_Hu += elem
print newstring_Hu
newstring_Fr = ""
for elem in string_Fr:
    if elem != '-':
        newstring_Fr += elem



alignment_matrix_Hum_local_Con = student.compute_alignment_matrix(newstring_Hu,consensusseq,scoring_matrix,True)
score1, str_Hu_Con, str_Con_Hu = student.compute_global_alignment(newstring_Hu,consensusseq,scoring_matrix,
                                                                  alignment_matrix_Hum_local_Con)

alignment_matrix_Fr_local_Con = student.compute_alignment_matrix(newstring_Fr,consensusseq,scoring_matrix,True)
score2, str_Fr_Con, str_Con_Fr = student.compute_global_alignment(newstring_Fr,consensusseq,scoring_matrix,
                                                                  alignment_matrix_Fr_local_Con)
def cal_percentage(str1, str2):
    count = 0
    num = len(str1)
    for i in range(num):
        if str1[i] == str2[i]:
            count += 1
    return float(count)/num

print cal_percentage(str_Hu_Con,str_Con_Hu), cal_percentage(str_Fr_Con, str_Con_Fr)

def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    scoring_distribution ={}
    while num_trials:
        shuffledlist = list(seq_y)
        random.shuffle(shuffledlist)
        rand_y = "".join(shuffledlist)

        loc_align_matric = student.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        score = max(col for row in loc_align_matric
                    for col in row)

        if score in scoring_distribution.keys():
            scoring_distribution[score] += 1
        else:
            scoring_distribution[score] = 1

        num_trials -= 1
    return scoring_distribution


num_trials = 1000
dict = generate_null_distribution(seq_x,seq_y,scoring_matrix, num_trials)
normalized_dict = [float(value)/num_trials for value in dict.values()]
print normalized_dict

plt.bar(dict.keys(), normalized_dict)
plt.title('score distribution of local alignment')
plt.xlabel('Score')
plt.ylabel('Fraction of total trials')
plt.show()