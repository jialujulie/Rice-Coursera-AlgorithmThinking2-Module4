"""project 4 - dynamic programming for sequence alignment"""
#build scoring matrix
def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """scoring matrix"""
    new_alphabet = set(alphabet)
    new_alphabet.add('-')
    matrix = dict()
    #matrix = {'-': {'-': 0}}
    #matrix['-'] = dict()
    #matrix['-']['-'] = off_diag_score

    for char_i in new_alphabet:
        matrix[char_i] = dict()
        for char_j in new_alphabet:
            if char_i == char_j and char_i != '-':
                matrix[char_i][char_j] = diag_score
            elif char_i == '-' or char_j == '-':
                matrix[char_i][char_j] = dash_score
            else:
                matrix[char_i][char_j] = off_diag_score

    return matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """alignment matrix to store all alignment score"""

    index_m = len(seq_x)
    index_n = len(seq_y)

    alighment_matrix = list([] for i in range(index_m+1))
    #print alighment_matrix

    alighment_matrix[0].append(0)

    for index_i in range(index_m):

        score = alighment_matrix[index_i][0] + scoring_matrix[seq_x[index_i]]['-']

        if (not global_flag) and score < 0:

            alighment_matrix[index_i+1].append(0)
        else:
            alighment_matrix[index_i+1].append(score)


        #alighment_matrix[index_i+1].append(alighment_matrix[index_i][0] + scoring_matrix[seq_x[index_i]]['-'])

    for index_j in range(index_n):

        score = alighment_matrix[0][index_j] + scoring_matrix['-'][seq_y[index_j]]

        if (not global_flag) and score < 0:
            alighment_matrix[0].append(0)
        else:
            alighment_matrix[0].append(score)


        #alighment_matrix[0].append(alighment_matrix[0][index_j] + scoring_matrix['-'][seq_y[index_j]])

    for index_i in range(index_m):
        #alighment_matrix[index_i] = list()
        for index_j in range(index_n):
            row_char = seq_x[index_i]
            col_char = seq_y[index_j]


            #score1 = alighment_matrix[index_i][index_j+1] + scoring_matrix[row_char]['-']
            #score2 = alighment_matrix[index_i][index_j] + scoring_matrix[row_char][col_char]
            #score3 = alighment_matrix[index_i+1][index_j] + scoring_matrix['-'][col_char]

            maxscore = max(alighment_matrix[index_i][index_j+1] + scoring_matrix[row_char]['-'],
                           alighment_matrix[index_i][index_j] + scoring_matrix[row_char][col_char],
                           alighment_matrix[index_i + 1][index_j] + scoring_matrix['-'][col_char])

            if (not global_flag) and maxscore < 0:
                alighment_matrix[index_i+1].append(0)
            else:
                alighment_matrix[index_i+1].append(maxscore)


    return alighment_matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """get the global alignment"""
    index_i = len(seq_x)
    index_j = len(seq_y)

    align_x = str()
    align_y = str()

    while index_i and index_j:

        if alignment_matrix[index_i][index_j] == alignment_matrix[index_i-1][index_j-1] + scoring_matrix[seq_x[index_i-1]][seq_y[index_j-1]]:
            align_x = seq_x[index_i-1] + align_x
            align_y = seq_y[index_j-1] + align_y
            index_i -= 1
            index_j -= 1
        elif alignment_matrix[index_i][index_j] == alignment_matrix[index_i-1][index_j] + scoring_matrix[seq_x[index_i-1]]['-']:
            align_x = seq_x[index_i-1] + align_x
            align_y = '-' + align_y
            index_i -= 1
        else:
            align_x = '-' + align_x
            align_y = seq_y[index_j-1] + align_y
            index_j -= 1

    while index_i:
        align_x = seq_x[index_i - 1] + align_x
        align_y = '-' + align_y
        index_i -= 1
    while index_j:
        align_x = '-' + align_x
        align_y = seq_y[index_j - 1] + align_y
        index_j -= 1

    return (alignment_matrix[-1][-1], align_x, align_y)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """get the local alignment string"""

    max_score, max_index = max((x, (i,j))
                               for i, row in enumerate(alignment_matrix)
                               for j, x in enumerate(row))

    index_i = max_index[0]
    index_j = max_index[1]

    align_x = str()
    align_y = str()

    while index_i and index_j and alignment_matrix[index_i][index_j]:

        if alignment_matrix[index_i][index_j] == alignment_matrix[index_i-1][index_j-1] + scoring_matrix[seq_x[index_i-1]][seq_y[index_j-1]]:
            align_x = seq_x[index_i-1] + align_x
            align_y = seq_y[index_j-1] + align_y
            index_i -= 1
            index_j -= 1
        elif alignment_matrix[index_i][index_j] == alignment_matrix[index_i-1][index_j] + scoring_matrix[seq_x[index_i-1]]['-']:
            align_x = seq_x[index_i-1] + align_x
            align_y = '-' + align_y
            index_i -= 1
        else:
            align_x = '-' + align_x
            align_y = seq_y[index_j-1] + align_y
            index_j -= 1

    while index_i and alignment_matrix[index_i][index_j]:
        align_x = seq_x[index_i - 1] + align_x
        align_y = '-' + align_y
        index_i -= 1
    while index_j and alignment_matrix[index_i][index_j]:
        align_x = '-' + align_x
        align_y = seq_y[index_j - 1] + align_y
        index_j -= 1

    return (max_score, align_x, align_y)


#print build_scoring_matrix(set(['A','T','C','G']),6,2,-4)

#print compute_alighnment_matrix('A','A',scoringmatrix,True)