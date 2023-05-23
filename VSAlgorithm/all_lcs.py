def lcs(str1, str2):
    a = len(str1)
    b = len(str2)
    string_matrix = [[0 for i in range(b + 1)] for i in range(a + 1)]
    for i in range(1, a + 1):
        for j in range(1, b + 1):
            if i == 0 or j == 0:
                string_matrix[i][j] = 0
            elif str1[i - 1] == str2[j - 1]:
                string_matrix[i][j] = 1 + string_matrix[i - 1][j - 1]
            else:
                string_matrix[i][j] = max(string_matrix[i - 1][j], string_matrix[i][j - 1])
    index = string_matrix[a][b]
    res = [""] * index
    i = a
    j = b
    while i > 0 and j > 0:
        if str1[i - 1] == str2[j - 1]:
            res[index - 1] = str1[i - 1]
            i -= 1
            j -= 1
            index -= 1
        elif string_matrix[i - 1][j] > string_matrix[i][j - 1]:
            i -= 1
        else:
            j -= 1
    return res


def all_lcs(str1, str2):
    a = len(str1)
    b = len(str2)
    string_matrix = [[0 for i in range(b + 1)] for i in range(a + 1)]
    for i in range(1, a + 1):
        for j in range(1, b + 1):
            if i == 0 or j == 0:
                string_matrix[i][j] = 0
            elif str1[i - 1] == str2[j - 1]:
                string_matrix[i][j] = 1 + string_matrix[i - 1][j - 1]
            else:
                string_matrix[i][j] = max(string_matrix[i - 1][j], string_matrix[i][j - 1])
    # for l in string_matrix:
    #     print(l)
    return get_lcs_from_matrix(string_matrix, a, b, str1, str2)
    # index = string_matrix[a][b]
    # res = [""] * index
    # i = a
    # j = b
    # while i > 0 and j > 0:
    #     if str1[i - 1] == str2[j - 1]:
    #         res[index - 1] = str1[i - 1]
    #         i -= 1
    #         j -= 1
    #         index -= 1
    #     elif string_matrix[i - 1][j] > string_matrix[i][j - 1]:
    #         i -= 1
    #     else:
    #         j -= 1
    # return res


def get_lcs_from_matrix(string_matrix, i, j, str1, str2):
    res = set()
    if i > 0 and j > 0:
        if str1[i - 1] == str2[j - 1]:
            for r in get_lcs_from_matrix(string_matrix, i - 1, j - 1, str1, str2):
                s = r
                s += str1[i - 1]
                res.add(s)
        # elif string_matrix[i - 1][j] == string_matrix[i][j - 1]:
        #     return get_lcs_from_matrix(string_matrix, i - 1, j, str1, str2).union(get_lcs_from_matrix(string_matrix, i,
        #                                                                                               j - 1, str1,
        #                                                                                               str2))
        elif string_matrix[i - 1][j] > string_matrix[i][j - 1]:
            return get_lcs_from_matrix(string_matrix, i - 1, j, str1, str2)
        else:
            return get_lcs_from_matrix(string_matrix, i, j - 1, str1, str2)
        return res
    else:
        return {''}


if __name__ == '__main__':
    # ress = all_lcs('abcabcaa', 'acbacba')
    ress = all_lcs(
        'TTCGCCAAGACGGCTTATCACCTCTAAAGAGTGCGTGCACTTCGTCGAAAAGTGTGACGAGATCTGAATATGTTGCGAGAAACAG',
        'CCTACCAGTGGAATCTACTTCTGAGAGTGAGCTTGCACCGTTGGCCCCGAGTTTTGGACGAGACTGAGATATATGCAGAGGAGTGT')
    for r in ress:
        print(r)
