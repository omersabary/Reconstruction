# A Dynamic Programming based Python program for edit
# distance problem


def edit_distance_ops(str1, str2, format_list=False):
    m = len(str1)
    n = len(str2)
    # Create a table to store results of subproblems
    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]

    # Fill d[][] in bottom up manner
    for i in range(m + 1):
        for j in range(n + 1):

            # If first string is empty, only option is to
            # insert all characters of second string
            if i == 0:
                dp[i][j] = j  # Min. operations = j

            # If second string is empty, only option is to
            # remove all characters of second string
            elif j == 0:
                dp[i][j] = i  # Min. operations = i

            # If last characters are same, ignore last char
            # and recur for remaining string
            elif str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]

            # If last character are different, consider all
            # possibilities and find minimum
            else:
                dp[i][j] = 1 + min(dp[i][j - 1],  # Insert
                                   dp[i - 1][j],  # Remove
                                   dp[i - 1][j - 1])  # Replace

    # get edit operations
    ops = []  # ops contains edit operations (insert, delete, replace) which converts s1 to s2
    i = m
    j = n
    while True:
        if i == 0 or j == 0:
            if i == 0:
                for c in str2[:j]:
                    if format_list:
                        ops.append('I' + c)
                    else:
                        ops.append("Insert " + c + " in string1")
            elif j == 0:
                for c in str1[:i]:
                    if format_list:
                        ops.append('D' + c)
                    else:
                        ops.append("Delete " + c + " in string1")

            break
        if str1[i - 1] == str2[j - 1]:
            if format_list:
                ops.append('H' + str1[i - 1])
            else:
                ops.append("Halt " + str1[i - 1] + " in string1")
            i = i - 1
            j = j - 1
        elif dp[i][j] == dp[i - 1][j - 1] + 1:  # replace
            if format_list:
                ops.append('R' + str1[i - 1] + str2[j - 1])
            else:
                ops.append("Replace " + str1[i - 1] + " in string1 to " + str2[j - 1] + " from string2")
            i = i - 1
            j = j - 1
        elif dp[i][j] == dp[i - 1][j] + 1:  # deletion
            if format_list:
                ops.append('D' + str1[i - 1])
            else:
                ops.append("Delete in string1 " + str1[i - 1])
            i = i - 1
        elif dp[i][j] == dp[i][j - 1] + 1:  # deletion
            if format_list:
                ops.append('I' + str2[j - 1])
            else:
                ops.append("Insert in string1 " + str2[j - 1])
            j = j - 1
        else:
            raise Exception("WTF")
    ops.reverse()  # FIXME push to start to improve complexity
    return dp[m][n], ops


if __name__ == '__main__':
    s1 = "TTCGCA"
    s2 = "-C-C-A"

    l, o = edit_distance_ops(s1, s2, format_list=True)
    print("d is   : ", l)
    print("ops are: ", o)
