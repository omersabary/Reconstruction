from collections import Counter


def majority(chars: list):
    s = sum([int(x) for x in chars])
    return '0' if s <= len(chars) / 2 else '1'


def b(cluster: list):
    if len(cluster)==0:
        return ""
    le = len(cluster[0])
    #print(len(cluster))
    #cluster = [c for c in cluster if len(c) > 1 ]
    for i, c in enumerate(cluster):
        #if len(c) == 0:
            #print(i)
        #    continue
        if len(c) != le:
            print(len(c), le)
            print(cluster[0])
            print(cluster[1])
            print(c)
            print("what")
        assert len(c) == le
    res = ""
    for i in range(le):
        # res += majority([c[i] for c in cluster])
        res += mode([c[i] for c in cluster])
    return res


def dH(x, y):  # hamming distance
    # return bin(int(x, 2) ^ int(y, 2)).count('1')
    return sum(ch1 != ch2 for ch1, ch2 in zip(x, y))


def delta_semi_match(y, x, d):
    if len(x) != len(y):
        print(len(x), len(y))
    assert len(x) == len(y)
    return dH(x, y) < d * len(x)


def mode(len_list):
    oc_table = Counter(len_list)
    return oc_table.most_common(1)[0][0]


def alg(m, y, l, delta, r, gamma):
    assert m > 0
    assert m == len(y)
    assert l > 0
    assert 0 <= delta <= 1
    assert r > 0
    # 1
    # indices array, equivalent to pointers array in the original algorithm.
    P = [0 for i in range(m)]
    current_M = [i for i in range(m)]
    next_M = []
    # i = 0  # loop index
    x = ""  # result word
    # 2
    while len(current_M) >= gamma * m:
        # a
        x += b([y[j][P[j]:P[j] + l] for j in current_M])
        # b
        next_M = []

        for j in range(m):
            # c
            if j in current_M:
                if abs(len(y[j]) - P[j] - l) >= l and delta_semi_match(y[j][P[j]:P[j] + l], x[-l:], delta):
                    next_M.append(j)
                P[j] += l
            # d
            elif j not in current_M:
                if P[j] >= len(y[j]):
                    continue
                flag = False
                for k in range(P[j] - r * l, P[j] + r * l + 1):
                    if (k + l) >= len(y[j]) or k < 0:
                        break
                    if delta_semi_match(y[j][k:k + l], x[-l:], delta):
                        flag = True
                        P[j] = k + l
                        if abs(len(y[j]) - k - l) >= l:
                            next_M.append(j)
                        break
                if not flag:
                    P[j] += l
        # e
        current_M = next_M
    # 3
    l_tag = mode([abs(len(y[j]) - P[j]) for j in range(m)])
    if l_tag == 0:
        return x
    current_M = [j for j in range(m) if len(y[j]) - P[j] == l_tag]
    x += b([y[j][P[j]:P[j] + l_tag] for j in current_M])
    return x


def run_alg(cluster, orig_len, p=0.3):
    return alg(len(cluster), cluster, 5, (1. + p) / 2.0, 2, 3 / 4)

# if __name__ == '__main__':
#     # y^1, y^2, ..., y^m
#     cluster = ['AGGA', 'TGGTA', 'GTA']
#     gamma = 3 / 4
#     l: int = 2  # window size
#     delta = (1 + 0) / 2  # delta (for hamming distance) parameter
#     r = 2  # flexibility parameter
#
#     res = alg(len(cluster), cluster, l, delta, r, gamma)
#     print(res)
