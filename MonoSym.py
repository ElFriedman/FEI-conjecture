import SymFracs
import math
import sys
import heapq


# w is the index of the last 1 in the n+1 bit positions from 1 to n+1
# aka the threshhold
# if w is 0, all -1's. if w is n+1, all 1's
def find_top_k_Cs(n, k=1):
    c_table = []
    for w in range(0,n+2):
        bitmap = []
        for j in range(0, w):
            bitmap.append(1)
        for j in range(0, n+1-w):
            bitmap.append(-1)
        sym = SymFracs.symmetric(bitmap, n)
        assert sym.weight() == 1
        c_table.append((sym.C(), sym, w))
    print(f"The top {k} bitmaps are as follows:")
    top_k = heapq.nlargest(k, c_table, key = lambda tup: tup[0])
    for (i, (c,s,w)) in enumerate(top_k):
        print(f"Rank {i}: threshold {w} leads to C = {c}")
        s.sym_report_abrv()
        print("")


def rank(n, w):
    return min(w, n+1-w)


def concave(n):
    c_table = []
    for w in range(1,n+1): # ignore constant functions with w = 0, w = n+1
        bitmap = []
        for j in range(0, w):
            bitmap.append(1)
        for j in range(0, n+1-w):
            bitmap.append(-1)
        sym = SymFracs.symmetric(bitmap, n)
        assert sym.weight() == 1
        c_table.append((sym.C(), sym, w))
    top_k = heapq.nlargest(n, c_table, key = lambda tup: tup[0])
    last_rank = -1
    for (i, (c,s,w)) in enumerate(top_k):
        curr_rank = rank(n, w)
        if curr_rank < last_rank:
            return False
        else:
            last_rank = curr_rank
    return True


def concave_until(max_n):
    for i in range(1, max_n+1):
        if not concave(i):
            print(f"not concave at {i} bits")
            return
    print(f"concave up through {max_n} bits")


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    # print("best C with n <= " + str_n)
    # find_best_C_until(int(str_n))
    print(f"Studying monotone symmetric boolean functions on up to {max_n} variables")
    concave_until(max_n)
    find_top_k_Cs(67, 70)
