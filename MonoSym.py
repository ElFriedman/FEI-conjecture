import SymFracs
import math
import sys
import heapq
import matplotlib.pyplot as plt


def make_C_table(n, shortcut=False):
    """
    Study monotone symmetric Boolean functions on n variables. Each such function
    is characterized by a single "threshold parameter" w, the Hamming weight at 
    which the function changes from 1 (false) to -1 (true). Note: 0 <= w <= n+1.
    Return a table of 3-tuples (c,s,w), where w is the threshold, s is the 
    corresponding Symnomial object, and c is the entropy/influence ratio.
    """
    c_table = []
    for w in range(0, (n+2)//2):   # we can ignore the rest, by symmetry
        if shortcut:
            sym = SymFracs.monotone_symmetric(w, n)
        else:
            bitmap = [1] * w + [-1] * (n+1-w)
            sym = SymFracs.symmetric(bitmap, n)
        assert sym.weight() == 1
        c_table.append((sym.C(), sym, w))
    return c_table


def find_top_k_Cs(n, k=1, shortcut=False):
    """
    Among monotone symmetric functions on n variables, find the top k c-values
    and print basic info about the corresponding functions.
    """
    print(f"The top {k} bitmaps are as follows:")
    top_k = heapq.nlargest(k, make_C_table(n, shortcut), key = lambda tup: tup[0])
    for (i, (c,s,w)) in enumerate(top_k):
        print(f"Rank {i}: threshold {w} leads to C = {c}")
        s.sym_report_abrv()
        print("")


def rank(n, w):
    return min(w, n+1-w)


def concave(n):
    c_table = []
    # in the loop below, ignore the constant function (w = 0) and use symmetry 
    # to ignore all w > (n+1)//2: thresholds w and (n+1)-w are equivalent
    for w in range(1,(n+1)//2): 
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


def simple_study():
    """
    Temporarily, not being used.
    """
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print(f"Studying monotone symmetric boolean functions on up to {max_n} variables")
    concave_until(max_n)
    find_top_k_Cs(67, 30)


if __name__ == "__main__":
    """
    If called as "python <program> <n>", prepare a table of (w, C) tuples 
    for monotone symmetric functions on n variables. We might then plot these
    values to how C changes as we increase w, with n being very large.
    """
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    wcs = [(w,c) for (c,s,w) in make_C_table(n, shortcut=True)[1:]]
    print("\n".join([f"{w} {c}" for (w,c) in wcs]))
    ws, cs = zip(*wcs)
    plt.plot(ws, cs)
    plt.xlabel('Threshold')
    plt.ylabel('Ent/Inf ratio')
    plt.show()
