# Symmetric Polynomials with fractions

import math
import fractions
import sys
import heapq


class Symnomial:
    def __init__(self, w_arr, n):
        self.w_arr = w_arr
        self.n = n

    def weight(self):
        weight = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            weight += math.comb(self.n, k) * p
        return weight

    def entropy(self):
        ent = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            if p > 0:
                ent -= math.comb(self.n, k) * p * math.log2(p)
        return ent

    def influence(self):
        inf = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            inf += math.comb(self.n, k) * p * k
        return inf

    def C(self):
        return self.entropy() / self.influence() if self.influence() > 0 else 0

    def __str__(self):
        out = ""
        for k in range(0, self.n + 1):
            out += " " + str(k) + ": " + str(self.w_arr[k]) + "\n"
        return out

    def __repr__(self):
        return self.__str__()

    def sym_report_full(self):
        print("Fourier weight by level:")
        print(self)
        print("Total weight: " + str(self.weight()))
        print("Entropy: " + str(self.entropy()))
        print("Influence: " + str(self.influence()))
        print("C: " + str(self.C()))

    def sym_report_abrv(self):
        print("Total weight: " + str(self.weight()))
        print("Entropy: " + str(self.entropy()))
        print("Influence: " + str(self.influence()))
        print("C: " + str(self.C()))


def kravchuk(m, x, n):
    """
    return evaluation of Kravchuk polynomial $K_m(x; n)$ for integer x.
    """
    return sum([(-1)**j * math.comb(x, j) * math.comb(n-x, m-j) for j in range(m+1)])


def symmetric(bitmap, n):  # len must be n + 1
    coef_arr = []
    for s in range(0, n + 1):  # compute Fourier coeff at subsets S of size s
        numerator = sum([kravchuk(m, s, n) * bitmap[m] for m in range(n+1)])
        coef_arr.append(fractions.Fraction(int(numerator), 2 ** n))
    return Symnomial(coef_arr, n)


def monotone_symmetric(thresh, n):
    coef_arr = []
    # compute the level-0 Fourier coeff
    numerator = sum([math.comb(n, m) for m in range(thresh)])
    coef_arr.append(fractions.Fraction(int(numerator), 2 ** (n - 1)) - 1)
    # compute the rest of the Fourier coeffs, using Kravchuk sum identity
    for s in range(1, n + 1):
        numerator = kravchuk(thresh - 1, s - 1, n - 1)
        coef_arr.append(fractions.Fraction(int(numerator), 2 ** (n - 1)))
    return Symnomial(coef_arr, n)


def show_bitmap(bs):
    return " ".join(["-" if b < 0 else "+" for b in bs])


def find_best_C_until(max_n):
    best_C = 0
    best_sym = 0
    best_bitmap = []
    for n in range(1, max_n + 1):
        for i in range(2 ** n):
            x_str = ('{0:b}'.format(i).rjust((n + 1), '0'))   # Use binary repr of i as "bitmap"
            bitmap = [-1 if b == '1' else 1 for b in x_str]   # Convert from "0/1" to "±1" notation
            sym = symmetric(bitmap, n)
            assert sym.weight() == 1
            if sym.C() > best_C:
                best_C = sym.C()
                best_sym = sym
                best_bitmap = bitmap
    print("Best bitmap is", best_bitmap)
    best_sym.sym_report_full()


def find_top_k_Cs(n, k=1):
    c_table = []
    for i in range(2 ** n):   # need only iterate up to here, due to symmetry
        x_str = ('{0:b}'.format(i).rjust((n + 1), '0'))   # Use binary repr of i as "bitmap"
        bitmap = [-1 if b == '1' else 1 for b in x_str]   # Convert from "0/1" to "±1" notation
        sym = symmetric(bitmap, n)
        assert sym.weight() == 1
        c_table.append((sym.C(), sym, bitmap))
    print(f"The top {k} bitmaps are as follows:")
    top_k = heapq.nlargest(k, c_table, key = lambda tup: tup[0])
    for (i, (c,s,b)) in enumerate(top_k):
        print(f"Rank {i}: bitmap {show_bitmap(b)} leads to C = {c}")
        s.sym_report_abrv()
        print("")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    # print("best C with n <= " + str_n)
    # find_best_C_until(int(str_n))
    print(f"Studying symmetric boolean functions on {n} variables")
    find_top_k_Cs(n, k)
