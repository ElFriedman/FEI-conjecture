# to add:
# evaluate at a value
# AND, OR, etc. of functions


import math
import fractions
from collections import Counter


class Monomial:
    def __init__(self, x_str, coef, n):
        self.x_str = x_str
        self.size = self.size()
        self.coef = coef
        self.p = self.coef ** 2
        self.n = n

    def size(self):
        size = 0
        for b in self.x_str:
            size += int(b)
        return size

    def scale(self, x):
        self.coef *= x
        self.p = self.coef ** 2

    def entropy(self):
        return -self.p * math.log2(self.p)

    def influence(self):
        return self.p * self.size

    def mult(self, other):
        new_str = ""
        for i in range(0, self.n):
            new_str += "0" if self.x_str[i] == other.x_str[i] else "1"
        new_coef = self.coef * other.coef
        return Monomial(new_str, new_coef, self.n)

    def add(self, other):
        return Monomial(self.x_str, self.coef + other.coef, self.n)

    def __str__(self):
        return self.x_str + "; size " + str(self.size) + "; coef " + str(self.coef)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return Monomial(self.x_str, self.coef, self.n)

    def eval(self, bits):
        val = self.coef
        for i in range(0, self.n):
            if self.x_str[i] == '1':
                val *= bits[i]
        return val




def mult_monomials(mono_arr):
    n = mono_arr[0].n
    new_str = ""
    for i in range(0, n):
        b = 0
        for mono in mono_arr:
            b += int(mono.x_str[i])
        b %= 2
        new_str += str(b)
    new_coef = 1
    for mono in mono_arr:
        new_coef *= mono.coef
    return Monomial(new_str, new_coef, n)


class Polynomial:
    def __init__(self, mono_dict):
        self.mono_dict = mono_dict

    def scale(self, x):
        for key in self.mono_dict:
            self.mono_dict[key].scale(x)

    def weight(self):
        total = 0
        for key in self.mono_dict:
            total += self.mono_dict[key].p
        return total

    def entropy(self):
        ent = 0
        for key in self.mono_dict:
            ent += self.mono_dict[key].entropy()
        return ent

    def influence(self):
        inf = 0
        for key in self.mono_dict:
            inf += self.mono_dict[key].influence()
        return inf

    def C(self):
        return self.entropy()/self.influence() if self.influence() > 0 else 0

    def mult(self, other):
        new_dict = {}
        for key1 in self.mono_dict:
            for key2 in other.mono_dict:
                n1 = self.mono_dict[key1]
                n2 = other.mono_dict[key2]
                prod = n1.mult(n2)
                key = prod.x_str
                if key in new_dict:
                    new_dict[key] = new_dict[key].add(prod)
                else:
                    new_dict[key] = prod
        return Polynomial(new_dict)

    def levels(self):
        levels = {}
        for key in self.mono_dict:
            mono = self.mono_dict[key]
            if mono.size in levels:
                levels[mono.size] += mono.p
            else:
                levels[mono.size] = mono.p
        return levels

    def __str__(self):
        out = ""
        for key in self.mono_dict:
            out += self.mono_dict[key].__str__() + "\n"
        return out

    def __repr__(self):
        return self.__str__()

    def poly_report_full(self):
        print("Fourier dist:")
        print(self)
        print("Total weight: " + str(self.weight()))
        print("Entropy: " + str(self.entropy()))
        print("Influence: " + str(self.influence()))
        print("C: " + str(self.C()))
        print("Distribution by levels: ")
        print(self.levels())

    def poly_report_abrv(self):
        print("Total weight: " + str(self.weight()))
        print("Entropy: " + str(self.entropy()))
        print("Influence: " + str(self.influence()))
        print("C: " + str(self.C()))

    def eval(self, bits):
        out = 0
        for key in self.mono_dict:
            out += self.mono_dict[key].eval(bits)
        return out


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


def mult_polynomials(arr):
    poly = arr.pop(0)
    while len(arr) > 0:
        poly = poly.mult(arr.pop(0))
    return poly


def add_polynomials(arr):
    new_dict = {}
    for poly in arr:
        for key in poly.mono_dict:
            if key in new_dict:
                new_dict[key] = new_dict[key].add(poly.mono_dict[key])
            else:
                new_dict[key] = poly.mono_dict[key].copy()
    return Polynomial(new_dict)


def monos_to_poly(arr):
    new_dict = {}
    for mono in arr:
        key = mono.x_str
        if key in new_dict:
            new_dict[key] = new_dict[key].add(mono)
        else:
            new_dict[key] = mono
    return Polynomial(new_dict)


def new_root(s, n):
    n1 = Monomial(mono_k_str(s, 1, n), fractions.Fraction(1, 2), n)
    n2 = Monomial(mono_k_str(s+1, 1, n), fractions.Fraction(1, 2), n)
    n3 = Monomial(mono_k_str(s+2, 1, n), fractions.Fraction(1, 2), n)
    n4 = Monomial(mono_k_str(s, 3, n), fractions.Fraction(-1, 2), n)
    return monos_to_poly([n1, n2, n3, n4])


def mono_k_str(s, k, n):
    return ("0" * s) + ("1" * k) + ("0" * (n-s-k))


def poly_maj3(poly1, poly2, poly3):
    term1 = add_polynomials([poly1, poly2, poly3])
    term1.scale(1/2)
    term2 = mult_polynomials([poly1, poly2, poly3])
    term2.scale(-1/2)
    return add_polynomials([term1, term2])


def basic_maj_test():
    maj3_1 = new_root(0, 3)
    print("\n---Maj3_1---")
    maj3_1.poly_report_full()

    r1 = new_root(0, 9)
    r2 = new_root(3, 9)
    r3 = new_root(6, 9)
    maj3_2 = poly_maj3(r1, r2, r3)
    print("\n---Maj3_2---")
    maj3_2.poly_report_full()


def maj3_k(k):
    for d in range(1, k):
        print("\n~~~depth of " + str(d))
        polys = []
        for i in range(0, 3**(d-1)):
            polys.append(new_root(3*i, 3**d))
        while len(polys) > 1:
            new_polys = []
            print("Processing " + str(len(polys)) + " polynomials")
            for j in range(0, int(len(polys)/3)):
                new_polys.append(poly_maj3(polys[3*j], polys[3*j+1], polys[3*j+2]))
            polys = new_polys
        polys[0].poly_report_abrv()


def symmetric(bitmap, n):  # len must be n + 1
    coef_arr = []
    for s in range(0, n+1):  # size of subset S
        coef = 0
        for k in range(0, n+1):  # level within sum
            sum_k = 0
            for j in range(0, s+1):  # count of negative terms chosen
                sum_k += ((-1) ** j) * math.comb(n-k, s-j) * math.comb(k, j)
            coef += (math.comb(n, k)/math.comb(n, s)) * sum_k * bitmap[k]
        coef_arr.append(fractions.Fraction(int(coef), 2**n))
    return Symnomial(coef_arr, n)


def find_best_C(max_n):
    best_C = 0
    best_sym = 0
    best_bitmap = []
    for n in range(1, max_n+1):
        for i in range(2 ** (n+1)):
            x_str = ('{0:b}'.format(i).rjust((n+1), '0'))
            bitmap = []
            for b in x_str:
                if b == '1':
                    bitmap.append(-1)
                else:
                    bitmap.append(1)
            sym = symmetric(bitmap, n)
            if sym.weight() != 1:
                print("Error: weight not equal to 1")
            if sym.C() > best_C:
                best_C = sym.C()
                best_sym = sym
                best_bitmap = bitmap
    print(best_bitmap)
    best_sym.sym_report_full()

print("best C n <= 10")
# find_best_C(10)

n = 10
print("\n\nC of AND with n = " + str(n))
bmap = []
bmap.append(1)
for i in range(n):
    bmap.append(-1)
# bmap = [-1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1]
sym = symmetric(bmap, n)
sym.sym_report_full()


