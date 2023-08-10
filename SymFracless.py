import math


class Symnomial:
    def __init__(self, w_arr, n):
        self.w_arr = w_arr
        self.n = n

    def weight(self):
        weight = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            weight += math.comb(self.n, k) * p
        return weight/(2**(2*self.n))

    def entropy(self):
        ent = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            n2 = 2 * self.n
            if p > 0:
                ent -= math.comb(self.n, k) * p * (math.log2(p) - n2)
        return ent/(2 ** n2)

    def influence(self):
        inf = 0
        for k in range(0, self.n + 1):
            p = self.w_arr[k] ** 2
            n2 = 2 * self.n
            inf += math.comb(self.n, k) * p * k
        return inf/(2 ** n2)

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


def symmetric(bitmap, n):  # len must be n + 1
    coef_arr = []
    for s in range(0, n+1):  # size of subset S
        coef = 0
        for k in range(0, n+1):  # level within sum
            sum_k = 0
            for j in range(0, k+1):  # count of negative terms chosen
                sum_k += ((-1) ** j) * math.comb(n-s, k-j) * math.comb(s, j)
            coef += sum_k * bitmap[k]
        coef_arr.append(int(coef))
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

# print("best C n <= 10")
# find_best_C(10)

print("correct AND of 20")
farr = [1048574, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 2, -2]
s = Symnomial(farr, 20)
s.sym_report_full()
n = 20
print("\n\nC of AND with n = " + str(n))
bmap = []
for i in range(n):
    bmap.append(1)
bmap.append(-1)
# bmap = [-1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1]

sym = symmetric(bmap, n)
sym.sym_report_full()


