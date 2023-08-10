# to add:
# evaluate at a value
# AND, OR, etc. of functions


import fractions
from Monomial import Monomial


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





