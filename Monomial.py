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


def mult_monomial_arr(mono_arr):
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
