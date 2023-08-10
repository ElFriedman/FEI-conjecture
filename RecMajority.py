import Monomial
from Monomial import Monomial
import Polynomial
from Polynomial import Polynomial
import fractions


def new_root(s, n):
    n1 = Monomial(mono_k_str(s, 1, n), fractions.Fraction(1, 2), n)
    n2 = Monomial(mono_k_str(s+1, 1, n), fractions.Fraction(1, 2), n)
    n3 = Monomial(mono_k_str(s+2, 1, n), fractions.Fraction(1, 2), n)
    n4 = Monomial(mono_k_str(s, 3, n), fractions.Fraction(-1, 2), n)
    return Polynomial.monos_to_poly([n1, n2, n3, n4])


def mono_k_str(s, k, n):
    return ("0" * s) + ("1" * k) + ("0" * (n-s-k))


def poly_maj3(poly1, poly2, poly3):
    term1 = Polynomial.add_polynomials([poly1, poly2, poly3])
    term1.scale(1/2)
    term2 = Polynomial.mult_polynomials([poly1, poly2, poly3])
    term2.scale(-1/2)
    return Polynomial.add_polynomials([term1, term2])


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
