import random
from math import pow
from math import tan
from math import pi


class Point:
    def __init__(self, i, value):
        self.id = i
        self.value = value
        self.tan_dict = {}
        self.err_dict = {}

    def add_tan(self, fun_index, tan_val):
        self.tan_dict[fun_index] = tan_val

    def add_err(self, fun_index, err):
        self.err_dict[fun_index] = err


points_array = []


def machine_precision():
    m = 0
    u = 1.0
    while 1 + u != 1:
        m += 1
        u = pow(10, -m)
    return u * 10


def associativity(u):
    x = 1.0
    y = float(u / 10)
    z = float(u / 10)

    if (x + y) + z == x + (y + z):
        print("+c associative")
    else:
        print("+c not associative")

    a = u/100000000
    b = u/100000000
    c = u/4.0
    if (a * b) * c == a * (b * c):
        print("*c associative")
    else:
        print("*c not associative")


def err(aprox_val, val):
    return abs(aprox_val - tan(val))


def T(i, a):
    match i:
        case 1:
            return a
        case 2:
            return (3 * a) / (3 - pow(a, 2))
        case 3:
            return (15 * a - pow(a, 3)) / (15 - 6 * pow(a, 2))
        case 4:
            return (105 * a - 10 * pow(a, 3)) / (105 - 45 * pow(a, 2) + pow(a, 4))
        case 5:
            return (945 * a - 105 * pow(a, 3) + pow(a, 5)) / (945 - 420 * pow(a, 2) + 15 * pow(a, 4))
        case 6:
            return (10395 * a - 1260 * pow(a, 3) + 21 * pow(a, 5)) / (
                    10395 - 4725 * pow(a, 2) + 210 * pow(a, 4) - pow(a, 6))
        case 7:
            return (135135 * a - 17325 * pow(a, 3) + 378 * pow(a, 5) - pow(a, 7)) / (
                    135135 - 62370 * pow(a, 2) + 3150 * pow(a, 4) - 28 * pow(a, 6))
        case 8:
            return (2027025 * a - 270270 * pow(a, 3) + 6930 * pow(a, 5) - 36 * pow(a, 7)) / (
                    2027025 - 945945 * pow(a, 2) + 51975 * pow(a, 4) - 630 * pow(a, 6) + pow(a, 8))
        case 9:
            return (34459425 * a - 4729725 * pow(a, 3) + 135135 * pow(a, 5) - 990 * pow(a, 7) + pow(a, 9)) / (
                    34459425 - 16216200 * pow(a, 2) + 945945 * pow(a, 4) - 13860 * pow(a, 6) + 45 * pow(a, 8))


if __name__ == "__main__":
    u = machine_precision()
    print(u)
    associativity(u)
    for i in range(0, 10001):
        value = random.uniform(-pi / 2, pi / 2)
        p = Point(i, value)
        points_array += [p]
        for j in range(1, 10):
            aprox_tan = T(j, value)
            p.add_tan(j, aprox_tan)
            p.add_err(j, err(aprox_tan, value))
        p.err_dict = dict(sorted(p.err_dict.items(), key=lambda item: item[1]))

    top_3 = dict.fromkeys(range(1, 10), 0)
    for point in points_array:
        keys = iter(point.err_dict)
        for j in range(0, 3):
            idx = next(keys)
            top_3[idx] += 1

    top_3 = sorted(top_3.items(), key=lambda item: item[1], reverse=True)
    print("Most accurate functions, in the following order:")
    for i in range(0,9):
        print(f"T{top_3[i][0]}")
