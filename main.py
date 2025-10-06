from fractions import Fraction
from typing import List, Tuple, Optional

# ---------- 1.1 Groups, Rings and Fields ----------
# В Python мы работаем с полем Fraction => все операции корректны по умолчанию.

# ---------- 1.2 Euclidean Division and Pseudo-Division ----------

def poly_trim(p: List[Fraction]) -> List[Fraction]:
    """Удаляет ведущие нули, оставляя хотя бы один элемент."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def PolyDivide(A, B):
    """
    Евклидово деление полиномов: A = B * Q + R, deg(R) < deg(B) или R = [0].
    A, B — списки Fraction, B ≠ [0].
    """
    if not B or (len(B) == 1 and B[0] == 0):
        raise ZeroDivisionError("Division by zero polynomial")
    
    # Копируем и удаляем ведущие нули
    A = [Fraction(x) for x in A]
    B = [Fraction(x) for x in B]
    A = poly_trim(A[:])
    B = poly_trim(B[:])

    if len(A) < len(B):
        return ([Fraction(0)], A[:])

    lc_B = B[-1]
    Q = [Fraction(0)] * (len(A) - len(B) + 1)
    R = A[:]

    while len(R) >= len(B) and not (len(R) == 1 and R[0] == 0):
        deg_diff = len(R) - len(B)
        coef = R[-1] / lc_B
        Q[deg_diff] = coef

        # Вычитаем coef * x^deg_diff * B из R
        for i in range(len(B)):
            R[deg_diff + i] -= coef * B[i]

        # Удаляем ведущие нули
        while len(R) > 1 and R[-1] == 0:
            R.pop()

    return (poly_trim(Q), poly_trim(R))
# ---------- 1.3 The Euclidean Algorithm ----------

def Euclidean(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    """Алгоритм Евклида для НОД полиномов."""
    a = poly_trim([Fraction(x) for x in a])
    b = poly_trim([Fraction(x) for x in b])
    while b != [Fraction(0)]:
        _, r = PolyDivide(a, b)
        a, b = b, r
    # Делаем моническим
    lc = a[-1]
    return [coef / lc for coef in a]

def ExtendedEuclidean(a: List[Fraction], b: List[Fraction]) -> Tuple[List[Fraction], List[Fraction], List[Fraction]]:
    """Расширенный алгоритм Евклида: s*a + t*b = gcd(a,b)"""
    a = poly_trim([Fraction(x) for x in a])
    b = poly_trim([Fraction(x) for x in b])
    if b == [Fraction(0)]:
        lc = a[-1]
        return ([Fraction(1)/lc], [Fraction(0)], [coef/lc for coef in a])
    s1, s2 = [Fraction(1)], [Fraction(0)]
    t1, t2 = [Fraction(0)], [Fraction(1)]
    while b != [Fraction(0)]:
        q, r = PolyDivide(a, b)
        a, b = b, r
        s1, s2 = s2, [x - y for x, y in zip(s1, poly_mul(q, s2))]
        t1, t2 = t2, [x - y for x, y in zip(t1, poly_mul(q, t2))]
    lc = a[-1]
    s1 = [coef / lc for coef in s1]
    t1 = [coef / lc for coef in t1]
    g = [coef / lc for coef in a]
    return (s1, t1, g)

def HalfExtendedEuclidean(a: List[Fraction], b: List[Fraction]) -> Tuple[List[Fraction], List[Fraction]]:
    """Полу-расширенный алгоритм Евклида: s*a = gcd(a,b) (mod b)"""
    a = poly_trim([Fraction(x) for x in a])
    b = poly_trim([Fraction(x) for x in b])
    s1, s2 = [Fraction(1)], [Fraction(0)]
    while b != [Fraction(0)]:
        q, r = PolyDivide(a, b)
        a, b = b, r
        s1, s2 = s2, [x - y for x, y in zip(s1, poly_mul(q, s2))]
    return (s1, a)

# ---------- Вспомогательные функции для полиномов ----------

def poly_mul(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    if not a or not b or (len(a) == 1 and a[0] == 0) or (len(b) == 1 and b[0] == 0):
        return [Fraction(0)]
    res = [Fraction(0)] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        for j in range(len(b)):
            res[i + j] += a[i] * b[j]
    return poly_trim(res)

def poly_add(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    n = max(len(a), len(b))
    res = [Fraction(0)] * n
    for i in range(len(a)):
        res[i] += a[i]
    for i in range(len(b)):
        res[i] += b[i]
    return poly_trim(res)

def poly_sub(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    n = max(len(a), len(b))
    res = [Fraction(0)] * n
    for i in range(len(a)):
        res[i] += a[i]
    for i in range(len(b)):
        res[i] -= b[i]
    return poly_trim(res)

# ---------- 1.4 Resultants and Subresultants ----------
# Для краткости и соответствия главе 1 реализуем только необходимое для 1.7

# ---------- 1.5 Polynomial Remainder Sequences ----------
# Не требуется напрямую для 1.7, но используется в Yun

# ---------- 1.6 Primitive Polynomials ----------

def poly_content(p: List[Fraction]) -> Fraction:
    """Содержание полинома — НОД коэффициентов (в ℚ всегда 1, но для целых — нет)."""
    # В ℚ содержание всегда 1, но для совместимости с целыми:
    if not p:
        return Fraction(0)
    from math import gcd
    def lcm(a, b):
        return a * b // gcd(a, b)
    nums = []
    dens = []
    for c in p:
        if c != 0:
            nums.append(c.numerator)
            dens.append(c.denominator)
    if not nums:
        return Fraction(0)
    g = abs(nums[0])
    for x in nums[1:]:
        g = gcd(g, abs(x))
    l = dens[0]
    for x in dens[1:]:
        l = lcm(l, x)
    return Fraction(g, l)

def poly_primitive_part(p: List[Fraction]) -> List[Fraction]:
    """Примитивная часть полинома."""
    cont = poly_content(p)
    if cont == 0:
        return [Fraction(0)]
    return [coef / cont for coef in p]

# ---------- 1.7 Squarefree Factorization ----------

def poly_derivative(p: List[Fraction]) -> List[Fraction]:
    """Производная полинома."""
    if len(p) <= 1:
        return [Fraction(0)]
    return [Fraction(i) * p[i] for i in range(1, len(p))]

def Squarefree(A: List[Fraction]) -> List[List[Fraction]]:
    """
    Алгоритм Юня для squarefree-факторизации.
    Возвращает [A1, A2, ..., Am] такие что A = A1 * A2^2 * ... * Am^m,
    и все Ai — squarefree и попарно взаимно просты.
    """
    if not A or A == [Fraction(0)]:
        return []
    A = poly_trim([Fraction(x) for x in A])
    c = poly_content(A)
    S = poly_primitive_part(A)
    S_prime = poly_derivative(S)
    G = Euclidean(S, S_prime)
    if G == [Fraction(1)]:
        return [S]
    Y = poly_sub(S_prime, poly_mul(poly_derivative(G), [Fraction(len(G) - 1)]))
    Z = poly_sub(poly_mul(S, [Fraction(len(S) - 1)]), poly_mul(G, poly_derivative(S)))
    H = Euclidean(Y, Z)
    A1 = poly_div_exact(S, H)
    factors = [A1]
    k = 1
    while G != [Fraction(1)]:
        k += 1
        Y = poly_sub(S_prime, poly_mul(poly_derivative(G), [Fraction(len(G) - 1)]))
        Z = poly_sub(poly_mul(S, [Fraction(len(S) - 1)]), poly_mul(G, poly_derivative(S)))
        H = Euclidean(Y, Z)
        Ak = poly_div_exact(G, H)
        factors.append(Ak)
        G = H
        S_prime = poly_derivative(G)
    # Умножаем первый множитель на содержание
    factors[0] = poly_mul([c], factors[0])
    return factors

def poly_div_exact(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    """Точное деление a / b (предполагается, что делится без остатка)."""
    q, r = PolyDivide(a, b)
    if r != [Fraction(0)]:
        raise ValueError("Polynomial division is not exact")
    return q

# ---------- Вспомогательные функции из книги ----------

def deg(poly: List[Fraction]) -> int:
    return len(poly) - 1

def lc(poly: List[Fraction]) -> Fraction:
    return poly[-1] if poly else Fraction(0)

# ---------- Тесты (опционально) ----------
if __name__ == "__main__":
    # Пример из книги: A = x^8 + 6x^6 + 12x^4 + 8x^2
    A = [Fraction(0), Fraction(0), Fraction(8), Fraction(12), Fraction(0), Fraction(6), Fraction(0), Fraction(1)]
    factors = Squarefree(A)
    print("Squarefree factors:")
    for i, f in enumerate(factors, 1):
        print(f"A{i} =", f[::-1])
    # Ожидается: A = x * (x^2 + 2)^3 → [x, x^2+2]