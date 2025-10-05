from fractions import Fraction
import math
class CustomError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return '{0} '.format(self.message)
        else:
            return 'MyCustomError has been raised'

# ---------- Вспомогательные функции для работы с полиномами над ℤ ----------

def poly_trim(p):
    """Удаляет ведущие нули, оставляя хотя бы один элемент."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_sum(a, b):
    if len(a) < len(b):
        a, b = b, a
    res = a.copy()
    for i in range(len(b)):
        res[i] += b[i]
    return res

def poly_multiply_naive(a, b):
    """Наивное умножение полиномов над ℤ."""
    if not a or not b:
        return [0]
    res = [0] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        for j in range(len(b)):
            res[i + j] += a[i] * b[j]
    return res

def poly_subtract(a, b):
    """Вычитание полиномов над ℤ: a - b."""
    n = max(len(a), len(b))
    res = [0] * n
    for i in range(len(a)):
        res[i] += a[i]
    for i in range(len(b)):
        res[i] -= b[i]
    return poly_trim(res)

# ---------- NTT (для одного модуля) ----------

def ntt_bit_reverse_copy(arr, n):
    j = 0
    res = [0] * n
    for i in range(n):
        res[i] = arr[j]
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
    return res

def ntt_transform(arr, mod, root, invert=False):
    n = len(arr)
    arr = ntt_bit_reverse_copy(arr, n)
    length = 2
    while length <= n:
        wlen = pow(root, (mod - 1) // length, mod)
        if invert:
            wlen = pow(wlen, mod - 2, mod)
        for i in range(0, n, length):
            w = 1
            for j in range(length // 2):
                u = arr[i + j]
                v = (arr[i + j + length // 2] * w) % mod
                arr[i + j] = (u + v) % mod
                arr[i + j + length // 2] = (u - v + mod) % mod
                w = (w * wlen) % mod
        length <<= 1
    if invert:
        inv_n = pow(n, mod - 2, mod)
        for i in range(n):
            arr[i] = (arr[i] * inv_n) % mod
    return arr

def poly_multiply_ntt(a, b, mod, root):
    if not a or not b:
        return [0]
    n = 1
    while n < len(a) + len(b) - 1:
        n <<= 1
    fa = (a + [0] * (n - len(a)))[:]
    fb = (b + [0] * (n - len(b)))[:]
    fa = ntt_transform(fa, mod, root, False)
    fb = ntt_transform(fb, mod, root, False)
    for i in range(n):
        fa[i] = (fa[i] * fb[i]) % mod
    res = ntt_transform(fa, mod, root, True)
    return poly_trim(res[:len(a) + len(b) - 1])

def PolyMultiply(A, B):
    if len(A) * len(B) < 2000:
        return poly_multiply_naive(A, B)
    else:
        return poly_multiply_ntt(A, B, mod=998244353, root=3)
# ---------- Деление полиномов в поле F_p (монический делитель не требуется) ----------

def poly_divide_mod(A, B, mod, root):
    """Деление A / B в поле F_p. Возвращает (Q, R)."""
    A = [a % mod for a in A]
    B = [b % mod for b in B]
    A = poly_trim(A)
    B = poly_trim(B)
    if B == [0]:
        raise ZeroDivisionError
    if len(A) < len(B):
        return ([0], A)
    # Обратный элемент к старшему коэффициенту B
    lc_b = B[-1]
    if lc_b == 0:
        raise ValueError("Leading coefficient is zero")
    inv_lc_b = pow(lc_b, mod - 2, mod)
    Q = [0] * (len(A) - len(B) + 1)
    R = A[:]
    for i in range(len(A) - len(B), -1, -1):
        if len(R) < len(B):
            break
        coef = (R[-1] * inv_lc_b) % mod
        Q[i] = coef
        # Вычитаем coef * x^i * B из R
        to_sub = [0] * i + [(coef * b) % mod for b in B]
        R = [(R[j] - to_sub[j]) % mod for j in range(len(to_sub))]
        if len(R) > len(to_sub):
            R = R[:len(to_sub)]
        # Удаляем ведущие нули
        while len(R) > 1 and R[-1] == 0:
            R.pop()
    Q = poly_trim(Q)
    R = poly_trim(R)
    return (Q, R)

# ---------- CRT (Китайская теорема об остатках) ----------

def crt(rems, mods):
    """
    Восстанавливает число x по остаткам: x ≡ rems[i] (mod mods[i])
    Все mods[i] — попарно простые.
    Возвращает x mod (prod mods), и произведение модулей.
    """
    x = 0
    M = 1
    for m in mods:
        M *= m
    for r, m in zip(rems, mods):
        Mi = M // m
        inv = pow(Mi, -1, m)
        x = (x + r * Mi * inv) % M
    return x, M

def crt_poly(quotients_mod, mods):
    """Восстанавливает полином над ℤ по его значениям по модулям."""
    if not quotients_mod:
        return [0]
    max_len = max(len(q) for q in quotients_mod)
    # Дополняем все частные до одинаковой длины
    for q in quotients_mod:
        while len(q) < max_len:
            q.append(0)
    result = []
    for i in range(max_len):
        rems = [q[i] for q in quotients_mod]
        val, M = crt(rems, mods)
        # Приводим к симметричному представлению [-M/2, M/2]
        if val > M // 2:
            val -= M
        result.append(val)
    return poly_trim(result)

# ---------- Набор модулей, подходящих для NTT ----------

NTT_MODS = [
    (998244353, 3),
    (1004535809, 3),
    (985661441, 3),
]

# ---------- Деление полиномов над ℤ с использованием NTT + CRT ----------

def PolyDivide(A, B):
    """
    Выполняет евклидово деление полиномов A / B над полем.
    A, B — списки коэффициентов [a0, a1, ..., an], a0 + a1*x + ...
    Возвращает (Q, R) такие что A = B*Q + R и deg(R) < deg(B) или R == [0].
    """
    A = [Fraction(x) for x in A]
    B = [Fraction(x) for x in B]

    # Удаляем ведущие нули (но оставляем хотя бы один элемент)
    def trim(p):
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p

    A = trim(A)
    B = trim(B)

    if B == [0]:
        raise ZeroDivisionError("Division by zero polynomial")

    if len(A) < len(B):
        return ([0], A)

    # Инициализируем частное и остаток
    Q = [0] * (len(A) - len(B) + 1)
    R = A[:]

    # Обратный порядок: начинаем со старших степеней
    while len(R) >= len(B) and not (len(R) == 1 and R[0] == 0):
        deg_diff = len(R) - len(B)
        lc_R = R[-1]
        lc_B = B[-1]
        # Деление в поле: lc_R / lc_B
        coef = lc_R / lc_B
        Q[deg_diff] = coef

        # Вычитаем coef * x^deg_diff * B из R
        for i in range(len(B)):
            R[deg_diff + i] -= coef * B[i]

        # Удаляем ведущие нули в R
        R = trim(R)

    # Удаляем ведущие нули в Q
    Q = trim(Q)

    return (Q, R)

# ---------- Глава 1. Algebraic Preliminaries ----------

def Euclidean(a, b):
    """Euclidean algorithm"""
    while b != [0]:
        q, r = PolyDivide(a, b)
        a = b
        b = r
    return a

def FirstHalfExtendedEuclidean(a, b):
    """Half extended Euclidean algorithm"""
    a1 = [1] 
    b1 = [0]
    while b != [0]:
        q, r = PolyDivide(a, b)
        a = b
        b = r
        r1 = poly_subtract(a1, PolyMultiply(q, b1))
        a1 = b1.copy()
        b1 = r1.copy()
    return (a1, a)

def FirstExtendedEuclidean(a, b):
    """Extended Euclidean algorithm - \"half/full\" diophantine version"""
    s, g = FirstHalfExtendedEuclidean(a, b)
    t, r = PolyDivide(poly_subtract(g, PolyMultiply(s, a)), b)
    return (s, t, g)

def HalfExtendedEuclidean(a, b, c):
    """Half extended Euclidean algorithm - diophantine version"""
    s, t, g = FirstExtendedEuclidean(a, b)
    q, r = PolyDivide(c, g)
    if r != [0]:
        raise CustomError('c isn\'t ideal generated by a and b')
    s = PolyMultiply(q, s)
    t = PolyMultiply(q, t)
    if s != [0] and nu(s) >= nu(b):
        q, r = PolyDivide(s, b)
        s = r
        t = poly_sum(t, PolyMultiply(q, a))
    return (s, t)

def ExtendedEuclidean(a, b, c):
    """Extended Euclidean algorithm - \"half/full\" diophantine version"""
    s,t = HalfExtendedEuclidean(a, b, c)
    t, r = PolyDivide(poly_subtract(c, PolyMultiply(s, a)), b)
    return (s, t)

def FirstPartialFraction(a, d):
    """Partial fraction decomposition"""
    a0, r = PolyDivide(a, math.prod(d))
    if len(d) == 1: return (a0, r)
    a1, t = ExtendedEuclidean(d[1:], d[0], r)
    b0, a = FirstPartialFraction(t, d[1:])
    return (poly_sum(a0, b0), a)

def PartialFraction(a, d, e): #TODO
    pass
# ---------- Вспомогательные функции (для совместимости) ----------

def deg(poly):
    return len(poly) - 1

def lc(poly):
    return poly[-1]

def nu(p):
    if len(p) == 1:
        return [abs(p[0])]
    else:
        return deg(p)
# ---------- Пример использования ----------
if __name__ == "__main__":
    # (x^4 - x^3 -19x^2 -11x + 30) / (x - 1)
    A = [30, -11, -19, -1, 1]
    B = [-1, 1]
    Q, R = PolyDivide(A, B)
    print("Q =", Q)  # Ожидается: [-30, -19, 0, 1] → -30 -19x + 0x^2 + x^3
    print("R =", R)  # Ожидается: [0]
    print(ExtendedEuclidean([9], [13], [-1]))