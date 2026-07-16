"""
This file implements the section 3.8.2 of Falcon's documentation.
"""

from OWL.hawk.ntrugen.fft import add_fft, adj_fft, div_fft, fft, ifft, mul_fft

q = 1


def karatsuba(a, b, n):
    """
    Karatsuba multiplication between polynomials.
    The coefficients may be either integer or real.
    """
    if n == 1:
        return [a[0] * b[0], 0]
    n2 = n // 2
    a0 = a[:n2]
    a1 = a[n2:]
    b0 = b[:n2]
    b1 = b[n2:]
    ax = [a0[i] + a1[i] for i in range(n2)]
    bx = [b0[i] + b1[i] for i in range(n2)]
    a0b0 = karatsuba(a0, b0, n2)
    a1b1 = karatsuba(a1, b1, n2)
    axbx = karatsuba(ax, bx, n2)
    for i in range(n):
        axbx[i] -= a0b0[i] + a1b1[i]
    ab = [0] * (2 * n)
    for i in range(n):
        ab[i] += a0b0[i]
        ab[i + n] += a1b1[i]
        ab[i + n2] += axbx[i]
    return ab


def karamul(a, b):
    """
    Karatsuba multiplication, followed by reduction mod (x ** n + 1).
    """
    n = len(a)
    ab = karatsuba(a, b, n)
    return [ab[i] - ab[i + n] for i in range(n)]


def galois_conjugate(a):
    """
    Galois conjugate of an element a in Q[x] / (x ** n + 1).
    Here, the Galois conjugate of a(x) is simply a(-x).
    """
    n = len(a)
    return [((-1) ** i) * a[i] for i in range(n)]


def field_norm(a):
    """
    Project an element a of Q[x] / (x ** n + 1) onto Q[x] / (x ** (n // 2) + 1).
    Only works if n is a power-of-two.
    """
    n2 = len(a) // 2
    ae = [a[2 * i] for i in range(n2)]
    ao = [a[2 * i + 1] for i in range(n2)]
    ae_squared = karamul(ae, ae)
    ao_squared = karamul(ao, ao)
    res = ae_squared[:]
    for i in range(n2 - 1):
        res[i + 1] -= ao_squared[i]
    res[0] += ao_squared[n2 - 1]
    return res


def lift(a):
    """
    Lift an element a of Q[x] / (x ** (n // 2) + 1) up to Q[x] / (x ** n + 1).
    The lift of a(x) is simply a(x ** 2) seen as an element of Q[x] / (x ** n + 1).
    """
    n = len(a)
    res = [0] * (2 * n)
    for i in range(n):
        res[2 * i] = a[i]
    return res


def bitsize(a):
    """
    Compute the bitsize of an element of Z (not counting the sign).
    The bitsize is rounded to the next multiple of 8.
    This makes the function slightly imprecise, but faster to compute.
    """
    val = abs(a)
    res = 0
    while val:
        res += 8
        val >>= 8
    return res


def reduce(f, g, upper_f, upper_g):
    """
    Reduce (F, G) relatively to (f, g).

    This is done via Babai's reduction.
    (F, G) <-- (F, G) - k * (f, g), where k = round((F f* + G g*) / (f f* + g g*)).
    Corresponds to algorithm 7 (Reduce) of Falcon's documentation.
    """
    n = len(f)
    base_bits = max(
        53, bitsize(min(f)), bitsize(max(f)), bitsize(min(g)), bitsize(max(g))
    )

    f_adjust = [elt >> (base_bits - 53) for elt in f]
    g_adjust = [elt >> (base_bits - 53) for elt in g]
    fa_fft = fft(f_adjust)
    ga_fft = fft(g_adjust)

    while 1:
        # Because we work in finite precision to reduce very large polynomials,
        # we may need to perform the reduction several times.
        current_bits = max(
            53,
            bitsize(min(upper_f)),
            bitsize(max(upper_f)),
            bitsize(min(upper_g)),
            bitsize(max(upper_g)),
        )
        if current_bits < base_bits:
            break

        upper_f_adjust = [elt >> (current_bits - 53) for elt in upper_f]
        upper_g_adjust = [elt >> (current_bits - 53) for elt in upper_g]
        upper_fa_fft = fft(upper_f_adjust)
        upper_ga_fft = fft(upper_g_adjust)

        den_fft = add_fft(
            mul_fft(fa_fft, adj_fft(fa_fft)), mul_fft(ga_fft, adj_fft(ga_fft))
        )
        num_fft = add_fft(
            mul_fft(upper_fa_fft, adj_fft(fa_fft)),
            mul_fft(upper_ga_fft, adj_fft(ga_fft)),
        )
        k_fft = div_fft(num_fft, den_fft)
        k = ifft(k_fft)
        k = [round(elt) for elt in k]
        if all(elt == 0 for elt in k):
            break
        # The two next lines are the costliest operations in ntru_gen
        # (more than 75% of the total cost in dimension n = 1024).
        # There are at least two ways to make them faster:
        # - replace Karatsuba with Toom-Cook
        # - mutualized Karatsuba, see ia.cr/2020/268
        # For simplicity reasons, we didn't implement these optimisations here.
        fk = karamul(f, k)
        gk = karamul(g, k)
        for i in range(n):
            upper_f[i] -= fk[i] << (current_bits - base_bits)
            upper_g[i] -= gk[i] << (current_bits - base_bits)
    return upper_f, upper_g


def xgcd(b, n):
    """
    Compute the extended GCD of two integers b and n.
    Return d, u, v such that d = u * b + v * n, and d is the GCD of b, n.
    """
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return b, x0, y0


def ntru_solve(f, g):
    """
    Solve the NTRU equation for f and g.
    Corresponds to NTRUSolve in Falcon's documentation.
    """
    n = len(f)

    if n == 1:
        f0 = f[0]
        g0 = g[0]
        d, u, v = xgcd(f0, g0)
        if d != 1:
            raise ValueError
        return [-q * v], [q * u]
    fp = field_norm(f)
    gp = field_norm(g)
    upper_fp, upper_gp = ntru_solve(fp, gp)
    upper_f = karamul(lift(upper_fp), galois_conjugate(g))
    upper_g = karamul(lift(upper_gp), galois_conjugate(f))
    upper_f, upper_g = reduce(f, g, upper_f, upper_g)
    return upper_f, upper_g


def is_invertible_2(x):
    return (sum(x) % 2) == 1
