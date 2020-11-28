from flint import *

q = 9223372036854775837

coeffn = {
    100: 1668737891791324096,
    1000: 2150352211794749766,
    10000: 5784169777460408978,
    100000: 8571646336776701591,   # 0.8 s
    1000000: 3646050264881874559,   # 8.5 s
    10000000: 6365781791077436955,   # 108 s
    100000000: 4042838613910496565,   # 1260/1348 s, virt/peak/res/peak(MB): 1012.96 21675.76 825.20 15090.02
}


def klooster(a,b,m):
    pi = acb.pi()
    j = acb(1j)
    s = 0
    for x in range(m):
        if fmpz(x).gcd(m) == 1:
            if x == 0 and m == 1:
                y = 0
            else:
                y = int(1 / nmod(x, m))
            s += (2*pi*j*(a*x + b*y)/m).exp()
    s = s.real
    return s

def klooster2(a,b,m):
    if m == 1:
        return arb(1)
    if m == 2:
        return arb(-1)**(a+b)
    pi = arb.pi()
    s = 0
    a %= m
    b %= m
    for x in range(m//2+1):
        if fmpz(x).gcd(m) == 1:
            if x == 0 and m == 1:
                y = 0
            else:
                y = int(1 / nmod(x, m))
            v = fmpq(a*x + b*y, m)
            s += (2*pi*(a*x + b*y)/m).cos()
    s = 2 * s.real
    return s

def klooster3(a,b,m):
    fac = fmpz(m).factor()
    r = arb(1)
    for p, e in fac:
        m1 = p**e
        m2 = m // p**e
        if m2 == 1:
            n1 = 0
        else:
            n1 = int(1 / nmod(m1, m2))
        n2 = int(1 / nmod(m2, m1))
        r *= klooster2((n2*a) % m1, (n2*b) % m1, m1)
        a = (n1*a) % m2
        b = (n1*b) % m2
        m = m2
    return r


for i in range(10000):
    a = randint(-30,30)
    b = randint(-30,30)
    n = randint(1,30)
    S1 = klooster(a,b,n)
    S2 = klooster3(a,b,n)
    print(a,b,n)
    assert abs(S1-S2) < 1e-10

def coeff(n):
    from time import clock
    ctx.prec = 53
    bits = int(float(4 * arb.pi() * arb(n).sqrt() / arb(2).log()) + 40)
    k = 1
    s = 0
    ctx.prec = prec = bits
    c1 = clock()
    pi = arb.pi()
    import math
    while 1:
        term_prec = max(10, 4*math.pi*n**0.5/k/math.log(2)) + 100
        print(k, prec, term_prec)
        ctx.prec = term_prec
        term = (4*pi*arb(n).sqrt()/k).bessel_i(1)
        term_mag = float(term.log())
        term *= klooster(n,-1,k)
        term /= k
        term = term * (2*pi) / arb(n).sqrt()
        ctx.prec = prec
        s += term
        # todo: quickbound
        N = k + 1
        ctx.prec = 53
        err = 36*arb(N)**0.75 * (4*pi*arb(n).sqrt()/N).bessel_i(1)
        err *= 2*pi/arb(n).sqrt()
        ctx.prec = prec
        if k % 10 == 0:
            print(k, s.str(5), err.str(5, radius=False))
        if err < q:
            s2 = s + arb(0,err)
            v = (s2 / q).floor().unique_fmpz()
            if v is not None:
                c2 = clock()
                print(n, "IN TIME", c2 - c1)
                return v * q + fmpz(coeffn[n])
        k += 1

print(coeff(100))

print(coeff(1000))

print(coeff(10000))

print(coeff(100000))

print(coeff(1000000))

print(coeff(10**7))

print(coeff(10**8))


"""
klooster(10,-1,5)

def coeff(n,N):
    s = 0
    pi = arb.pi()
    err = 36*arb(N)**0.75 * (4*pi*arb(n).sqrt()/N).bessel_i(1)
    err *= 2*pi/arb(n).sqrt()
    for k in range(1,N+1):
        t = 4*pi*arb(n).sqrt()/k
        t = t.bessel_i(1) * klooster(n,-1,k) / k
        s += t
        print(k, 2*pi*s/arb(n).sqrt(), err)

ctx.dps = 50
coeff(2, 20)


def coeffer(n,N):
    s = 0
    pi = arb.pi()
    err = 36*arb(N)**0.75 * (4*pi*arb(n).sqrt()/N).bessel_i(1)
    err *= 2*pi/arb(n).sqrt()
    err2 = 36*arb(N)**0.75 * (4*pi*arb(n).sqrt()/N) * 0.501
    err2 *= 2*pi/arb(n).sqrt()
    err3 = 1425 / arb(N)**0.25
    print(err)
    print(err2)
    print(err3)

coeffer(1e12, 1e10)


    print("snurrigt", err / (9*arb(N)**0.25) * N**(1.067/log(log(N))))


from flint import *
from pyca import *

def powoir(a,n):
    if n == 0:
        return a**0
    if n == 1:
        return a
    if n == 2:
        return a*a
    if n % 2 == 0:
        return powoir(a,n//2)**2
    else:
        return powoir(a,n//2)**2 * a

def klooster(a,b,m):
    s = 0
    w = (2*pi*1j/m).exp()
    for x in range(m):
        if fmpz(x).gcd(m) == 1:
            if x == 0 and m == 1:
                y = 0
            else:
                y = int(1 / nmod(x, m))
            s += powoir(w, a*x + b*y)
    return s

klooster(10,-1,5)


"""

