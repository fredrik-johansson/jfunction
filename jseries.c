/*
Code for computing coefficients in the j-function.

This is ad-hoc for now; with some cleanup, could be merged
into Flint and Arb.
*/

#include "flint/profiler.h"
#include "flint/nmod_poly.h"
#include "arb.h"
#include "arb_hypgeom.h"
#include <math.h>

static void
_eta_four(mp_ptr c, slong N, nmod_t mod)
{
    slong k1, n1, k2, n2;
    mp_limb_t v, two;

    _nmod_vec_zero(c, N);
    two = 2 % mod.n;

    /* P * R */
    for (k1 = 0, n1 = 0; n1 < N; n1 += 3 * k1 + 1, k1++)
    {
        for (k2 = 0, n2 = 0, v = 1; n1 + n2 < N; )
        {
            if ((k1 + k2) % 2)
                c[n1 + n2] = nmod_sub(c[n1 + n2], v, mod);
            else
                c[n1 + n2] = nmod_add(c[n1 + n2], v, mod);

            n2 += k2 + 1;
            k2++;
            v = nmod_add(v, two, mod);
        }
    }

    /* Q * R */
    for (k1 = 1, n1 = 2; n1 < N; n1 += 3 * k1 + 2, k1++)
    {
        for (k2 = 0, n2 = 0, v = 1; n1 + n2 < N; )
        {
            if ((k1 + k2) % 2)
                c[n1 + n2] = nmod_sub(c[n1 + n2], v, mod);
            else
                c[n1 + n2] = nmod_add(c[n1 + n2], v, mod);

            n2 += k2 + 1;
            k2++;
            v = nmod_add(v, two, mod);
        }
    }
}

void
_nmod_poly_j_series(mp_ptr a, slong len, int last_only, nmod_t mod)
{
    mp_ptr b, t, u;
    slong i, j;
    mp_limb_t d, d2, d3, q;

    t = _nmod_vec_init(len);

    /* t = eta^-8 */
    _eta_four(t, len, mod);    
    _nmod_poly_mullow(a, t, len, t, len, len, mod);
    _nmod_poly_inv_series(t, a, len, len, mod);

    /* a = E4 */
    _nmod_vec_zero(a, len);
    for (i = 1; i < len; i++)
    {
        d = ((ulong) i) % mod.n;
        d2 = nmod_mul(d, d, mod);
        d3 = nmod_mul(d2, d, mod);

        for (j = i; j < len; j += i)
            a[j] = nmod_add(a[j], d3, mod);
    }
    a[0] = 1;
    _nmod_vec_scalar_mul_nmod(a + 1, a + 1, len - 1, 240 % mod.n, mod);

    u = _nmod_vec_init(len);

    /* u = E4 * eta^-8 */
    _nmod_poly_mullow(u, t, len, a, len, len, mod);

    /* t = u^2 */
    _nmod_poly_mullow(t, u, len, u, len, len, mod);

    /* a = u^3 */
    if (last_only)
    {
        int nlimbs = _nmod_vec_dot_bound_limbs(len, mod);

        NMOD_VEC_DOT(d, i, len, t[i], u[len - 1 - i], mod, nlimbs);
        a[len - 1] = d;
    }
    else
    {
        _nmod_poly_mullow(a, t, len, u, len, len, mod);
    }

    flint_free(t);
    flint_free(u);
}

/* Note: this is apparently not optimal over Z. It is better to
   construct Delta^24 and E_4^3, delaying the division until
   the end. The following code might use somewhat less
   memory though. */
void
_fmpz_poly_j_series(fmpz * a, slong len, int last_only)
{
    fmpz * b, * t, * u;
    slong i, j;
    fmpz_t d, d2, d3, q;

    t = _fmpz_vec_init(len);
    fmpz_init(d);
    fmpz_init(d2);
    fmpz_init(d3);
    fmpz_init(q);

    /* t = eta^-8 */
    _fmpz_poly_eta_qexp(t, 4, len);
    _fmpz_poly_mullow(a, t, len, t, len, len);
    _fmpz_poly_inv_series(t, a, len, len);

    /* a = E4 */
    _fmpz_vec_zero(a, len);
    for (i = 1; i < len; i++)
    {
        fmpz_set_ui(d, i);
        fmpz_mul(d2, d, d);
        fmpz_mul(d3, d2, d);

        for (j = i; j < len; j += i)
            fmpz_add(a + j, a + j, d3);
    }
    fmpz_one(a + 0);
    _fmpz_vec_scalar_mul_ui(a + 1, a + 1, len - 1, 240);

    u = _fmpz_vec_init(len);

    /* u = E4 * eta^-8 */
    _fmpz_poly_mullow(u, t, len, a, len, len);

    /* t = u^2 */
    _fmpz_poly_mullow(t, u, len, u, len, len);

    /* a = u^3 */
    if (last_only)
    {
        fmpz_zero(d);
        for (i = 0; i < len; i++)
            fmpz_addmul(d, t + i, u + len - 1 - i);
        fmpz_set(a + len - 1, d);
    }
    else
    {
        _fmpz_poly_mullow(a, t, len, u, len, len);
    }

    _fmpz_vec_clear(t, len);
    _fmpz_vec_clear(u, len);

    fmpz_clear(d);
    fmpz_clear(d2);
    fmpz_clear(d3);
    fmpz_clear(q);
}

/* Estimate size of the coefficient c_n (in fact, this is
   a rigorous upper bound). */
static double
cn_size_bits(ulong n)
{
    return 18.129440567308775239 * sqrt(n) + 1;
}

/* Compute c_n using a multimodular approach. */
void
fmpz_modular_j_coeff_multi_mod(fmpz_t res, ulong n)
{
    slong bits;
    mp_ptr t;
    mp_limb_t p;
    fmpz_t mod;
    nmod_t nmod;

    bits = cn_size_bits(n);

    t = _nmod_vec_init(n + 2);
    fmpz_init(mod);

    fmpz_zero(res);
    fmpz_one(mod);

    p = UWORD(1) << (FLINT_BITS - 1);

    while (1)
    {
        p = n_nextprime(p, 1);
        nmod_init(&nmod, p);
        TIMEIT_ONCE_START
        _nmod_poly_j_series(t, n + 2, 1, nmod);
        TIMEIT_ONCE_STOP
        fmpz_CRT_ui(res, res, mod, t[n + 1], p, 0);
        fmpz_mul_ui(mod, mod, p);
        if (fmpz_bits(mod) > bits)
            break;
    }

    fmpz_clear(mod);
    _nmod_vec_clear(t);
}

/* Stupid algorithm to compute S(a,b,m), for testing purposes. */
void
arb_kloosterman_sum_naive(arb_t res, ulong a, ulong b, ulong m, slong prec)
{
    ulong x, y;
    fmpq_t u;
    arb_t t;
    nmod_t mod;

    arb_zero(res);
    fmpq_init(u);
    arb_init(t);
    nmod_init(&mod, m);

    a %= m;
    b %= m;

    for (x = 0; x < m; x++)
    {
        if (n_gcd(m, x) == 1)
        {
            if (m == 1 && x == 0)
                y = 0;
            else
                y = n_invmod(x, m);

            fmpq_set_si(u, 2 * nmod_add(nmod_mul(a, x, mod), nmod_mul(b, y, mod), mod), m);
            arb_cos_pi_fmpq(t, u, prec);
            arb_add(res, res, t, prec);
        }
    }

    fmpq_clear(u);
    arb_clear(t);
}

/* Slightly less stupid algorithm to compute S(a,b,m). */
/* More can be done here. */
void
arb_kloosterman_sum_direct(arb_t res, ulong a, ulong b, ulong m, slong prec)
{
    ulong x, y;
    fmpq_t u;
    arb_t t;
    nmod_t mod;

    if (m == 1)
    {
        arb_one(res);
        return;
    }

    if (m == 2)
    {
        arb_set_si(res, (a+b) % 2 ? -1 : 1);
        return;
    }

    arb_zero(res);
    fmpq_init(u);
    arb_init(t);
    nmod_init(&mod, m);

    a %= m;
    b %= m;

    for (x = 0; x <= m / 2; x++)
    {
        if (n_gcd(m, x) == 1)
        {
            if (m == 1 && x == 0)
                y = 0;
            else
                y = n_invmod(x, m);

            fmpq_set_si(u, 2 * nmod_add(nmod_mul(a, x, mod), nmod_mul(b, y, mod), mod), m);
            arb_cos_pi_fmpq(t, u, prec);
            arb_add(res, res, t, prec);
        }
    }

    arb_mul_2exp_si(res, res, 1);
    fmpq_clear(u);
    arb_clear(t);
}

/* Kloosterman sum using factorization. Todo: use exact formula
   for prime powers; detect when the sum is zero. */
void
arb_kloosterman_sum(arb_t res, ulong a, ulong b, ulong m, slong prec)
{
    n_factor_t fac;
    ulong m1, m2, n1, n2, p, e;
    slong i;
    arb_t t;

    if (m <= 2)
    {
        arb_kloosterman_sum_direct(res, a, b, m, prec);
        return;
    }

    arb_init(t);
    n_factor_init(&fac);
    n_factor(&fac, m, 1);

    arb_one(res);

    for (i = 0; i < fac.num; i++)
    {
        p = fac.p[i];
        e = fac.exp[i];
        m1 = n_pow(p, e);
        m2 = m / m1;

        n1 = (m2 == 1) ? 0 : n_invmod(m1 % m2, m2);
        n2 = n_invmod(m2 % m1, m1);

        arb_kloosterman_sum_direct(t, n_mulmod2(a, n2, m1), n_mulmod2(n2, b, m1), m1, prec);
        arb_mul(res, res, t, prec);

        a = n_mulmod2(n1, a, m2);
        b = n_mulmod2(n1, b, m2);
        m = m2;
    }

    arb_clear(t);
}

/* Compute a congruence for c_n with modulus of size at least bits.
   Uses congruences for the primes 2, 3, 5 and 7, and
   a power series as a fallback. */
void
fmpz_modular_j_coeff_mod(fmpz_t res, fmpz_t mod, ulong n, slong bits)
{
    mp_ptr t;
    mp_limb_t p;
    nmod_t nmod;
    fmpz_t c, d, e, f;
    ulong a, m;

    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(f);

    fmpz_zero(res);
    fmpz_one(mod);

    if (n % 2 == 0)
    {
        a = 0;
        m = n;
        while (m % 2 == 0)
        {
            m /= 2;
            a++;
        }

        fmpz_one(mod);
        fmpz_mul_2exp(mod, mod, 3 * a + 13);

        fmpz_set_ui(res, m);
        fmpz_divisor_sigma(res, res, 7);
        fmpz_set_ui(c, 3);
        fmpz_pow_ui(c, c, a - 1);
        fmpz_mul(res, res, c);
        fmpz_mul_2exp(res, res, 3 * a + 8);
        fmpz_neg(res, res);
        fmpz_mod(res, res, mod);
    }

    if (n % 3 == 0)
    {
        a = 0;
        m = n;
        while (m % 3 == 0)
        {
            m /= 3;
            a++;
        }

        fmpz_set_ui(d, m);
        fmpz_divisor_sigma(d, d, 1);

        fmpz_set_ui(e, 10);
        fmpz_pow_ui(e, e, a - 1);
        fmpz_mul(d, d, e);

        fmpz_set_ui(e, 3);
        fmpz_pow_ui(e, e, 2*a+3);
        fmpz_mul(d, d, e);

        if (m % 3 == 1)
            fmpz_neg(d, d);

        fmpz_set_ui(e, 3);
        fmpz_pow_ui(e, e, 2*a + 6);

        fmpz_set_ui(f, m);
        fmpz_invmod(f, f, e);
        fmpz_mul(d, d, f);

        fmpz_CRT(res, res, mod, d, e, 0);
        fmpz_mul(mod, mod, e);
    }

    if (n % 5 == 0)
    {
        a = 0;
        m = n;
        while (m % 5 == 0)
        {
            m /= 5;
            a++;
        }

        fmpz_set_ui(d, m);
        fmpz_divisor_sigma(d, d, 1);
        fmpz_mul_ui(d, d, m);
        fmpz_set_ui(e, 3);
        fmpz_pow_ui(e, e, a - 1);
        fmpz_mul(d, d, e);
        fmpz_set_ui(e, 5);
        fmpz_pow_ui(e, e, a + 1);
        fmpz_mul(d, d, e);
        fmpz_neg(d, d);

        fmpz_set_ui(e, 5);
        fmpz_pow_ui(e, e, a + 2);

        fmpz_CRT(res, res, mod, d, e, 0);
        fmpz_mul(mod, mod, e);
    }

    if (n % 7 == 0)
    {
        a = 0;
        m = n;
        while (m % 7 == 0)
        {
            m /= 7;
            a++;
        }

        fmpz_set_ui(d, m);
        fmpz_divisor_sigma(d, d, 3);
        fmpz_mul_ui(d, d, m);

        fmpz_set_ui(e, 5);
        fmpz_pow_ui(e, e, a - 1);
        fmpz_mul(d, d, e);
        fmpz_set_ui(e, 7);
        fmpz_pow_ui(e, e, a);
        fmpz_mul(d, d, e);
        fmpz_neg(d, d);

        fmpz_set_ui(e, 7);
        fmpz_pow_ui(e, e, a + 1);

        fmpz_CRT(res, res, mod, d, e, 0);
        fmpz_mul(mod, mod, e);
    }

    p = UWORD(1) << FLINT_MIN(FLINT_BITS - 1, bits);
    p = FLINT_MAX(p, 12);

    while (fmpz_bits(mod) < bits)
    {
        p = n_nextprime(p, 1);
        nmod_init(&nmod, p);

        t = _nmod_vec_init(n + 2);
        _nmod_poly_j_series(t, n + 2, 1, nmod);
        fmpz_CRT_ui(res, res, mod, t[n + 1], p, 0);
        _nmod_vec_clear(t);
        fmpz_mul_ui(mod, mod, p);
    }

    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(f);
}

/* Bound I_1(x) */
void
mag_bessel_i1(mag_t res, mag_t x)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);
    if (mag_cmp_2exp_si(x, 2) < 0)
    {
        /* 0.5*x*(1+x^2/4) */
        mag_mul(t, x, x);
        mag_mul_2exp_si(t, t, -2);
        mag_add_ui(t, t, 1);
        mag_mul(t, t, x);
        mag_mul_2exp_si(res, t, -1);
    }
    else
    {
        /* e^x / sqrt(2pix) */
        mag_rsqrt(t, x);
        mag_exp(u, x);
        mag_mul(res, t, u);
        mag_mul_ui(res, res, 26146);
        mag_mul_2exp_si(res, res, -16);
    }
}

/* Error bound for the Petersson-Rademacher series. */
static void
error_bound(mag_t bound, ulong n, ulong N)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);

    /* 72*pi/sqrt(n) * N^(3/4) * I_1(4*pi*sqrt(n)/N) */

    mag_const_pi(t);
    mag_mul_2exp_si(t, t, 2);
    mag_set_ui(u, n);
    mag_sqrt(u, u);
    mag_mul(t, t, u);
    mag_div_ui(t, t, N);
    mag_bessel_i1(t, t);

    mag_set_ui(u, N);
    mag_sqrt(u, u);
    mag_sqrt(u, u);
    mag_pow_ui(u, u, 3);
    mag_mul(t, t, u);

    mag_mul_ui(t, t, 72);
    mag_const_pi(u);
    mag_mul(t, t, u);
    mag_set_ui_lower(u, n);
    mag_rsqrt(u, u);
    mag_mul(t, t, u);

    mag_set(bound, t);

    mag_clear(t);
    mag_clear(u);
}

/* Hybrid algorithm: compute c_n given r = c_n mod M, 0 <= r < M. */
void
fmpz_modular_j_coeff_given_residue(fmpz_t res, ulong n, const fmpz_t r, const fmpz_t M)
{
    arb_t s, s2, t, u, v, c;
    fmpz_t d;
    mag_t err, mod_mag;
    double term_mag;
    ulong k;
    slong prec, prec2, prec3, actual_initial_bits;

    arb_init(c);
    arb_init(s);
    arb_init(s2);
    arb_init(t);
    arb_init(u);
    arb_init(v);
    mag_init(err);
    mag_init(mod_mag);

    fmpz_init(d);

    actual_initial_bits = fmpz_bits(M);

    prec = 18.129440567309 * sqrt(n) + 10;
    prec3 = prec;

    for (k = 1; ; k++)
    {
        term_mag = 18.129440567309 * sqrt(n) / k;
        prec2 = term_mag + 10 + FLINT_BIT_COUNT(n);

        arb_kloosterman_sum(u, n % k, (k-1) % k, k, prec2);
        arb_div_ui(u, u, k, prec2);

        /* Todo: when u contains 0, skip the Bessel function
           (only need upper bound). */

        arb_sqrt_ui(t, n, prec2);
        arb_const_pi(v, prec2);
        arb_mul(t, t, v, prec2);
        arb_mul_2exp_si(t, t, 2);
        arb_div_ui(t, t, k, prec2);
        arb_one(c);

        arb_hypgeom_bessel_i(t, c, t, prec2);
        arb_mul(u, u, t, prec2);
        arb_add(s2, s2, u, prec3);

        /* Only update main sum and check for convergence
           for a sparse subset of k (not tuned). */
        if ((k < 300 && (k % 4 == 0)) || (k % 64 == 0))
        {
            arb_add(s, s, s2, prec);
            arb_zero(s2);
            prec3 = prec2 + 100;

            error_bound(err, n, k);
            mag_set_fmpz_lower(mod_mag, M);

            if (mag_cmp(err, mod_mag) < 0)
            {
                arb_const_pi(t, prec);
                arb_mul(v, s, t, prec);
                arb_sqrt_ui(t, n, prec);
                arb_div(v, v, t, prec);
                arb_mul_2exp_si(v, v, 1);
                arb_add_error_mag(v, err);

                arb_sub_fmpz(v, v, r, prec);
                arb_div_fmpz(v, v, M, prec);

                if (arb_get_unique_fmpz(d, v))
                {
                    if (n >= 1e5 && 0)
                    {
                        printf("%lu  bound, mag: ", k);
                        mag_printd(err, 5); printf("  "); mag_printd(mod_mag, 5); printf("\n");
                    }

                    fmpz_mul(d, d, M);
                    fmpz_add(res, d, r);
                    break;
                }
            }
        }
    }

    fmpz_clear(d);

    arb_clear(c);
    arb_clear(s);
    arb_clear(s2);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    mag_clear(err);
    mag_clear(mod_mag);
}

void
fmpz_modular_j_coeff_hybrid(fmpz_t res, ulong n, slong initial_bits)
{
    fmpz_t r, M;

    fmpz_init(r);
    fmpz_init(M);

    fmpz_modular_j_coeff_mod(r, M, n, initial_bits);
    fmpz_modular_j_coeff_given_residue(res, n, r, M);

    fmpz_clear(r);
    fmpz_clear(M);
}

void
test_kloosterman_sum()
{
    flint_rand_t state;
    flint_randinit(state);
    arb_t t, u;
    ulong a, b, m;
    slong iter;

    arb_init(t);
    arb_init(u);

    for (iter = 0; iter < 10000; iter++)
    {
        a = n_randint(state, 1000);
        b = n_randint(state, 1000);
        m = 1 + n_randint(state, 300);

        arb_kloosterman_sum_naive(t, a, b, m, 100);
        arb_kloosterman_sum(u, a, b, m, 100);

        if (!arb_overlaps(t, u))
            flint_abort();
    }

    arb_clear(t);
    arb_clear(u);
}

void
search_primes(ulong nmin, ulong nmax)
{
    fmpz_t c, r, M;
    slong len;
    mp_ptr t;
    ulong n, P, cm;
    nmod_t mod;
    int prime;
    fmpz_init(c);
    fmpz_init(r);
    fmpz_init(M);

    len = nmax + 2;

    P = UWORD(614889782588491410);
    nmod_init(&mod, P);

    t = _nmod_vec_init(len);
    _nmod_poly_j_series(t, len, 0, mod);

    for (n = nmin; n <= nmax; n++)
    {
        cm = t[n + 1];

        if (n_gcd(cm, P) == 1)
        {
            printf("Candidate: %lu\n", n);

            fmpz_set_ui(M, P);
            fmpz_set_ui(r, cm);
            fmpz_modular_j_coeff_given_residue(c, n, r, M);
            prime = fmpz_is_probabprime(c);
            printf("is prime: %d\n", prime);
            if (prime)
            {
                printf("\n");
                fmpz_print(c);
                printf("\n\n");
            }
        }
    }

    _nmod_vec_clear(t);

    fmpz_clear(c);
    fmpz_clear(r);
    fmpz_clear(M);
}

int main()
{
    /* test_kloosterman_sum(); */
    /* search_primes(457871, 457871); */

    /* search_primes(1, 1.5e7); */
}
