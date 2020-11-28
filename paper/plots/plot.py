from pylab import *

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True

ns = []
err1000 = []
bound1000 = []
err1000000 = []
bound1000000 = []

for line in open("exp1000.txt").readlines():
    n, err, bound = line.split()
    n = int(n)
    if n < 20 or (n > 1000 and n % 3 != 0):
        continue
    ns.append(float(n))
    from mpmath import pi, besseli, sqrt
    n = float(n)
    bound = 72*pi/sqrt(1000)*n**0.75*besseli(1, 4*pi*sqrt(1000)/n)
    err1000.append(float(err))
    bound1000.append(float(bound))

for line in open("exp1000000.txt").readlines():
    n, err, bound = line.split()
    n = int(n)
    if n < 20 or (n > 1000 and n % 3 != 0):
        continue
    n = float(n)
    bound = 72*pi/sqrt(1000000)*n**0.75*besseli(1, 4*pi*sqrt(1000000)/n)
    err1000000.append(float(err))
    bound1000000.append(float(bound))

subplot(121)

title("$n = 10^3$")
plot(ns, bound1000, label="Bound", linewidth=3)
loglog(ns, err1000, label="Actual error")
ylim([1e-3,1e4])
xlim([10,1e4])
axhline(0.5, color="black")
xlabel("$N$")

subplot(122)

title("$n = 10^6$")
plot(ns, bound1000000, linewidth=3, label="Bound")
loglog(ns, err1000000, label="Actual error")
legend()
ylim([1e-6,1e6])
xlim([1e2,1e4])
axhline(0.5, color="black")
xlabel("$N$")

tight_layout()

fig = plt.gcf()
fig.set_size_inches(6, 3.5)

savefig("error.pdf", dpi=200)
savefig("error.eps", dpi=200)

