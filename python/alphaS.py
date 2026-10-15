import matplotlib.pyplot as plt
import math
import numpy as np

def alphaSOld(temp):
  return 6.0 * math.pi / (27.0 * math.log(temp / 0.022))

def alphaS(temp):
  return 6.0 * math.pi / (27.0 * math.log(3.2* temp / 0.2))

def alphaSfixed(temp):
  return 0.3

Ts = np.linspace(100,900,100)

plt.figure(figsize=(5,5))
plt.plot(Ts, [alphaS(T/1000) for T in Ts], label=r"$\alpha_S$")
#plt.plot(Ts, [alphaSOld(T/1000) for T in Ts], label=r"old $\alpha_S$")
plt.plot(Ts, [alphaSfixed(T/1000) for T in Ts], label=r"fixed $\alpha_S$")
plt.ylim(0,0.5)
plt.legend()
plt.vlines(154,0,0.5,'black','dashed')
plt.title(r"Running $\alpha_S$ Over Relevant Temperatures")
plt.xlabel(r"$T$ [MeV]")
plt.ylabel(r"$\alpha_S$")
plt.savefig('plots/others/alphaS.pdf', dpi=1000)