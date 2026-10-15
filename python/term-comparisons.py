import math
import matplotlib.pyplot as plt
import numpy as np

def alphaS(temp):
  return 6.0 * math.pi / (27.0 * math.log(temp / 0.022))

def gS(temp):
  return math.sqrt(4.0 * math.pi * alphaS(temp))

def f(x):
   return x/(1+math.exp(x))

def g(x):
   C22 = 0.041 / x - 0.3615 + 1.01 * math.exp(-1.35 * x)
   Cab = math.sqrt(1.5) * (0.548 / pow(x, 1.5) * math.log(12.28 + 1.0 / x) + 0.133 * x / math.sqrt(1.0 + x / 16.27))
  
   return (0.5 * math.log(2.0 * x) + C22 + Cab) * f(x)

def fterm(x,T):
   return math.log(math.sqrt(3)/gS(T))*f(x)

xs = np.linspace(0.01,10,100)

# plot the data
plt.figure(figsize=(5,5))
plt.plot(xs, [fterm(x,0.16) for x in xs], label="f for 0.16")
plt.plot(xs, [fterm(x,0.20) for x in xs], label="f for 0.20")
plt.plot(xs, [fterm(x,0.25) for x in xs], label="f for 0.25")
plt.plot(xs, [0.09*f(x) for x in xs], label="AMY f")
plt.plot(xs, [g(x) for x in xs], label="g")
plt.plot(xs, [fterm(x,0.20) + g(x) for x in xs], label="sum for 0.20")
plt.plot(xs, [0.09*f(x) + g(x) for x in xs], label="sum for AMY f")
plt.yscale("log")
plt.legend()
plt.xlabel(r"$x$")
plt.ylabel(r"term strength")
plt.savefig('plots/term-comparison.png', dpi=1000)