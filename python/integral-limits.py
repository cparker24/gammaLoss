import matplotlib as plt
import numpy as np
import pandas as pd


df = pd.read_csv("../build/170pdRd3p.out", header=None, delimiter=' ')
df.columns = ["$p$ (GeV)", "HRG", r"$ \pi + \rho + a_1$", "$\omega$", "strangeness", "meson bremsstrahlung", "QGP"]
df["Hadron Gas Total"] = df["HRG"] + df[r"$ \pi + \rho + a_1$"] + df["$\omega$"] + df["strangeness"] + df["meson bremsstrahlung"]
df["meson tree channels"] = df[r"$ \pi + \rho + a_1$"] + df["$\omega$"] + df["strangeness"]*1.0
df["QGP"] = df["QGP"]/df["$p$ (GeV)"]

print(df)
x = df["$p$ (GeV)"].to_numpy()

plot = df.plot(title="Hadron Gas Rates 170 MeV", x=0, y=[6], sharex=True)
plot.fill_between(x,df["QGP"].to_list(),where=(0.5 <= x)&(x <= 2.5), alpha=0.2)
plot.set_ylabel("$dR/d^3p$ (fm$^{-4}$GeV$^{-3}$)")
plot.set_yscale("log")
plot.set_xbound(0.2,3)
#plot.set_ybound(1e-9,1e-3)
pltimage = plot.get_figure()
pltimage.savefig("plots/QGP 170 MeV.png", dpi=300)