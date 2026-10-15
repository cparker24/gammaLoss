import matplotlib as plt
import numpy as np
import pandas as pd


df = pd.read_csv("../build/200strangepdRd3p.out", header=None, delimiter=' ')
df.columns = ["$p$ (GeV)", r"$\pi+K*$", r"$\pi+K$", r"$\rho+K$", r"$K+K*$"]
df["Total"] = df[ r"$\pi+K*$"] + df[r"$\pi+K$"] + df[r"$\rho+K$"] + df[r"$K+K*$"]

print(df)

plot = df.plot(title="Hadron Gas Strange Rates 200 MeV", x=0, y=[1,2,3,4,5], sharex=True)
plot.set_ylabel("$p*dR/d^3p$ (fm$^{-4}$GeV$^{-2}$)")
plot.set_yscale("log")
plot.set_xbound(0,3)
plot.set_ybound(1e-11,3e-3)
pltimage = plot.get_figure()
pltimage.savefig("Hadron Gas Strange Rates 200 MeV.png")