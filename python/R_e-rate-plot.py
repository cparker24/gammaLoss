import matplotlib as plt
import numpy as np
import pandas as pd


df = pd.read_csv("../build/R_e.out", header=None, delimiter=' ')
df.columns = ["$T$ [GeV]", "normal", "$\Psi$CE", "enhancement", "Hadronic phase", "QGP"]

print(df)

plot = df.plot(title="Emission Rate", x=0, y=[1,2,4,5], sharex=True)
plot.set_ylabel("$R_e$ [GeV$^4$]")
plot.set_yscale("log")
pltimage = plot.get_figure()
pltimage.savefig("plots/others/Emission-rate.pdf")