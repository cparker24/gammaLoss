import matplotlib as plt
import numpy as np
import pandas as pd

def RatePlot(p):
    df = pd.read_csv("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/"+p+"dRd3p.out", header=None, delimiter=' ')
    df.columns = ["$T$ [GeV]", "normal", r"$\Psi$CE", "enhancement", "Hadronic phase", "QGP"]

    print(df)

    plot = df.plot(title="Interaction Rate for a "+p[:3]+" GeV Photon", x=0, y=[1,2,4,5], sharex=True)
    plot.set_ylabel("$dR/dp^3$ [GeV]")
    plot.set_yscale("log")
    plot.vlines(0.154,0,0.5,'black','dashed')
    pltimage = plot.get_figure()
    pltimage.savefig("plots/others/"+p[:3]+"rate.pdf")

RatePlot("0.500000")
RatePlot("1.000000")
RatePlot("2.000000")