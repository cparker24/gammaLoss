import matplotlib as plt
import numpy as np
import pandas as pd

def MakePlot(p):
    df = pd.read_csv("../build/"+p+"pdRd3p.out", header=None, delimiter=' ')
    df.columns = ["$p$ [GeV]", r"in medium $\rho$", r"$\omega$ t channel", r"meson bremsstrahlung", r"$\pi \omega$", r"$\rho \omega$", r"$\pi \rho$", r"$\pi + \rho + a_1$", "strangeness", "Hadron Gas Total", "QGP: NLO", "QGP: LO"]
    #df["Hadron Gas Total"] = df["HRG"] + df[r"$ \pi + \rho + a_1$"] + df["$\omega$"] + df["strangeness"] + df["meson bremsstrahlung"]
    #df["meson tree channels"] = df[r"$ \pi + \rho + a_1$"] + df["$\omega$"] + df["strangeness"]*1.0

    print(df)

    plot = df.plot(title="Photon Interaction Rates "+p+" MeV", x=0, y=[1,2,3,4,5,6,7,8,9,10], sharex=True)
    plot.set_ylabel("$p*d^3R/dp^3$ [fm$^{-4}$GeV$^{-2}$]")
    plot.set_yscale("log")
    plot.set_xbound(0.2,3)
    plot.set_ybound(1e-15,1e-3)
    plot.tick_params(axis='y',which='both',direction='in')
    ax2 = plot.secondary_yaxis("right")
    ax2.tick_params(axis="y", direction="in")
    pltimage = plot.get_figure()
    pltimage.savefig("plots/others/Rates "+p+" MeV.png", dpi=3000)

MakePlot("154")
MakePlot("170")
MakePlot("200")