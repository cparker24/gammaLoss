import matplotlib as plt
import numpy as np
import pandas as pd

def RatePlot():
    df = pd.read_csv("../build/AMY-NLO.out", header=None, delimiter=' ')
    df.columns = ["$p$ [GeV]", "C-LO", "C-NLO", "QGP: NLO", "QGP: Base", "delCcoll", "Csoftsc"]
    df["NLO/LO"] = df["QGP: NLO"]/df["QGP: Base"]
    df["x"] = df["$p$ [GeV]"] / 0.225
    df["C+delCcoll"] = df["C-LO"] + df["delCcoll"]
    df["C+Csoftsc"] = df["C-LO"] + df["Csoftsc"]

    print(df)

    plot = df.plot(title="AMY Rate", x=8, y=[7], sharex=True)
    plot.set_ylabel("$(LO + NLO)/LO$")
    plot.set_ybound(0.,2.)
    pltimage = plot.get_figure()
    pltimage.savefig("plots/others/NLOcomparison.png")

    plot = df.plot(title="C Rate", x=8, y=[1,2,9,10], sharex=True)
    plot.set_ybound(0.,5.)
    plot.set_xbound(0.5,10.)
    plot.tick_params(axis='y',which='both',direction='in')
    ax2 = plot.secondary_yaxis("right")
    ax2.tick_params(axis="y", direction="in")
    pltimage = plot.get_figure()
    pltimage.savefig("plots/others/NLO-C-comparison.png")
    

RatePlot()