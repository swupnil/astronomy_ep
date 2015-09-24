from matplotlib import rc
rc("font", family="serif", size=12)
rc("text", usetex=True)
rc("./weaklensing.tex")

import daft

pgm = daft.PGM([4.7, 5.25], origin=[-1.35, 2.15])
#intercept
pgm.add_node(daft.Node("A", r"$A_i$", 0, 6))
pgm.add_node(daft.Node("muA", r"$\mu_A$", 0, 7))
pgm.add_node(daft.Node("tauA", r"$\tau_A$", -1, 6))
pgm.add_node(daft.Node("B", r"$B_i$", 0, 4.25))
pgm.add_node(daft.Node("muB", r"$\mu_B$", -1, 5))
pgm.add_node(daft.Node("tauB", r"$\tau_D$", -1, 4.25))
pgm.add_node(daft.Node("alpha", r"$\alpha_{i}$", 0, 5))

#slope
pgm.add_node(daft.Node("C", r"$C_i$", 2, 6))
pgm.add_node(daft.Node("muC", r"$\mu_C$", 2, 7))
pgm.add_node(daft.Node("tauC", r"$\tau_C$", 3, 6))
pgm.add_node(daft.Node("D", r"$D_i$", 2, 4.25))
pgm.add_node(daft.Node("muD", r"$\mu_D$", 3, 5))
pgm.add_node(daft.Node("tauD", r"$\tau_D$", 3, 4.25))
pgm.add_node(daft.Node("beta", r"$\beta_{i}$", 2, 5))

#variance
pgm.add_node(daft.Node("sigma", r"$\sigma_i^2$", 2, 3.5))
pgm.add_node(daft.Node("muS", r"$\mu_{\sigma}$", 3, 3.5))
pgm.add_node(daft.Node("tauS", r"$\tau_{\sigma}$", 2, 2.5))

#light intensity
pgm.add_node(daft.Node("fuv",r"$f_{i}$", 1, 6, observed=True))

pgm.add_node(daft.Node("obs", r"$y_{ij}$", 1, 5, observed=True))
pgm.add_node(daft.Node("x", r"$x_{ij}$", 1, 4, observed=True))
pgm.add_plate(daft.Plate([0.5, 3.65, 1, 1.75],
                         label=r"$J$"))
pgm.add_plate(daft.Plate([-0.5, 3, 3, 3.5],
                         label=r"Longitude $I$"))

#arrows
pgm.add_edge("fuv", "alpha")
pgm.add_edge("fuv", "beta")

pgm.add_edge("muA", "A")
pgm.add_edge("tauA", "A")
pgm.add_edge("A", "alpha")
pgm.add_edge("muB", "B")
pgm.add_edge("tauB", "B")
pgm.add_edge("B", "alpha")
pgm.add_edge("alpha", "obs")

pgm.add_edge("muC", "C")
pgm.add_edge("tauC", "C")
pgm.add_edge("C", "beta")
pgm.add_edge("muD", "D")
pgm.add_edge("tauD", "D")
pgm.add_edge("D", "beta")
pgm.add_edge("beta", "obs")

pgm.add_edge("muS", "sigma")
pgm.add_edge("tauS", "sigma")
pgm.add_edge("sigma", "obs")

pgm.add_edge("x", "obs")

pgm.render()
pgm.figure.savefig("LightIntensity_Model.png",dpi=1000)