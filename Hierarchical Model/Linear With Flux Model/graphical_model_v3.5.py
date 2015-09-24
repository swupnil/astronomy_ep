from matplotlib import rc
rc("font", family="serif", size=12)
rc("text", usetex=True)
rc("./weaklensing.tex")

import daft

pgm = daft.PGM([4.7, 4.75], origin=[-1.35, 2.15])
#intercept
pgm.add_node(daft.Node("A1", r"$A_1$", -1, 6))
pgm.add_node(daft.Node("tauA", r"$\tau_\alpha$", -1, 5))
pgm.add_node(daft.Node("A2", r"$A_2$", -1, 4))
pgm.add_node(daft.Node("alpha", r"$\alpha_{i}$", 0, 5))

#slope
pgm.add_node(daft.Node("B1", r"$B_1$", 3, 6))
pgm.add_node(daft.Node("tauB", r"$\tau_\beta$", 3, 5))
pgm.add_node(daft.Node("B2", r"$B_2$", 3, 4))
pgm.add_node(daft.Node("beta", r"$\beta_{i}$", 2, 5))

#variance
pgm.add_node(daft.Node("sigma", r"$\sigma_i^2$", 2, 3.5))
pgm.add_node(daft.Node("muS", r"$\mu_{\sigma}$", 1.5, 2.5))
pgm.add_node(daft.Node("tauS", r"$\tau_{\sigma}$", 2.5, 2.5))

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

pgm.add_edge("tauA", "alpha")
pgm.add_edge("A1", "alpha")
pgm.add_edge("A2", "alpha")
pgm.add_edge("alpha", "obs")

pgm.add_edge("tauB", "beta")
pgm.add_edge("B1", "beta")
pgm.add_edge("B2", "beta")
pgm.add_edge("beta", "obs")

pgm.add_edge("muS", "sigma")
pgm.add_edge("tauS", "sigma")
pgm.add_edge("sigma", "obs")

pgm.add_edge("x", "obs")

pgm.render()
pgm.figure.savefig("LightIntensity_Model.png",dpi=1000)