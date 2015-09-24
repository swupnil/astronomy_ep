# plot 1-6 degree polynomials for a particular bin
poly.plot(90, 12)

# plot 4th degree fit in a particular bin and alvin's plots over
comparison.plotter(150)
comparison.plotter(270)
comparison.plotter(90)
comparison.plotter(30)

# plot all simulation curves with particular bin's data overlayed
overlay.plotter(30, "down")
overlay.plotter(150, "up")
overlay.plotter(30, "down", log = F)
overlay.plotter(150, "up", log = F)
overlay.plotter(300, "up", log = F)

# fit third degree polynomial to each simulation bin
par(mfrow=c(3,3))
for(i in 350:359){
  sim.fitter(degrees = i, inc = 1, P = 7)
}

# fifth degree works well except for 15-84
# dip around 50 in simulations..

# try to discover some good initial values
poly.fitter(150, .1, P=4)
a = c(180, 350, 5, .5, 30, 7, 0, 40, .5)
lines((0:500)/10,exp.curve((0:500)/10,a), col ="red")
a = c(320, 2, 2.5, 350, 5, .5, 30, 7)
lines((0:500)/10,exp.curve2((0:500)/10,a), col ="red")

poly.fitter(270, .1, P=4)
a = c(180, 750, 5, .5, 50, 20, 1, 55, .5)
lines((0:500)/10,exp.curve((0:500)/10,a), col ="red")
a = c(480, 4, 1, 750, 5, .5, 50, 20)
lines((0:500)/10,exp.curve2((0:500)/10,a), col ="red")
