X <- scan("parameters_all.txt", what="", sep=" ")

is.x <- TRUE
group <- -1
curve.y <- curve.x <- curve.bin <- c()
for(i in 1:length(X)){
  if(X[i] != "" & X[i] != "group"){
    if (X[i] == "]"){
      is.x <- !is.x
      cat(group, "...")
    } else if(is.x & !is.na(as.numeric(X[i]))) {
      curve.x <- c(curve.x, as.numeric(X[i]))
      curve.bin <- c(curve.bin, group)
    } else if(!is.x & !is.na(as.numeric(X[i]))) {
      curve.y <- c(curve.y, as.numeric(X[i]))
    }
  } else if (X[i] == "group"){
    group <- group + 1
  }
}
sims = data.frame(x = curve.x, y = curve.y, bin = curve.bin)
rm(group,is.x,X,curve.x,curve.y,curve.bin)

alvin_plotter = function(type=1, part = "up"){
  if(part == "up"){ range <- 8:64 } else { range <- 65:120 }
  if(type == 1){
    lines(sims$x[sims$bin==1][range], sims$y[sims$bin==1][range])
  } else{
    plot(sims$x[sims$bin==1][range], sims$y[sims$bin==1][range], xlim = c(0,45), ylim = c(0,10000), bty = "l", type = "l", xlab="i100", ylab="FUV")
  }
  for (i in 2:358){
    lines(sims$x[sims$bin==i][range], sims$y[sims$bin==i][range], col = i)
  }
}

alvin_log_plotter = function(type=1, part = "up"){
  if(part == "up"){ range <- 8:64 } else { range <- 65:120 }
  if(type == 1){
    lines(log10(sims$x[sims$bin==1][range]), log10(sims$y[sims$bin==1][range]))
  } else{
    plot(log10(sims$x[sims$bin==1][range]), log10(sims$y[sims$bin==1][range]), xlim = c(-3,3), ylim = c(1,5), bty = "l", type = "l", xlab="log10(i100)", ylab="log10(FUV)")
  }
  for (i in 2:358){
    lines(log10(sims$x[sims$bin==i][range]), log10(sims$y[sims$bin==i][range]), col = i)
  }
}

comparison.plotter = function(theta){
  par(mfrow=c(2,2))
  alvin_plotter(2)
  poly.fitter(theta, .1, P=4, zoomed = FALSE)
  alvin_plotter(1)
  poly.fitter(theta, .1, P=4)
  poly.fitter(theta, .1, P=4, zoomed = TRUE)
  alvin_plotter(1)
  par(mfrow=c(1,1))
}

overlay.plotter = function(theta, part = "down", log = T) {
  if(log){
    alvin_log_plotter(type = 2, part = part)
    points(log10(i100[i100*FUV>0 & long>theta & long <theta + .1]), log10(FUV[i100*FUV>0 & long>theta & long < theta +.1]), col = "black")
  } else {
    alvin_plotter(type = 2, part = part)
    points((i100[i100*FUV>0 & long>theta & long <theta + .1]), (FUV[i100*FUV>0 & long>theta & long < theta +.1]), col = "black")
  }
}
