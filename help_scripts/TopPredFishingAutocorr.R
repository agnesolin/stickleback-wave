
#### how similar are values across space ####

spat.cor1 = spline.correlog(df_sub$X, df_sub$Y, df_sub$totalTopPred,  resamp = 10, na.rm = T)
spat.cor2 = spline.correlog(df_sub$X, df_sub$Y, df_sub$fishing,  resamp = 10, na.rm = T)


jpeg("suppFigures/PredFishAutocorr.jpg", width = 18.5, height = 24, units = "cm", res = 300)

par(mfrow = c(2,1))

plot(spat.cor1,  xlab = "Distance (km)", ylab = "Correlation predation pressure", xaxt="n")
axis(1, at = seq(0, 500000, 100000), labels = seq(0, 500, 100), las=1)

interc = spat.cor1$real$x.intercept
abline(v = interc, lty = 3)
text(interc, 0.8, pos = 4,
     paste(round(interc/1000), "km"))
text(0,0.95, "a", cex = 1.3)

plot(spat.cor2,  xlab = "Distance (km)", ylab = "Correlation fishing", xaxt="n")
axis(1, at = seq(0, 500000, 100000), labels = seq(0, 500, 100), las=1)

interc = spat.cor2$real$x.intercept
abline(v = interc, lty = 3)
text(interc, 0.8, pos = 4,
     paste(round(interc/1000), "km"))
text(0,0.95, "b", cex = 1.3)

dev.off()
