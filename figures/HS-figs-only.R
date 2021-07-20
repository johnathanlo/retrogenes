load("HS-sim-data.RData")

###supp fig 1###
hist(ref_dist_2a,
     xlim = c(0,30),
     main = "",
     xlab = "Retrogenes from X chromosome (human)",
     ylab = "Count",
     col = "green",
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(density_2a,
        col = rgb(0,0,1,.2))
lines(c(totalX_2b, totalX_2b),c(0, 1000), lwd = 3, col = "blue")
text(30*.8, max(hist_ootX$counts)*.8, paste("p = ",p_val_2c), pos = 4, cex = 1.5)

###fig 2###
hist(ref_dist_3a,
     breaks = 15,
     main = "",
     xlab = "Distance to parental gene (human)",
     ylab = "Count",
     col = "green",
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(c(ts_3b, ts_3b),c(-0.05, 420), lwd = 3, col = "blue")
text(((max(hist_3a$breaks)-min(hist_3a$breaks))*.8) + min(hist_3a$breaks),
     max(hist_3a$counts)*.8, paste("p = ",pval_3c), pos = 4, cex = 1.5)

###fig 3###
hist(sim_dists_4a,
     breaks = 15,
     main = "",
     xlab = "Avg distance between parent and daughter(s) (human)",
     ylab = "Count",
     col = "green",
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(c(ts_4b, ts_4b),c(-0.05, 400), lwd = 3, col = "blue")
text(((max(hist_4a$breaks)-min(hist_4a$breaks))*.8) + min(hist_4a$breaks), 
     max(hist_4a$counts)*.8, paste("p = ",pval_4c), pos = 4, cex = 1.5)


###fig 4###
hist(ref_dist_5a, 
     breaks = 15,
     main = "", 
     xlab = "Mean minimum distance (human)", 
     ylab = "Count",
     col = "green",
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(c(ts_5b, ts_5b),c(-0.05, 400), lwd = 3, col = "blue")
text(((max(hist_5a$breaks)-min(hist_5a$breaks))*.8) + min(hist_5a$breaks), 
     max(hist_5a$counts)*.8, paste("p = ",pval_5c), pos = 4, cex = 1.5)
