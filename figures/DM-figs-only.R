load("DM-sim-data.RData")

###supp fig 1###
hist(ref_dist_test_ootX,
     xlim = c(0,30),
     xlab = "Retrogenes from X chromosome (D. melanogaster)",
     ylab = "Count",
     main = "",
     col = "green",
     cex.lab = 1.5,
     cex.axis = 1.5)
lines(c(totalX_ootX, totalX_ootX),c(0, 250), lwd = 3, col = "blue")
text(30*.8, max(hist_test_ootX$counts)*.8, paste("p = ",p_val_ootX), pos = 4, cex = 1.5)

###fig 2###
hist(sim_dists_3a,
     breaks = 15,
     main = "",
     xlab = "Distance to parental gene (D. melanogaster)",
     ylab = "Count",
     col = "green",
     cex.lab = 1.5,
     cex.axis = 1.5)
lines(c(ts_3b, ts_3b),c(-0.05, 300), lwd = 3, col = "blue")
text(((max(hist_3a$breaks)-min(hist_3a$breaks))*.8) + min(hist_3a$breaks),  
     max(hist_3a$counts)*.8, paste("p = ",pval_3c), pos = 4, cex = 1.5)

###fig 3###
hist(sim_dists_4a, 
     breaks = 15,
     main = "", 
     xlab = "Avg distance between parent and daughter(s) (D. melanogaster)", 
     ylab = "Count",
     col = "green",
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(c(ts_4b, ts_4b),c(0, 400), lwd = 3, col = "blue")
text(((max(hist_4a$breaks)-min(hist_4a$breaks))*.8) + min(hist_4a$breaks), max(hist_4a$counts)*.8, paste("p = ",pval_4c), pos = 4, cex = 1.5)


###fig 4###
hist(ref_dist_5a, 
     breaks = 15,
     main = "", 
     xlab = "Mean minimum distance (D. melanogaster)", 
     ylab = "Count",
     col = "green" ,
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(c(ts_5b, ts_5b),c(0, 400), lwd = 3, col = "blue")
text(((max(hist_5a$breaks)-min(hist_5a$breaks))*.8) + min(hist_5a$breaks), 
     max(hist_5a$counts)*.8, paste("p = ",pval_5c), pos = 4, cex = 1.5)
