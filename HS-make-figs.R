
load("HS-final-data.RData") #final dataset from 7/14/2021
runs=1000 #set # runs for simulations and bootstraps
set.seed(6)
library(dplyr)
library(rentrez)
library(chromPlot)
library(GenomicFeatures)
library(ggplot2)
library(reshape2)


##distance function
find_dist <- function(chr, loc1, loc2){
  if(chr %in% c("X", "Y")){
    chr = "23"
  }
  rec_chrom <- human_rec_map[human_rec_map$chr == chr, ]
  loc1_match <- human_rec_map[which.min(abs(loc1 - rec_chrom$phys.loc)),]
  loc2_match <- human_rec_map[which.min(abs(loc2 - rec_chrom$phys.loc)),]
  rec_dist <- as.numeric(abs(loc1_match$gen.loc - loc2_match$gen.loc))
  return(rec_dist)
}

#total number of partner genes
all_partners <- c()
for(i in 1:length(retrogenes_final)){
  partners_i <- retrogenes_final[[i]][[1]]$partnergenes
  all_partners <- c(all_partners, partners_i)
}
length(unique(all_partners))


###HS ootX


#build reference distribution for set of just protein coding retrogenes
ref_dist_2a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  sim_parents <- sample(c(1:22, "X", "Y"), 
                        size = length(retrogenes_final_pc), 
                        prob = human_genes_per_chrom$chrom.size,
                        replace = TRUE)
  ref_dist_2a[i] <- sum(sim_parents == "X")
}
density_2a <- density(ref_dist_2a, adjust = 1.4)
hist_ootX <- hist(ref_dist_2a)
plot(density_2a)

#test statistic
totalX_2b <- 0
for(i in 1:length(retrogenes_final_pc)){
  if(retrogenes_final_pc[[i]][[1]]$parentalchrom == "X"){
    totalX_2b <- totalX_2b +1
  }
}
#calculate p-values and plot on histogram 
p_val_2c <- sum(ref_dist_2a >= totalX_2b)/runs * 2 #test statistic less than median; multiply by 2 b/c test should be two-tailed


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


###########################fig 2#####################
######################################################

#make reference distribution 
ref_dist_3a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  retro_dists <- vector(mode = "logical", length = length(retrogenes_final_pc_noX))
  for(j in 1:length(retrogenes_final_pc_noX)){
    parent_chrom <- retrogenes_final_pc_noX[[j]][[1]]$parentalchrom
    parent_loc <- retrogenes_final_pc_noX[[j]][[1]]$parentalloc
    retro_chrom <- sample(row.names(human_genome), size = 1, prob = human_genome$chrom.size)
    retro_loc <- sample(1:human_genome$chrom.size[row.names(human_genome) == retro_chrom], size = 1)
    if(parent_chrom == retro_chrom){
      retro_dists[j] <- ifelse(find_dist(parent_chrom, parent_loc, retro_loc)/100 >.5,
                               .5, 
                               find_dist(parent_chrom, parent_loc, retro_loc)/100)
    }else{
      retro_dists[j] <- .5
    }
    
  }
  ref_dist_3a[i] <- sum(retro_dists)/length(retrogenes_final_pc_noX)
}
density_3a <- density(ref_dist_3a, adjust = 1.5)
plot(density_3a)
hist_3a <- hist(ref_dist_3a, breaks = 15)

#test stat

retro_dists_test <- vector(mode = "logical", length = length(retrogenes_final_pc_noX))
for(i in 1:length(retrogenes_final_pc_noX)){
  parent_chrom <- retrogenes_final_pc_noX[[i]][[1]]$parentalchrom
  parent_loc <- retrogenes_final_pc_noX[[i]][[1]]$parentalloc
  retro_chrom <- retrogenes_final_pc_noX[[i]][[1]]$retrochrom
  retro_loc <- retrogenes_final_pc_noX[[i]][[1]]$retroloc
  if(parent_chrom == retro_chrom){
    retro_dists_test[i] <- ifelse(find_dist(parent_chrom, parent_loc, retro_loc)/100 >.5,
                                  .5, 
                                  find_dist(parent_chrom, parent_loc, retro_loc)/100)
  }else{
    retro_dists_test[i] <- .5
  }
}
ts_3b <- sum(retro_dists_test)/length(retrogenes_final_pc_noX)

#calc p-value and plot
pval_3c <- sum(ref_dist_3a<=ts_3b)/runs*2

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




###########################fig 3######################
######################################################
sim_dists_4a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  avgdist <- 0
  n <- 0
  for(j in 1:length(retrogenes_by_parentals_pc)){
    tot_retro <- sum(retrogenes_by_parentals_pc[[j]]$retroid %in% retrogenes_pc$retrodb_id)
    if(tot_retro && retrogenes_by_parentals_pc[[j]]$parentalchrom != "X"){
      sim_rets <- sample(row.names(human_genome), 
                         prob = human_genome$chrom.size, 
                         size = tot_retro, 
                         replace = TRUE)
      sim_ret_locs <- vector(mode = "numeric", length=tot_retro)
      for(k in 1:tot_retro){
        sim_ret_locs[k] <- sample(1:human_genome$chrom.size[row.names(human_genome) == sim_rets[k]],
                                  size = 1)
      }
      
      ret_locs_same <- sim_ret_locs[sim_rets == retrogenes_by_parentals_pc[[j]]$parentalchrom]
      dists <- vector(mode = "numeric", length = length(ret_locs_same))
      if(length(ret_locs_same)){
        dists <- sapply(chr = retrogenes_by_parentals_pc[[j]]$parentalchrom, 
                        loc1 = retrogenes_by_parentals_pc[[j]]$parentalloc, 
                        X = ret_locs_same, 
                        FUN = find_dist)
        dists <- ifelse(dists>50, yes = .50, no = dists/100)
      }
      totdist <- (sum(dists) + (tot_retro-length(dists))*.5)
      avgdist <- avgdist + totdist/tot_retro
      n <- n+1
    }
  }
  sim_dists_4a[i] = avgdist/n
}
density_4a <- density(sim_dists_4a)
plot(density_4a)
hist_4a <- hist(sim_dists_4a, breaks = 15)

#test stat
totdist <- 0
n <- 0
for(i in 1:length(retrogenes_by_parentals_pc)){
  tot_retro <- sum(retrogenes_by_parentals_pc[[i]]$retroid %in% retrogenes_pc$retrodb_id)
  if(tot_retro && retrogenes_by_parentals_pc[[i]]$parentalchrom != "X"){
    retro_chroms <- retrogenes_by_parentals_pc[[i]]$retrochrom
    retro_locs <- retrogenes_by_parentals_pc[[i]]$retroloc
    
    
    ret_locs_same <- retro_locs[retro_chroms == retrogenes_by_parentals_pc[[i]]$parentalchrom]
    dists <- vector(mode = "numeric", length = length(ret_locs_same))
    if(length(ret_locs_same)){
      dists <- sapply(chr = retrogenes_by_parentals_pc[[i]]$parentalchrom, 
                      loc1 = retrogenes_by_parentals_pc[[i]]$parentalloc, 
                      X = ret_locs_same, 
                      FUN = find_dist)
      dists <- ifelse(dists>50, yes = .50, no = dists/100)
    }
    totdist_parental <- (sum(dists) + (tot_retro-length(dists))*.5)
    totdist <- totdist + totdist_parental/tot_retro
    
    n <- n+1
  }
}
ts_4b <- totdist/n


#calc p-value, make fig
pval_4c <- sum(sim_dists_4a<=ts_4b)/runs*2

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




###########################fig 4######################
######################################################

ref_dist_5a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  delta_dists_sim <- vector(mode = "numeric", length = length(retrogenes_final_pc_noX))
  for(j in 1:length(retrogenes_final_pc_noX)){
    par_chrom <- retrogenes_final_pc_noX[[j]][[1]]$parentalchrom 
    par_loc <- retrogenes_final_pc_noX[[j]][[1]]$parentalloc
    part_chroms <- retrogenes_final_pc_noX[[j]][[1]]$partnerchroms
    part_locs <- retrogenes_final_pc_noX[[j]][[1]]$partnerlocs
    
    all_chroms <- c(par_chrom, part_chroms)
    all_locs <- c(par_loc, part_locs)
    
    sim_ret_chrom <- sample(row.names(human_genome), 
                            prob = human_genome$chrom.size,
                            size = 1)
    sim_ret_loc <- sample(1:human_genome_h$chrom.size[row.names(human_genome) == sim_ret_chrom], size = 1)
    all_chrom_same <- all_locs[all_chroms == sim_ret_chrom]
    dists <- vector(mode = "numeric", length = length(all_chrom_same))
    if(length(all_chrom_same)){
      dists <- sapply(chr = sim_ret_chrom, 
                      loc1 = sim_ret_loc, 
                      X = all_chrom_same, 
                      FUN = find_dist)
      dists <- ifelse(dists>50, yes = .50, no = dists/100)
      mindist <- min(dists)
    }else{
      mindist <- .5
    }
    delta_dists_sim[j] <- mindist
  }
  ref_dist_5a[i] <- mean(delta_dists_sim)
}
density_5a <- density(ref_dist_5a)
plot(density_5a)
hist_5a <- hist(ref_dist_5a, breaks = 15)

#test stat
dists_tot <- vector(mode = "numeric", length = length(retrogenes_final_pc_noX))
for(i in 1:length(retrogenes_final_pc_noX)){
  par_chrom <- retrogenes_final_pc_noX[[i]][[1]]$parentalchrom 
  par_loc <- retrogenes_final_pc_noX[[i]][[1]]$parentalloc
  part_chroms <- retrogenes_final_pc_noX[[i]][[1]]$partnerchroms
  part_locs <- retrogenes_final_pc_noX[[i]][[1]]$partnerlocs
  
  all_chroms <- c(par_chrom, part_chroms)
  all_locs <- c(par_loc, part_locs)
  
  ret_chrom <- retrogenes_final_pc_noX[[i]][[1]]$retrochrom
  ret_loc <- retrogenes_final_pc_noX[[i]][[1]]$retroloc
  all_chrom_same <- all_locs[all_chroms == ret_chrom]
  dists <- vector(mode = "numeric", length = length(all_chrom_same))
  if(length(all_chrom_same)){
    dists <- sapply(chr = ret_chrom, 
                    loc1 = ret_loc, 
                    X = all_chrom_same, 
                    FUN = find_dist)
    dists <- ifelse(dists>50, yes = .50, no = dists/100)
    mindist <- min(dists)
  }else{
    mindist <- .5
  }
  dists_tot[i] <- mindist
}
ts_5b <- mean(dists_tot)

#calc pvalues, plot fig

pval_5c <- sum(ref_dist_5a>=ts_5b)/runs*2


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
