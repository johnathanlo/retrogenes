
load("DM-final-data.RData") #final dataset made 7/14/21
set.seed(6)
library(dplyr)
library(rentrez)
library(chromPlot)
library(GenomicFeatures)
library(ggplot2)
library(reshape2)
#library(HTSanalyzeR)
library(org.Dm.eg.db)
library(MareyMap)

##distance function
data("Drosophila_melanogaster")
dm_map <- Drosophila_melanogaster
find_dist <- function(chr1, chr2, loc1, loc2){
  if(strsplit(chr1, split = "")[[1]][[1]] == strsplit(chr2, split = "")[[1]][[1]]){
    positions1 <- physicalPositions(dm_map[[chr1]])
    distances1 <- geneticDistances(dm_map[[chr1]])
    positions2 <- physicalPositions(dm_map[[chr2]])
    distances2 <- geneticDistances(dm_map[[chr2]])
    loc1_match <- distances1[which.min(abs(loc1 - positions1))] 
    loc2_match <- distances2[which.min(abs(loc2 - positions2))]
    rec_dist <- abs(loc1_match - loc2_match)/100
    return(rec_dist)
  }else{
    return(.5)
  }
  
}

#total number of partner genes
all_partners <- c()
for(i in 1:length(retrogenes_final)){
  partners_i <- retrogenes_final[[i]][[1]]$partnergenes
  all_partners <- c(all_partners, partners_i)
}
length(unique(all_partners))


#########sup fig 1 (ootX)
#######################################################
#build reference distribution for set of all retrogenes
ref_dist_test_ootX <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  sim_parents <- sample(c("X","Y", "2L", "2R", "3L", "3R", "4"), 
                        size = length(retrogenes_final_pc), 
                        prob = dm_genes_per_chrom$gene.num,
                        replace = TRUE)
  ref_dist_test_ootX[i] <- sum(sim_parents == "X")
}
density_test_ootX <- density(ref_dist_test_ootX)
hist_test_ootX <- hist(ref_dist_test_ootX)
#take a look at reference distributions
plot(density_test_ootX)

#ootx test statistic
totalX_ootX <- 0
for(i in 1:length(retrogenes_final_pc)){
  if(retrogenes_final_pc[[i]][[1]]$parentalchrom == "X"){
    totalX_ootX <- totalX_ootX + 1
  }
}
p_val_ootX <- sum(ref_dist_test_ootX >= totalX_ootX)/runs * 2

##plot
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


#########fig 2
#######################################################

sim_dists_3a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  totdist <- 0
  for(j in 1:length(retrogenes_final_pc_noX)){
    sim_ret <- sample(row.names(dm_genome), 
                      prob = dm_genome$chrom.size, 
                      size = 1, 
                      replace = TRUE)
    sim_ret_loc <- sample(1:dm_genome$chrom.size[row.names(dm_genome) == sim_ret],
                          size = 1)
    
    dist <- find_dist(sim_ret, retrogenes_final_pc_noX[[j]][[1]]$parentalchrom,
                      sim_ret_loc, retrogenes_final_pc_noX[[j]][[1]]$parentalloc)
    totdist <- totdist+dist
  }
  sim_dists_3a[[i]] <- totdist/length(retrogenes_final_pc_noX)
}
density_3a <- density(sim_dists_3a)
plot(density_3a)
hist_3a <- hist(sim_dists_3a, breaks = 15)

dists_3b <- 0
dists_close <- vector(mode = "logical", length = length(retrogenes_final_pc_noX))
for(i in 1:length(retrogenes_final_pc_noX)){
  retro_chrom <- retrogenes_final_pc_noX[[i]][[1]]$retrochrom
  parent_chrom <- retrogenes_final_pc_noX[[i]][[1]]$parentalchrom
  retro_loc <- retrogenes_final_pc_noX[[i]][[1]]$retroloc
  parent_loc <- retrogenes_final_pc_noX[[i]][[1]]$parentalloc
  
  dist <- find_dist(retro_chrom, parent_chrom, retro_loc, parent_loc)
  dists_3b <- dists_3b + dist
  if(dist < .2){
    dists_close[i] <- T
  }
}
ts_3b <- dists_3b/length(retrogenes_final_pc_noX)

pval_3c <- sum(sim_dists_3a<=ts_3b)/runs*2


##plot
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

#########fig 3
#######################################################

parental_noX <- 0
for(i in 1:length(retrogenes_by_parentals)){
  if(retrogenes_by_parentals[[i]]$parentalchrom != "X" && 
     sum(retrogenes_by_parentals[[i]]$retroid %in% retrogenes_pc$retrodb_id) == length(retrogenes_by_parentals[[i]]$retroid)){
    parental_noX <- parental_noX + 1
  }
}

sim_dists_4a <- vector(mode = "numeric", length = runs)
for(i in 1:runs){
  totdist <- 0
  n <- 0
  for(j in 1:length(retrogenes_by_parentals)){
    tot_retro <- sum(retrogenes_by_parentals[[j]]$retroid %in% retrogenes_pc$retrodb_id)
    if(retrogenes_by_parentals[[j]]$parentalchrom != "X" && tot_retro){
      sim_rets <- sample(row.names(dm_genome), 
                         prob = dm_genome$chrom.size, 
                         size = tot_retro, 
                         replace = TRUE)
      sim_ret_locs <- vector(mode = "numeric", length=tot_retro)
      parent_chrom <- retrogenes_by_parentals[[j]]$parentalchrom
      parent_loc <- retrogenes_by_parentals[[j]]$parentalloc
      
      for(k in 1:length(sim_rets)){
        sim_ret_locs[k] <- sample(1:dm_genome$chrom.size[row.names(dm_genome) == sim_rets[k]],
                                  size = 1)
      }
      
      sample_dist <- 0
      for(k in 1:tot_retro){
        dist <- find_dist(sim_rets[k], parent_chrom,
                          sim_ret_locs[k], parent_loc)
        sample_dist <- sample_dist + dist
      }
      totdist <- totdist + sample_dist/tot_retro
      n<-n+1
    }
  }
  sim_dists_4a[i] <- totdist/n
}
density_4a <- density(sim_dists_4a)
plot(density_4a)
hist_4a <- hist(sim_dists_4a, breaks = 15)

totdist_4b<- 0
k = 0
for(i in 1:length(retrogenes_by_parentals)){
  tot_retro <- sum(retrogenes_by_parentals[[i]]$retroid %in% retrogenes_pc$retrodb_id)
  if(length(retrogenes_by_parentals[[i]]$retrochrom) && retrogenes_by_parentals[[i]]$parentalchrom != "X" 
     && tot_retro>0){
    retro_chroms <- retrogenes_by_parentals[[i]]$retrochrom
    retro_ids <- retrogenes_by_parentals[[i]]$retroid
    retro_locs <- retrogenes_by_parentals[[i]]$retroloc
    parent_chrom <- retrogenes_by_parentals[[i]]$parentalchrom
    parent_loc <- retrogenes_by_parentals[[i]]$parentalloc
    
    avg_dist <- 0
    for(j in 1:tot_retro){
      if(retro_ids[j] %in% retrogenes_pc$retrodb_id){
        dist <- find_dist(retro_chroms[j], parent_chrom,
                          retro_locs[j], parent_loc)
        avg_dist <- avg_dist + dist/tot_retro
      }
    }
    totdist_4b <- totdist_4b + avg_dist
    k = k+1
  }
}

ts_4b <- totdist_4b/k

pval_4c <- sum(sim_dists_4a<=ts_4b)/runs*2

##plot
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

#########fig 4
#######################################################

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
    
    sim_ret_chrom <- sample(row.names(dm_genome), 
                            prob = dm_genome$chrom.size,
                            size = 1)
    sim_ret_loc <- sample(1:dm_genome_h$chrom.size[row.names(dm_genome) == sim_ret_chrom], size = 1)
    
    all_dists <- vector(mode = "numeric" , length = length(all_chroms))
    for(k in 1: length(all_chroms)){
      all_dists[k] <- find_dist(chr1 = sim_ret_chrom, chr2 = all_chroms[k], 
                                loc1 = sim_ret_loc, loc2 = all_locs[k])
    } 
    
    delta_dists_sim[j] <- min(all_dists, na.rm = T)
  }
  ref_dist_5a[i] <- mean(delta_dists_sim)
}
density_5a <- density(ref_dist_5a)
plot(density_5a)
hist_5a <- hist(ref_dist_5a, breaks = 15)


dists <- vector(mode = "numeric", length = length(retrogenes_final_pc_noX))
for(i in 1:length(retrogenes_final_pc_noX)){
  par_chrom <- retrogenes_final_pc_noX[[i]][[1]]$parentalchrom 
  par_loc <- retrogenes_final_pc_noX[[i]][[1]]$parentalloc
  part_chroms <- retrogenes_final_pc_noX[[i]][[1]]$partnerchroms
  part_locs <- retrogenes_final_pc_noX[[i]][[1]]$partnerlocs
  
  all_chroms <- c(par_chrom, part_chroms)
  all_locs <- c(par_loc, part_locs)
  
  ret_chrom <- retrogenes_final_pc_noX[[i]][[1]]$retrochrom
  ret_loc <- retrogenes_final_pc_noX[[i]][[1]]$retroloc
  all_dists <- vector(mode = "numeric" , length = length(all_chroms))
  for(j in 1:length(all_chroms)){
    all_dists[j] <- find_dist(chr1 = all_chroms[j], chr2 = ret_chrom,
                              loc1 = all_locs[j], loc2 = ret_loc)
  }
  
  dists[i] <- min(all_dists, na.rm = T)
}
ts_5b <- mean(dists, na.rm = T)

pval_5c <- sum(ref_dist_5a>=ts_5b)/runs*2

##plot
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