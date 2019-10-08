# Assuming we've already made sac_p, the phylogenetic object and it doesn't have negative branch lengths.

autoSummary = function(sac_p, num_sampled, startpi, updatepi, startr, gshape, gscale, sshape, sscale) {
  
  sac = inferTTree(sac_p, w.shape=gshape, w.scale=gscale, ws.shape=sshape, ws.scale=sscale, mcmcIterations=1000, thinning=1, startNeg=680/365, startOff.r=startr, startOff.p=0.5, startPi=startpi, updateNeg=TRUE, updateOff.r=TRUE, updateOff.p=FALSE, updatePi=updatepi, startCTree=NA, updateTTree=TRUE, optiStart=TRUE, dateT=Inf)	# MCMC of epidemiological and infectious parameters. 
  summary_data <- data.frame(a=c("Start Pi", "Update Pi", "Inferred Pi", "Total Transmissions", "Infected >= 1989", "Infected >= 1989 and Sampled", "Infected >= 1989 (IN) and Sampled", "Infected >= 1989 (ON) and Sampled", "Infected < 1989 and Sampled", "Infected ><1989 (IN) and Sampled", "Historical", "1 Secondary", "Multi Secondary"))
  summary_data[1,2] = startpi
  summary_data[2,2] = updatepi
  sac_last = sac[[length(sac)]]	# Epidemiological and infectious properties of the last iteration. 
  sac_ctree = sac_last$ctree 		# The (MCC?) combined tree associated with the last iteration
  summary_data[3,2] = sac_last$pi 	# The proportion of all patients from root to all tips that are sampled.
  sac_ttree = extractTTree(sac_ctree)
  summary_data[4,2] = length(sac_ttree$ttree[,1])		# The number of total individuals that were infected. 
  summary_data[5,2] = length(which(sac_ttree$ttree[,1]>1989))	# The number of individuals that were infected after 1989 i.e.: known to SAC. 
  summary_data[6,2] = length(which(sac_ttree$ttree[,1]>1989&!is.na(sac_ttree$ttree[,2])))	# The number of individuals as above and sampled.
  summary_data[7,2] = length(which(sac_ttree$ttree[,1]>1989&!is.na(sac_ttree$ttree[,2])&sac_ttree$ttree[,3]<num_sampled))	# The number of individuals as above and infected by another known SAC patient in the tree.
  summary_data[8,2] = length(which(sac_ttree$ttree[,1]>1989&!is.na(sac_ttree$ttree[,2])&sac_ttree$ttree[,3]>num_sampled))	# Opposite of the above.
  summary_data[9,2] = length(which(sac_ttree$ttree[,1]<1989&!is.na(sac_ttree$ttree[,2])))		# The number of individuals infected before 1989 and sampled. 
  summary_data[10,2] = length(which(sac_ttree$ttree[,1]<1989&!is.na(sac_ttree$ttree[,2])&sac_ttree$ttree[,3]<num_sampled))	# The number of individuals as above and infected by another known SAC patient in the tree.
  summary_data[11,2] = length(which(sac_ttree$ttree[,1]<1989&is.na(sac_ttree$ttree[,2])))	# The number of unsampled individuals that were infected before the SAC started tracking patients.
  hist_multi_trans = data.frame(table(table(sac_ttree$ttree[which(sac_ttree$ttree[,1]>1989),3])))
  summary_data[12,2] = hist_multi_trans[1,2]
  summary_data[13,2] = 0
  for (i in 2:length(hist_multi_trans)) {
    summary_data[13,2] = summary_data[13,2] + hist_multi_trans[i,2]
  }
  colnames(summary_data) = c("Summary", "Data")
  cat("\n")
  for (i in 1:13) {
    cat(summary_data[i,2], "\n")
  }
  return(sac_ctree)
}