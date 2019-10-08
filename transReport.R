transReport = function(c, num_samp) {
  lent = extractTTree(c)
  report = lent$ttree
  nam = lent$nam
  library(lubridate)
  na_vect = rep("NA", length(report[,1]) - num_samp)
  id_vect = c(nam, na_vect)
  tmp = as.data.frame(report)
  tmp[,1] = as.numeric(as.character(tmp[,1]))
  tmp[,2] = as.numeric(as.character(tmp[,2]))
  report2 = data.frame(matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp) + 2))
  for (r in 1:length(tmp[,1])) {
    report2[r,1] = r
    report2[r,2] = id_vect[r]
    report2[r,3] = format(date_decimal(tmp[r,1]), "%Y-%m-%d")
    report2[r,4] = format(date_decimal(tmp[r,2]), "%Y-%m-%d")
    if (0 < tmp[r,3] && tmp[r,3] <= num_samp) {
      id = tmp[r,3]
      report2[r,5] = nam[id]
    } else {
      report2[r,5] = tmp[r,3]
    }
  }
  file = paste(getwd(),"/predicted_transmission_table.tsv", sep = "")
  write.table(report2, file=file, quote=FALSE, sep='\t', row.names = FALSE, col.names = c("Patient No.","SAC Patient ID","Infection Date (Pred.)","Sampling Date","No. or ID of Infector"))
  return(lent)
}
