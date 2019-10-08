brLenChecker = function(intree) {
  for (i in (ceiling(nrow(intree$ptree)/2)+1):nrow(intree$ptree)) 
    for (j in 2:3) 
      if (intree$ptree[intree$ptree[i,j],1]-intree$ptree[i,1]<0) 
        cat(i, "\t", j, "\t", intree$ptree[intree$ptree[i,j],1], "\t", intree$ptree[i,1], "\n")
}