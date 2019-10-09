clusterGen = function(nfile, lfile, clustcol) {
  nodes = read.csv(nfile, sep = "\t", header = T)
  links = read.csv(lfile, sep = "\t", header = T)
  lty = c("dashed", "solid") # Known transmission is solid line, novel is dashed.
  net = graph_from_data_frame(d = links, vertices = nodes, directed=T)
  # shape = c("circle", "triangle", "square")
  # V(net)$shape = shape[V(net)$Transmission.Information + 1]
  V(net)$color = clustcol[V(net)$Sampled + 1]
  E(net)$edge.color = "gray80"
  # E(net)$edge.lty = lty[E(net)$Known + 1]
  E(net)$arrow.size = 0.75
  plot(net, edge.width = 1.5, edge.lty = lty[E(net)$Known + 1], vertex.size = 15, vertex.label=V(net)$Patient.ID, vertex.label.family = 2, vertex.label.font=1, vertex.label.color="black", vertex.label.cex=.7, vertex.label.degree = 0, vertex.label.dist = 4)
  return(net)
}