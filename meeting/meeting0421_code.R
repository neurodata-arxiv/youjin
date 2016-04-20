library(igraph)

## read a data
elegans <- read.graph("data/c.elegans.herm_pharynx_1.graphml", format = "graphml")
soma_pos <- V(elegans)$soma_pos 

png("figure/somapos.png")
igraph.options(vertex.size=4.5, vertex.label=NA,
               edge.arrow.size=0.5)
V(elegans)[soma_pos <= 0.1]$color <- "lightblue1"
V(elegans)[(soma_pos > 0.1) & (soma_pos <= 0.2)]$color <- "skyblue"
V(elegans)[(soma_pos > 0.2) & (soma_pos <= 0.5)]$color <- "skyblue3"
V(elegans)[(soma_pos > 0.5) & (soma_pos <= 0.7)]$color <- "dodgerblue"
V(elegans)[(soma_pos > 0.7)]$color <- "dodgerblue4"
E(elegans)$width <- 0.3
E(elegans)$color <- "ivory3"
plot(elegans,layout=layout.fruchterman.reingold,
     main = "Neuronal Network")
legend("topright", c("<=0.1", "<=0.2", "<=0.5", "<=0.7", ">0.7"),
       col = c("lightblue1", "skyblue", "skyblue3", "dodgerblue", "dodgerblue4"),
       pch = 20)
dev.off()