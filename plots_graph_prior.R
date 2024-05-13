
load("SimulationStudy/sparse_smallN/N=100_P=200/truePara.RData")
load("./SimStudy/truePara.RData")
p=200
path = ""
G = matrix(data = as.numeric( truePara$sigma != 0 ), nrow=p, ncol=p)
diag(G) = 0 
sigma = truePara$sigma

#true
true = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
diag(true) = 0  
image(true[1:20,1:20], col=c(0,1))


#empty
empty = matrix(0, nrow=p, ncol=p)
image(empty, col=c(0,1))


#partial uniform
partial_uniform = matrix(data = as.numeric( sigma != 0 ), nrow=p, ncol=p)
diag(partial_uniform) = 0  
vectorised = c(partial_uniform)
removed_last = FALSE

for (i in 1:length(vectorised)) {
  if (vectorised[i] == 1) {
    if (!removed_last) {
      vectorised[i] = 0
      removed_last = TRUE
    }
    else {
      removed_last = FALSE
    }
  }
}
  
  partial_uniform = matrix(vectorised, nrow=p, ncol=p)
  
  

image(partial_uniform[1:20, 1:20], col=c(0,1))


#partial_non_uniform
partial_non_uniform = G
partial_non_uniform[1:5, 1:5] = 0
image(partial_non_uniform[1:20, 1:20], col=c(0,1))


  


  
#noise
load(sprintf("./SimStudy/real_sim/noise_graph.RData", path))
noise_graph = G
image(noise_graph, col=c(0,1), xaxt="n", yaxt="n")
axis(1, at=c(0, 0.5, 1), labels=c("0", "100", "200"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "100", "200"))
title("Noisy graph (50%)")



#par(mfrow = c(2, 3))
layout(matrix(c(1,2,3,4,5), nrow=2, byrow=T), widths=c(1,1,1))
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2,6, byrow=T))
title.size = 2

noise_graph = G
image(noise_graph, col=c(0,1), xaxt="n", yaxt="n", cex=12)
axis(1, at=c(0, 0.5, 1), labels=c("0", "100", "200"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "100", "200"))
title("Noise", cex.main=title.size)

image(empty, col=c(0,1), xaxt="n", yaxt="n")
title("Empty", cex.main = title.size)
axis(1, at=c(0, 0.5, 1), labels=c("0", "100", "200"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "100", "200"))

#second row
image(true[1:20,1:20], col=c(0,1), xaxt="n", yaxt="n")
title("True", cex.main=title.size)
axis(1, at=c(0, 0.5, 1), labels=c("0", "10", "20"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "10", "20"))

image(partial_uniform[1:20, 1:20], col=c(0,1), xaxt="n", yaxt="n")
title("Partial uniform", cex.main=title.size)
axis(1, at=c(0, 0.5, 1), labels=c("0", "10", "20"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "10", "20"))

image(partial_non_uniform[1:20, 1:20], col=c(0,1), xaxt="n", yaxt="n")
title("Partial non-uniform", cex.main=title.size)
axis(1, at=c(0, 0.5, 1), labels=c("0", "10", "20"))
axis(2, at=c(0, 0.5, 1), labels=c("0", "10", "20"))

mtext("Overall Title", outer = TRUE, line = 1)

library(igraph)
G_p_uni = graph_from_adjacency_matrix(partial_uniform[1:15, 1:15], mode="undirected")
plot.igraph(G_p_uni)
#tkplot(G_p_uni)

G_true = graph_from_adjacency_matrix(true[1:15, 1:15], mode="undirected")
plot.igraph(G_true)

G_p_non_uni = graph_from_adjacency_matrix(partial_non_uniform[1:15, 1:15], mode="undirected")
plot.igraph(G_p_non_uni)

G_circle = graph_from_adjacency_matrix(G[1:15, 1:15], mode="undirected")
plot.igraph(G_circle)


avg = c()
for (i in 1:20) {
  load(sprintf("./SimStudy/datasets/dataset%d.RData", i))
  avg = c(avg, mean(dataset$status.train))
  
       
}
mean(avg)


