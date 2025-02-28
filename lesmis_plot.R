# Weighted Plot -----------------------------------------------------------

library(RSpectra)
library(igraph)
library(RColorBrewer)

real_data = read.graph("dataset/Realdataset/lesmis.gml", format = "gml")
A = as.matrix(as_adjacency_matrix(real_data, attr = 'value'))
n = nrow(A)

par(mar=c(0.2,0.2,0.2,0.2))

lesmis_weighted <- graph_from_adjacency_matrix(A, weighted = T, mode = "undirected")
lay <- layout_with_fr(lesmis_weighted)


set.seed(1024)

m = 5
score_5_label <- SCORE(A, m)$labels
rsc_5_label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster

# Clustering result of SCORE
deg <- degree(lesmis_weighted, mode="all")
ncolor <-  brewer.pal(m, "Set1")
V(lesmis_weighted)$size <- 5
V(lesmis_weighted)$color <- ncolor[score_5_label]
V(lesmis_weighted)$frame.color <- NA
V(lesmis_weighted)$label <- ""
plot(lesmis_weighted,layout=lay)

# Clustering result of RSC
new_color <- c("#E41A1C",  "#984EA3", "#FF7F00", "#377EB8", "#4DAF4A")
V(lesmis_weighted)$color <- new_color[rsc_5_label]
plot(lesmis_weighted,layout=lay)


m = 6
score_6_label <- SCORE(A, m)$labels
rsc_6_label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster

# Clustering result of SCORE
deg <- degree(lesmis_weighted, mode="all")
ncolor <-  brewer.pal(m, "Set1")
V(lesmis_weighted)$size <- 5
V(lesmis_weighted)$color <- ncolor[score_6_label]
V(lesmis_weighted)$frame.color <- NA
V(lesmis_weighted)$label <- ""
plot(lesmis_weighted, layout=lay)

# Clustering result of RSC
new_color <- c("#E41A1C", "#FFFF33", "#FF7F00","#4DAF4A",  "#984EA3", "#377EB8")
V(lesmis_weighted)$color <- new_color[rsc_6_label]
plot(lesmis_weighted, layout=lay)



m = 7
score_7_label <- SCORE(A, m)$labels
rsc_7_label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster

# Clustering result of SCORE
deg <- degree(lesmis_weighted, mode="all")
ncolor <-  brewer.pal(m, "Set1")
V(lesmis_weighted)$size <- 5
V(lesmis_weighted)$color <- ncolor[score_7_label]
V(lesmis_weighted)$frame.color <- NA
V(lesmis_weighted)$label <- ""

plot(lesmis_weighted,layout=lay)

# Clustering result of RSC
new_color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",  "#FFFF33", "#A65628")
V(lesmis_weighted)$color <- new_color[rsc_7_label]
plot(lesmis_weighted, layout=lay)




# Unweighted Plot ---------------------------------------------------------

real_data <- read.graph("dataset/Realdataset/lesmis.gml", format = "gml")
A <- as.matrix(as_adjacency_matrix(real_data, attr = 'value'))
n <- nrow(A)
A[which(A > 1)] <- 1

par(mar=c(0.2, 0.2, 0.2, 0.2))
lesmis_unweighted <- graph_from_adjacency_matrix(A, weighted = NULL, mode = "undirected")


set.seed(1024)

m <- 3
score_label <- SCORE(A, m)$labels
rsc_label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster

# Clustering result of SCORE
deg <- degree(lesmis_unweighted, mode="all")
ncolor <-  brewer.pal(m, "Set1")
V(lesmis_unweighted)$size <- 5
V(lesmis_unweighted)$color <- ncolor[score_label]
V(lesmis_unweighted)$frame.color <- NA
V(lesmis_unweighted)$label <- ""
plot(lesmis_unweighted, layout=lay)

# Clustering result of RSC
rsc_color <- c("#377EB8", "#E41A1C", "#4DAF4A")
V(lesmis_unweighted)$color <- rsc_color[rsc_label]
plot(lesmis_unweighted, layout=lay)





















