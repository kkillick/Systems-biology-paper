

#Using igraph package to generate random networks

#Generate 1000 random networks and calculate the average Clustering Coefficient (transitivity) and Average shortest path Length


#Real script

#Set up an empty array with 1000 rows and 3 columns
#Call this array graphRes
#Use for loop to repeat 1000 times
  #Generate a random network using the Erdos-Renyi model that consists of 435 nodes and 1025 edges (the same as is in our bov24_vs_BCG24hr (FDR0.01) network)
	#Use transivity function to calculate clustering coefficient of this network and place it in column 1 of the empty array graphRes
	#Use average.path.length function (part of the shortest path function) to calculate the average shortest path legth of this network and place it in column 2 of the empty array graphRes
#Write graphRes out as a .txt file





graphRes <- array (,c(1000,3))
  for (i in c(1:1000)) {
    g <- erdos.renyi.game(435, 1025, type=c("gnm"))
    graphRes[i,1] <- transitivity(g)
    graphRes[i,2] <- average.path.length(g, directed=FALSE, unconnected=TRUE)
    print (transitivity(g))
    print (average.path.length(g, directed=FALSE, unconnected=TRUE))
  }

graphRes
 CC <- mean(graphRes[,1])
 ASPL <- mean(graphRes[,2])


write.table(graphRes, file = "random_network.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
