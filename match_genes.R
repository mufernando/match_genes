all = read.csv ("all_genes.csv", header=T)
interest = read.csv ("genes_of_interest.csv", header=T)

match_genes = function (all, interest, opt = 10^4, graph = F) {
    #Initialize the data.frame that is going to be returned.
    match = data.frame(matrix(ncol=3, nrow=0))
    #The function is gonna use the median to calculate gene distances.
    all$md = apply(all[,c(3,4)], 1, function (x) median(x))
    #Same as above, but now for interest data.frame.
    interest$md = apply(interest[,c(3,4)], 1, function (x) median(x))
    #Here, we start the first loop; we are going to define the control genes one by one.
    for (i in 1:nrow(interest)) {
        #Subset the data.frame all containing only genes in the same chromosome as the interest gene.
        data = subset(all, all[,2] == interest[i, 2])
        #Remove genes that were a match already.
        data = subset(data, !(data[,1] %in% match[,2]))
        #Just in case the data.frame isn't ordered by position in the chromosome.
        data = data[order(data$md),]
        #Calculate distance betwen gene of interest median position and every gene availiable.
        data$dists = abs(interest[i,6] - data[,6])
        #Creates a data.frame with the interest gene ID, the ID of the gene closer to the optimal distance and the distance.
        m = data.frame(interest[i,1], data[which.min(abs(data$dists-opt)),][c(1,7)])
        #Binds the row created before to the existing match data.frame.
        match = rbind(match, m)
    }
    #Assign names to the columns of match.
    colnames(match) = c("interest_ID", "match_ID", "distance")
    #Checks if the parameter graph is True.
    if (graph == T) {
        #Starts a graphics device X.
        x11()
        #Plot the distance between genes in a boxplot.
        boxplot(match$distance)
        #Put a title to the plot.
        title("Distances between each pair gene of interest/control")
        #Traces a line in the optimal distance.
        abline(h=opt, col="red")
        #Closes the graphic device.
        dev.off()
    }
    #Return the match data.frame and the median distance gene of interest/control gene
    return (list(match, median(match$distance)))
}

shuffle_match = function (all, interest, opt = 10^4, times = 100) {
    #Creates an empty list.
    results = list()
    #Loop to run match_genes several times.
    for (i in 1:times) {
        #Shuffle the interest data.frame
        interest = interest[sample(nrow(interest)),]
        #Appends the results from match_genes to a pre-existing list called results
        results = append(results, list(match_genes(all, interest, opt)))
    }
    #Gets only the one result which got the median distance closest to the optimal.
    results = results[[which.min(abs(sapply(results, function(x) { x[[2]] })-opt))]][[1]]
    #Return results
    return (results)
}

shuffle_match(all, interest)
