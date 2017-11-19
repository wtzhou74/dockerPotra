print("Hello Potra")
PoTRA.Ftest <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)
  
  Ftest<-c()
  E.normal<-c()
  E.case<-c()
  length.pathway<-c()
  
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  for (x in 1:length(Pathway.database)){
    print(x)
    p0 <-Pathway.database[[x]]
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)    
    length.intersect<-length(intersect(unlist(nodelist),unlist(genelist)))
    
    length.pathway[x]<-length(nodelist)
    
    #collect expression data of genes for a specific pathway across normal and tumor samples.
    
    path<-data.frame(matrix(200,length.intersect,(Num.sample.normal+Num.sample.case)))
    a<- c()
    
    for (j in 1:length.intersect){
      a[j]<-intersect(unlist(nodelist),unlist(genelist))[j]  
      
      path[j,]<-mydata[which(genelist==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
    }
    
    ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
    
    cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
    
    cor.normal<-as.matrix(cor.normal) 
    
    cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
    
    cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
    cor.normal.adj[ is.na(cor.normal.adj)] <- 0
    diag(cor.normal.adj) <- 0
    graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
    E.normal[x]<-length(E(graph.normal))
    
    PR.normal<-page.rank(graph.normal,direct=FALSE)
    
    ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
    
    cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
    
    cor.case<-as.matrix(cor.case) 
    cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
    
    cor.case.adj[ cor.case.adj > 0.05 ] <- 0
    cor.case.adj[ is.na(cor.case.adj)] <- 0
    diag(cor.case.adj) <- 0
    graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
    E.case[x]<-length(E(graph.case))
    
    PR.case<-page.rank(graph.case,direct=FALSE)
    
    Ftest[x]<-var.test(PR.normal$vector,PR.case$vector)$p.value
  }
  return(list(p.value=Ftest,TheNumberOfEdges.normal=E.normal,TheNumberOfEdges.case=E.case,LengthOfPathway=length.pathway))
}









PoTRA <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database, PR.quantile) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)
  
  Fishertest<-c()
  TheNumOfHubGene.normal<-c()
  TheNumOfHubGene.case<-c()
  E.normal<-c()
  E.case<-c()
  length.pathway<-c()
  
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  for (x in 1:length(Pathway.database)){
    print(x)
    p0 <-Pathway.database[[x]]
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)    
    length.intersect<-length(intersect(unlist(nodelist),unlist(genelist)))
    
    length.pathway[x]<-length(nodelist)
    
    if (length.intersect<5){
      next
    }else{
      
      #collect expression data of genes for a specific pathway across normal and tumor samples.
      
      path<-data.frame(matrix(0,length.intersect,(Num.sample.normal+Num.sample.case)))
      a<- c()
      
      for (j in 1:length.intersect){
        a[j]<-intersect(unlist(nodelist),unlist(genelist))[j]  
        
        path[j,]<-mydata[which(genelist==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
      }
      
      ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
      
      cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
      
      cor.normal<-as.matrix(cor.normal) 
      
      cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
      
      cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
      cor.normal.adj[ is.na(cor.normal.adj)] <- 0
      diag(cor.normal.adj) <- 0
      graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
      E.normal[x]<-length(E(graph.normal))
      
      PR.normal<-page.rank(graph.normal,direct=FALSE)$vector
      
      ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
      
      cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
      
      cor.case<-as.matrix(cor.case) 
      cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
      
      cor.case.adj[ cor.case.adj > 0.05 ] <- 0
      cor.case.adj[ is.na(cor.case.adj)] <- 0
      diag(cor.case.adj) <- 0
      graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
      E.case[x]<-length(E(graph.case))
      
      PR.case<-page.rank(graph.case,direct=FALSE)$vector
      
      matrix.HCC<- matrix("NA",2*length.intersect,2)
      rownames(matrix.HCC)<-as.character(c(PR.normal,PR.case))
      colnames(matrix.HCC)<-c("Disease_status","PageRank")
      
      matrix.HCC[,1]<-c(rep("Normal",length.intersect), rep("Cancer",length.intersect))
      
      loc.largePR<-which(as.numeric(rownames(matrix.HCC))>=quantile(PR.normal,PR.quantile))
      loc.smallPR<-which(as.numeric(rownames(matrix.HCC))<quantile(PR.normal,PR.quantile))
      
      matrix.HCC[loc.largePR,2]<-"large_PageRank"
      matrix.HCC[loc.smallPR,2]<-"small_PageRank"
      
      table.HCC<-list(1,2)
      names(table.HCC)<-c("Disease_status","PageRank")
      
      table.HCC$Disease_status<-matrix("NA",2*length.intersect,2)
      table.HCC$PageRank<-matrix("NA",2*length.intersect,2)
      
      table.HCC$Disease_status<-matrix.HCC[,1]
      table.HCC$PageRank<-matrix.HCC[,2]
      
      cont.HCC<-table(table.HCC$Disease_status,table.HCC$PageRank)
      TheNumOfHubGene.normal[x]<-cont.HCC[2]
      TheNumOfHubGene.case[x]<-cont.HCC[1]
      
      Fishertest[x]<-fisher.test(cont.HCC)$p.value
    }
  }
  return(list(p.value=Fishertest,TheNumberOfHubGenes.normal=TheNumOfHubGene.normal,TheNumOfHubGene.case=TheNumOfHubGene.case,TheNumberOfEdges.normal=E.normal,TheNumberOfEdges.case=E.case))
}









overlayplot <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database) {
  
  require(BiocGenerics)
  require(graph)
  require(graphite)
  require(igraph)
  
  
  length.pathway<-c()
  
  humanReactome <- pathways("hsapiens", "reactome")
  humanBiocarta <- pathways("hsapiens", "biocarta")
  humanKEGG <- pathways("hsapiens", "kegg")
  
  plot.multi.dens <- function(s)
  {
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s))
    {
      junk.x = c(junk.x, density(s[[i]])$x)
      junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
    for(i in 1:length(s))
    {
      lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
    }
  }
  
  for (x in Pathway.database){
    p0 <-Pathway.database[[x]]
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)    
    length.intersect<-length(intersect(unlist(nodelist),unlist(genelist)))
    
    length.pathway[x]<-length(nodelist)
    
    #collect expression data of genes for a specific pathway across normal and tumor samples.
    
    path<-data.frame(matrix(200,length.intersect,(Num.sample.normal+Num.sample.case)))
    a<- c()
    
    for (j in 1:length.intersect){
      a[j]<-intersect(unlist(nodelist),unlist(genelist))[j]  
      
      path[j,]<-mydata[which(genelist==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
    }
    
    ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
    
    cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
    
    cor.normal<-as.matrix(cor.normal) 
    
    cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
    
    cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
    cor.normal.adj[ is.na(cor.normal.adj)] <- 0
    diag(cor.normal.adj) <- 0
    graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
    
    PR.normal<-page.rank(graph.normal,direct=FALSE)
    
    ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
    
    cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
    
    cor.case<-as.matrix(cor.case) 
    cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
    
    cor.case.adj[ cor.case.adj > 0.05 ] <- 0
    cor.case.adj[ is.na(cor.case.adj)] <- 0
    diag(cor.case.adj) <- 0
    graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
    
    PR.case<-page.rank(graph.case,direct=FALSE)
    
    PoTRA.plot<-plot.multi.dens(list(PR.normal,PR.case))
  }
  return(PoTRA.plot)
}









