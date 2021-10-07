## ---------------------------
##
## Script name: Species-level Kelp-Forest Food Web (extant species only)
##
## Purpose of script: Generates 3 versions of the kelp-forest food web, aggregated to species level. Also generates igraph objects that can be used in plotting and caculates basic metrics
##
## Author: Dr. Dana Morton
##
## Date Created: September 29, 2021
##
## Copyright (c) Dr. Dana Morton

## Email: dnmorton@colby.edu
##
## ---------------------------
##
## Notes: 
## Two input files are the nodes and links list from Dryad:https://datadryad.org/stash/dataset/doi:10.25349/D9JG70
## See column descriptor file for metadata: "0_Column_Descriptors.csv"
## See Morton et al. 2021 https://www.nature.com/articles/s41597-021-00880-4 for description of source dataset
##
## ---------------------------

## set working directory for Mac and PC
setwd("~/YOUR directory/")     
setwd("C:/Users/YOUR director/")    

#load packages needed
library(plyr)
library(dplyr)
library(igraph)
library(NetIndices)
library(reshape2)
library(tidyr)
library(NetSwan)


#load node and edge lists (life-stage resolution), downloaded from Dryad (see Notes above)
kelp.nodes.all <- read.csv("1_Nodes.csv", header = T) 
kelp.edge.all <- read.csv("2_Links.csv", header = T)

#Filter node and edge lists to create appropriate subwebs:
    #filter node list to remove extinct species, coded with "E" in the Presence.code column
    kelp.node.extant <- subset.data.frame(kelp.nodes.all, Presence.code != "E")
    
    #create free-living only and parasite-only node list
    free.nodes <- subset.data.frame(kelp.node.extant, Type != "symbiont")
    ps.nodes <- subset.data.frame(kelp.node.extant, Type == "symbiont")
    
    nodes <- kelp.node.extant
    kelp.edge <- kelp.edge.all
    
    #create node and edge lists for different versions of the web:
      #subset edge list to include same sp as node list (ie. remove extinct links)
      kelp.edge <- subset.data.frame(kelp.edge, consumers %in% nodes$Node.ID)
      kelp.edge <- subset.data.frame(kelp.edge, resources %in% nodes$Node.ID)
    
    #remove concomittant links, see metadata file ("0_Column_Descriptors.csv") for codes
      kelp.edge.no.conc <- subset.data.frame(kelp.edge, Consumer.Interaction.Code != 14) 
    
    #subset to create free-living edge list
      kelp.free.edge <- subset.data.frame(kelp.edge, consumers %in% free.nodes$Node.ID)
      kelp.free.edge <- subset.data.frame(kelp.free.edge, resources %in% free.nodes$Node.ID)



#Create web aggregated to species with predator-prey, parasite-host, and predator-parasite (concomitant links) 
    ###rename node and edge list to versions you will work with
    nodes <- kelp.node.extant
    kelp.edge.1 <- kelp.edge
    
    #calculate trophic level
    #convert edge list to adj matrix (use consuerSP and resourceSP columns as identifiers)
    edge.only.sp <- kelp.edge.1[,c(4,9)] 
    edge.only.sp <- unique(edge.only.sp)
    g.sp = graph.data.frame(edge.only.sp)
    kelp.mat.sp <- get.adjacency(g.sp)  #get adjacency matrix from graph object
    kelp.mat.sp <- as.matrix(kelp.mat.sp)  #convert to numeric matrix
    
    #calculate trophic level and add to node list
    kelp.matT.sp <- t(kelp.mat.sp)
    tl.sp <- TrophInd(kelp.matT.sp)  # calculate trophic level, needs numeric matrix
    tl.sp$sp.ID <- rownames(tl.sp)
    tl.sp$sp.ID <- as.numeric(tl.sp$sp.ID)
    
    #join TL to node list, this also drops any of nodes that aren't in the edge list, as only those have trophic levels
    nodes.tl.sp <- right_join(nodes, tl.sp, by = "sp.ID")  
    
    ##create species-level node list
    collapse_unique <- function(x) {
      paste(unique(x), collapse = ",")
    }
    all.species.nodes <- ddply(nodes.tl.sp,.(sp.ID),colwise(collapse_unique))  
    
    #reorder edgelist so resourcespID is first column, consumer SpID is second, for reading into igraph
    #Note: this edge list is not yet aggregated to species. creation of the igraph object will remove duplicate links. 
    
    all.kelp.edge.sp <- kelp.edge.1[,c(9,4,1:3,5:8,10:17)]
    
    #create igraph object, for graphing and analysis.
    kelp.g.sp.parasites.w.conc<- graph.data.frame(all.kelp.edge.sp,
                                                  directed = T,
                                                 vertices = all.species.nodes)
    

#Create web aggregated to species with predator-prey, parasite-host, but NO predator-parasite (concomitant links) 
    ###rename node and edge list to versions you will work with
    nodes <- nodes
    kelp.edge.1 <- kelp.edge.no.conc
    
    #calculating trophic level
    #convert edge list to adj matrix (use consuerSP and resourceSP columns as identifiers)
    edge.only.sp <- kelp.edge.1[,c(4,9)] 
    edge.only.sp <- unique(edge.only.sp)
    g.sp = graph.data.frame(edge.only.sp)
    kelp.mat.sp <- get.adjacency(g.sp)  #get adjacency matrix from graph object
    kelp.mat.sp <- as.matrix(kelp.mat.sp)  #convert to numeric matrix
    
    #calculate trophic level and add to node list
    kelp.matT.sp <- t(kelp.mat.sp)
    tl.sp <- TrophInd(kelp.matT.sp)  # calculate trophic level, needs numeric matrix
    tl.sp$sp.ID <- rownames(tl.sp)
    tl.sp$sp.ID <- as.numeric(tl.sp$sp.ID)
    
    #join TL to node list, this also drops any of nodes that aren't in the edge list, as only those have trophic levels
    nodes.tl.sp <- right_join(nodes, tl.sp, by = "sp.ID")  
    
    ##create sp level node list
    collapse_unique <- function(x) {
      paste(unique(x), collapse = ",")
    }
    species.nodes <- ddply(nodes.tl.sp,.(sp.ID),colwise(collapse_unique))  
    
    #reorder edgelist so resourcespID is first column, consumer SpID is second, for reading into igraph
    #Note: this edge list is not yet aggregated to species. creation of the igraph object will remove duplicate links. 
    kelp.edge.sp.no.conc <- kelp.edge.1[,c(9,4,1:3,5:8,10:17)]
    
    #create igraph object, for graphing and analysis
    kelp.g.sp.parasites.NO.conc<- graph.data.frame(kelp.edge.sp.no.conc,
                                                  directed = T,
                                                  vertices = species.nodes)



### #Create web aggregated to species with predator-prey links only (Free-living Species Only)
    ###rename node and edge list to versions you will work with
    nodes <- free.nodes
    kelp.edge.1 <- kelp.free.edge
    
    #calculating trophic level
    #convert edge list to adj matrix (using SPECIES ID instead of Node.ID)
    edge.only.sp <- kelp.edge.1[,c(4,9)] 
    edge.only.sp <- unique(edge.only.sp)
    g.sp = graph.data.frame(edge.only.sp)
    kelp.mat.sp <- get.adjacency(g.sp)  #get adjacency matrix from graph object
    kelp.mat.sp <- as.matrix(kelp.mat.sp)  #convert to numeric matrix
    
    #calculate trophic level and add to node list
    kelp.matT.sp <- t(kelp.mat.sp)
    tl.sp <- TrophInd(kelp.matT.sp)  # calculate trophic level, needs numeric matrix
    tl.sp$sp.ID <- rownames(tl.sp)
    tl.sp$sp.ID <- as.numeric(tl.sp$sp.ID)
    
    #join TL to node list, this also drops any of nodes that aren't in the edge list, as only those have trophic levels
    nodes.tl.sp <- right_join(nodes, tl.sp, by = "sp.ID")  
    ##create sp level node list
    collapse_unique <- function(x) {
      paste(unique(x), collapse = ",")
    }
    free.species.nodes <- ddply(nodes.tl.sp,.(sp.ID),colwise(collapse_unique))  
    
    #reorder edgelist so resourcespID is first column, consumer SpID is second, for reading into igraph
    #Note: this edge list is not yet aggregated to species. creation of the igraph object will remove duplicate links. 
    kelp.edge.sp.free <- kelp.edge.1[,c(9,4,1:3,5:8,10:17)]
    
    ##create igraph object, for graphing and analysis
    kelp.g.sp.free<- graph.data.frame(kelp.edge.sp.free,
                                                   directed = T,
                                                   vertices = free.species.nodes)
    
    

    
## Calculate metrics for web of interest
    ## Example: sp web w parasites and concomitant links
        #Load igraph object of web of interest
        kelp.g2 <- kelp.g.sp.parasites.w.conc
        
        #Node degree, generality, and vulnerability
        deg <- degree(kelp.g2, mode="all")
        vk <- degree(kelp.g2, mode="out")
        gk <- degree(kelp.g2, mode="in")
        
        deg.dat <- as.data.frame(deg)
        vuln.dat <- as.data.frame(vk)
        gen.dat <- as.data.frame(gk)
        
        degs.all <- cbind(deg.dat, vuln.dat)
        degs.all <- cbind(degs.all, gen.dat)
        degs.all$Sp.ID <- rownames(degs.all)
        
        #bind TL to degree list for plotting if desired. 
        degs.all$TL <- nodes.tl.sp$TL[match(degs.all$Sp.ID, nodes.tl.sp$sp.ID)]
        
        #convert graph object to adjacency matrix, and use GenInd() to calculate multiple summary metrics
        #GenInd() requires input of adjacency matrix, calculates multiple metrics (incl connectance)
        mat <- get.adjacency(kelp.g2,sparse = F)  
        M <- GenInd(mat)  
        
        #additional metrics of interest:
        D <- diameter(kelp.g2, directed=T, weights=NA)  #longest chain length
        Mean_Degree = mean(deg)
        SD_Degree = sd(deg)
        Mean_Generality = mean(gk)
        SD_Generality = sd(gk)
        Mean_Vulnerability = mean(vk)
        SD_Vulnerability = sd(vk)
        
        ## tabulating motif proportions
        #names following McLaughlin 2018. order of motifs in output in R help for triad_census
        motif.names <- c("M0", "M1", "M2", "S4", "S5", "S1", "D3", "D4", "S2", "S3", "D8", "D1", "D2", "D5", "D7", "D6")
        
        motifs.free<- as.data.frame(triad_census(kelp.g2))
        motifs.free$Motif <- motif.names
        
        
        
        
        
