AiccWeights <- function(trt,Variable){
  library(geiger)
  library(caper)
  #get tree with matching tips and star tree with no information
  tree <- trt$phy
  notree <- geiger::rescale(tree, model = "lambda", 0)
  data_to_fit <- trt$data[,Variable]
  names(data_to_fit) <- rownames(trt$data)
  #Compute candidate and null models
  nosigModel <- fitContinuous(notree, data_to_fit)
  brownianModel <- fitContinuous(tree, data_to_fit)
  lambdaModel <- fitContinuous(tree, data_to_fit, model = "lambda")
  OUModel <- fitContinuous(tree, data_to_fit, model = "OU")
  EBModel <- fitContinuous(tree, data_to_fit, model = "EB")
  #extract AICC
  zeroAICC <- nosigModel$opt$aicc
  bmAICC <- brownianModel$opt$aicc
  ouAICC <- OUModel$opt$aicc
  ebAICC <- EBModel$opt$aicc
  #compute AIC weights
  aicc <- c(zeroAICC, bmAICC, ouAICC, ebAICC)
  aiccD <- aicc - min(aicc)
  aw <- exp(-0.5 * aiccD)
  aiccW <- aw/sum(aw)
  names(aiccW) <- c("zeroAICC", "bmAICC", "ouAICC", "ebAICC")
  return(list(weights = aiccW, alpha = OUModel$opt$alpha))
}


AncestralReconstruction <- function(dataset,phylogeny){
  library(phytools)
  ancestra_fit <- anc.ML(trt$phy,tips_traits,model="OU")
  tips_traits<-as.matrix(trt$data)[,1]
  all_nodes_traits <- append(ancestra_fit$ace, tips_traits)
  plotBranchbyTrait(trt$phy, all_nodes_traits, mode="nodes")
}