noiser <-  function(unw_adj,mean_noise,var_noise,min_noise=0.0001){
  adj_dim <- dim(unw_adj)
  rnorm(adj_dim[1] * adj_dim[2],mean_noise,var_noise) %>% 
    pmax(min_noise) %>%
    pmin(1) %>% 
    matrix(nrow=adj_dim[1],ncol=adj_dim[2],byrow=TRUE) %>%
  "*"(unw_adj) %>%
    return
}

comparison_flow <- function(unw_adj,Mode,Rank,mean_noise=0.5,var_noise=0.5){
  species_names <- unw_adj %>% row.names
  unw_adj %>%
    noiser(mean_noise,var_noise) %>%
    compute_strain(Rank,mode=Mode) %>%
    setNames(species_names) %>%
    sort(decreasing = TRUE) %>%
  return
}

comparison_flow_meand <- function(unw_adj,Mode,Rank,mean_noise=0.5,var_noise=0.5){
  if(Mode == "inward") asfwe <- ASFWE_in
  if(Mode == "outward") asfwe <- ASFWE_out
  if(Mode == "total") asfwe <- ASFWE
  species_names <- unw_adj %>% row.names
  unw_adj %>%
    noiser(mean_noise,var_noise) %>%
    asfwe(Rank) %>% dist %>%
    as.matrix %>% rowMeans %>%
    setNames(row.names(unw_adj)) %>%
    return
}

comparison_flow_fd <- function(unw_adj,Mode,Rank,mean_noise=0.5,var_noise=0.5){
  species_names <- unw_adj %>% row.names
  unw_adj %>%
    noiser(mean_noise,var_noise) %>%
    compute_strain(Rank,mode=Mode) %>%
    setNames(species_names) %>%
    sort(decreasing = TRUE) %>%
    return
}