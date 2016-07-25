plot_adjacency <- function(adjacency_df,Web_name,adjacency_kind,Title){
  adjacency_df %>%
  filter(name_web %in% Web_name) %>%
  select_(adjacency_kind) %>%
  unnest() %>%
  as.matrix() %>%
  melt() %>%
  ggplot(aes(Var1,Var2, fill=value)) +
  geom_raster(show.legend = FALSE) +
  coord_fixed() +
  theme_minimal() +
  scale_fill_continuous(guide = FALSE) +
  theme(axis.line        = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank()) +
    ggtitle(Title)
}

plot_traits <- function(adjacency_df,Web_name,Traits_kind,Wardness){
  adjacency_df %>%
    filter(name_web %in% Web_name) %>%
    unnest_(Traits_kind) %>%
    filter(mode %in% Wardness) %>%
    ggplot(aes(trait.1,trait.2))  +
    stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) +
    geom_point(colour = "white")   +
    xlab("Coordinate 1") +
    ylab("Coordinate 2") +
    ggtitle(paste0(Web_name," trait space")) +
    theme_minimal() 
}


plot_uniqs <- function(adjacency_df,Web_name,Traits_kind,Wardness){
  expose_uniqueness(adjacency_df,Web_name,Traits_kind[1],"inward") %>%
    left_join(expose_uniqueness(adjacency_df,Web_name,Traits_kind[2],"inward"),
              by = c("node_name")) %>%
    ggplot(aes(MeanDist.x,MeanDist.y)) +
    geom_point() + geom_smooth()   +
    xlab(Traits_kind[1]) +
    ylab(Traits_kind[2]) +
    ggtitle(paste0(Web_name,", inward")) +
    theme_minimal()
}

plot_deg_uniq <- function(adjacency_df,Web_name,Traits_kind,Wardness){
  expose_uniqueness(adjacency_df,Web_name,Traits_kind[1],Wardness) %>%
    left_join(
      expose_degree(adjacency_df,Web_name,Traits_kind[2]),
      by = "node_name"
    ) %>%
    ggplot(aes_string("MeanDist",Traits_kind[2])) +
    geom_point() +
    xlab("Ecological uniqueness") +
    ylab("Degree Centrality") +
    ggtitle(paste0(Web_name," food web, outward")) +
    theme_minimal()
}

plot_summary_stat <- function(adjacency_df,spaces,variable){
  variables <- stringr::str_c(spaces,variable,sep = ".")
  
  adjacency_df %>%
    select_(spaces[1],spaces[2]) %>%
    unnest(.sep = ".") %>%
    ggplot(aes_string(variables[1],variables[2])) +
    geom_hline(yintercept = 0.5) +
    geom_vline(xintercept = 0.5) +
    geom_point() +
    xlab(paste0(variable," outward")) + xlim(0,1) +
    ylab(paste0(variable," inward")) + ylim(0,1) +
    ggtitle("Degree Centrality ~ Uniqueness") +
    theme_minimal()
}
