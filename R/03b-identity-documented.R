# ==========================================================================
# SVMC: pairwise genome identity (MASH)
# MASH command generation and pairwise distance/identity extraction.
# ==========================================================================

# Note: the simpler duplicate mash_commands() was removed; the general
# version (files/threads/sketch_size/msh_name args) is kept.

get_dists <- 
  function(distance_matrix){
    diag(distance_matrix) <- NA
    sequence_averages <- rowMeans(distance_matrix, na.rm = TRUE)
    all_distances <- as.vector(distance_matrix)
    all_distances <- 
      lapply(all_distances, function(x){x[!is.na(x)]})
    #!all_distances[!is.na(all_distances)]  # Remove NA values
    overall_mean_distance <-
      mean(unlist(all_distances))
    #lapply(all_distances, mean)
    #mean(all_distances)
    difference_from_mean <- sequence_averages - overall_mean_distance
    
    distance_table <- 
      tibble(
        Sequence = names(sequence_averages), 
        Sequence_mean_distance = sequence_averages, 
        Mean_distance = overall_mean_distance,
        Difference_from_Mean = difference_from_mean
      )
    return(distance_table)
  }

mash_commands <- function(species_name,
                          identity_dir,
                          sync_dir,
                          files = NULL,          # NULL => all *.txt in sync_dir
                          threads = 8,
                          sketch_size = 10000,
                          msh_name = "mash.msh") {
  
  # Build input list for `mash sketch`
  inputs <- if (is.null(files)) {
    sprintf("%s/*.txt", sync_dir)                         # glob all
  } else {
    # accept basenames or relative paths; always join to sync_dir and quote
    paste(shQuote(file.path(sync_dir, files)), collapse = " ")
  }
  
  msh_path <- file.path(identity_dir, msh_name)
  out_path <- file.path(identity_dir, sprintf("%s_mash.txt", species_name))
  
  sketch_cmd <- sprintf("mash sketch -p %d -o %s -s %d %s",
                        threads, msh_path, sketch_size, inputs)
  
  # self-distance of the (sub)set; if you want subset vs full, make a second sketch
  dist_cmd   <- sprintf("mash dist -p %d %s %s > %s",
                        threads, msh_path, msh_path, out_path)
  
  c(sketch_cmd, dist_cmd)
}

get_pairwise_identities <- 
  function(
    sync_dir, id_dir,
    species_name, clust_method = "complete"
  ){
    
    mash_file <- 
      sprintf(
        "%s/%s_mash.txt", 
        id_dir, species_name
      )
    # run mash if it hasn't been run
    if (!
        file.exists(
          mash_file
        )
    ){
      mash_bash <- 
        mash_commands(species_name, id_dir, sync_dir)
      system(mash_bash[1], ignore.stderr = T)
      system(mash_bash[2])
    }
    
    
    id_table <-
      fread(
        mash_file,
        select = c(1:3)
      ) %>%
      mutate(V1 = sub(".txt", "", basename(V1))) %>%
      mutate(V2 = sub(".txt", "", basename(V2))) %>%
      `colnames<-`(c("ref", "qry", "identity")) %>%
      mutate(identity = as.numeric(identity)) %>%
      as_tibble %>%
      mutate(
        identity = (1-identity)*100
      )

    id_table_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
      as.data.frame() %>%
      column_to_rownames("ref")
    
    ## for clustering purposes
    id_dist <- 
      100-as.matrix(id_table_wide) 
    id_dist[upper.tri(id_dist, diag = T)]<-NA
    
    id_dist <- 
      as.dist(id_dist)
    id_clust <- 
      hclust(id_dist, method = clust_method)
    clust_order <-
      id_clust$labels[id_clust$order]
    
    id_table_long <- 
      id_table_wide %>%
      rownames_to_column("ref") %>%
      gather(key = qry,value = "identity", na.rm = FALSE, -c(ref)) %>%
      group_by(ref) %>% 
      mutate(identity_mean = mean(identity)) %>%
      arrange(desc(identity_mean)) %>%
      ungroup
        
        #Recreate the wide matrices
    id_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
      as.data.frame() %>%
      column_to_rownames("ref") %>%
      as.matrix()
    
    id_wide_dist <- 
      100 - id_wide
    
    id_sld_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
      as.data.frame() %>%
      column_to_rownames("ref")
    
    id_sld_wide_dist <- 
      100 - id_sld_wide
    
    
    distance_matrix <- 
      id_wide_dist
    
    id_dists <- 
      get_dists(distance_matrix) %>%
      arrange(Sequence_mean_distance)
    
    id_sld_dists <-
      get_dists(id_sld_wide_dist) %>%
      arrange(Sequence_mean_distance)
    
    id_dist_lt <- 
      as.dist(distance_matrix)
    
    seq_order <- 
      colnames(distance_matrix)[
        hclust( id_dist_lt, method = "ward.D" )$order
      ]
    
    #
    mash_plot <- 
      ggplot(
        data = id_table_long, 
        aes(
          x= factor(ref, seq_order),
          y= factor(qry, seq_order),
          fill = identity#_sld#identity
        )
      ) + 
      geom_tile() + 
      theme_classic() +
      theme(
        #axis.text.x = element_text(angle=45, hjust = 1),
        panel.background=element_rect(fill="white"),
        plot.title = element_text(size=20,face="bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank()
      ) + 
      ggtitle(
        sprintf(
          "%s \n Sequence pairwise MASH similarity", 
          sub("_", " ", species_name)
        )
      ) + 
      theme(
        legend.title=element_text( size=18),
        axis.text=element_blank(),
        #element_text(size=14), 
        legend.text = element_text(size=18)
      ) +
      scale_fill_viridis(
        option = "B", limits = c(95,100)
      ) 
    
    mash_list <- 
      list(
        mash_tib = id_table,
        mash_plot = mash_plot,
        id_dists = id_dists
      )
    return(mash_list)
  }
