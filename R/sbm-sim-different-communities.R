# =============================================================================
# sbm-sim-different-communities.R  --  SBM simulation study (Supp. Tables S16/S17)
# -----------------------------------------------------------------------------
# Standalone supporting analysis: simulates stochastic block models across
# k = 2..6 communities and low/high between-community connection probabilities,
# then classifies them with the first-stage classifier.
#
# This script is NOT part of run_all.R. Run it on its own:
#     source(here::here("R", "00_setup.R"))
#     source(here::here("R", "sbm-sim-different-communities.R"))
# =============================================================================

if (!exists("dir_outputs")) source(here::here("R", "00_setup.R"))

## ---- Model route ------------------------------------------------------------
## Use the version-portable JSON booster + trained recipe, exactly as scripts 02,
## 05 and 06 do. This avoids the superseded xgb_final_trained_model.rds workflow,
## which fails to reload under current xgboost. Same booster, same trained recipe,
## so predictions are unchanged. predict_wf() returns a .pred_class column, which
## is what the rest of this script expects.
if (!exists("load_main_model")) source(here::here("R", "01_functions.R"))

main_model <- load_main_model()

predict_wf <- function(new_data) {
  tibble::tibble(
    .pred_class = factor(predict_main_model(main_model, new_data)$predicted_class,
                         levels = CLASS_ORDER)
  )
}


###--Normalized Laplacina function--##
normalized_laplacian=function(Graphs){
  laplacian_matrix(Graphs,normalized = T)
}



calcGraphFeatures <- function(Graphs=NULL) {
  
  features <- c(
    "order",                    # number of vertices
    "edges",                     # number of edges
    "connected",                # True / False
    "max_component",            # maximum component size (=order iff the graph is connected)
    "minDegree",                # minimum degree of any vertex
    "maxDegree",                # maximum degree of any vertex
    "mean_degree",                # average degree of any vertex
    "minCut",                   # minimum cut weight of the graph (might take a while to compute)
    "FiedlerValue",             # second-highest eigenvalue of the Laplacian matrix
    "Normalized_FiedlerValue",   # second-highest eigenvalue of the Normaized Laplacian matrix
    "closeness_centr",                # average inverse of distance between any pair of vertices
    "modularity",               # DEFINITION REQUIRED
    "diameter",                 # maximum distance between any two vertices (NAN if not connected)
    "betw_centr",              # max_{v} proportion of shortest paths going through vertex v
    "transitivity",             # aka Clustering Coefficient, is proportion of connected triples that form triangles: e.g., (a--b--c--a) when (a--b--c) is present.
    "threshold",                 # 1/max(eigen value of A)
    "spectral_radius"         # max (eigen value of A)
    
  )
  
  df <- as.data.frame(matrix(ncol=length(features),nrow=length(Graphs)))
  colnames(df)=features
  
  # Stuff that is simple to apply and needs no interim components:
  
  df$order = base::as.numeric(lapply(Graphs, gorder))
  df$edges = base::as.numeric(lapply(Graphs, gsize))
  df$connected = base::as.numeric(lapply(Graphs, is_connected))
  df$minCut = base::as.numeric(lapply(Graphs, min_cut))
  df$diameter = base::as.numeric(lapply(Graphs, diameter))
  df$transitivity = base::as.numeric(lapply(Graphs, transitivity))
  
  # stuff that needs interim things:
  degrees = lapply(Graphs,igraph::degree )
  df$minDegree = base::as.numeric(lapply(degrees, min))
  df$maxDegree = base::as.numeric(lapply(degrees, max))
  df$mean_degree = base::as.numeric(lapply(degrees, mean))
  
  # stuff that lapply doesn't like so has to be done in a loop:
  communities <-lapply(Graphs, cluster_walktrap) #lapply(Graphs, cluster_leading_eigen)
  Adj <- lapply(Graphs, as_adjacency_matrix)
  L <- lapply(Graphs, laplacian_matrix)
  Norm_Lap<-lapply(Graphs, normalized_laplacian)
  Fiedler.value=NULL
  norm.fiedler.value=NULL
  
  for (i in 1:length(Graphs)) {
    if (is.null(Graphs[[i]]$type)) { Graphs[[i]]$type = "untyped" }
    df$modularity[i] <- modularity(communities[[i]])
    df$spectral_radius[i] <- eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values[1]
    
    Fiedler.value[[i]]=eigen(L[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$FiedlerValue[i] <- Fiedler.value[[i]][length(Fiedler.value[[i]])-1]
    
    norm.fiedler.value[[i]]=eigen(Norm_Lap[[i]], symmetric=TRUE, only.values=TRUE)$values
    
    df$Normalized_FiedlerValue[i] <- norm.fiedler.value[[i]][length(norm.fiedler.value[[i]])-1]
    
    df$eigen_centr[i] <- centr_eigen(Graphs[[i]])$centralization
    df$deg_centr[i] <- centr_degree(Graphs[[i]])$centralization
    df$betw_centr[i] <- centr_betw(Graphs[[i]])$centralization
    
    df$max_component[i] <- max(components(Graphs[[i]])$csize)
    df$mean_eccentr[i]<-mean(eccentricity(Graphs[[i]]))
    df$radius[i]<-radius(Graphs[[i]])
    df$mean_path_length[i]<-mean_distance(Graphs[[i]])
    #df$trace[i]<-sum(diag(Adj[[i]]))
    df$graph_energy[i]<-sum(abs(eigen(Adj[[i]], symmetric=TRUE, only.values=TRUE)$values))
    df$min_triangle[i]= min(count_triangles(Graphs[[i]]))
    df$mean_triangle[i]= mean(count_triangles(Graphs[[i]]))
    df$sd_triangle[i]= sd(count_triangles(Graphs[[i]]))
    df$max_triangle[i]= max(count_triangles(Graphs[[i]]))
    df$num_triangle[i]= sum(count_triangles(Graphs[[i]]))
    df$deg_assort_coef[i]=assortativity_degree(Graphs[[i]])
    
    df$threshold[i] <- 1/(df$spectral_radius[i])
    
    if (df$connected[i]==TRUE) {
      df$closeness_centr[i] = mean(closeness(Graphs[[i]]))
    } else { # handle the case where G isn't connected
      df$closeness_centr[i] = -1
    }
  }
  return (df)
}


RunSimOnGraphFeatures<-function(Graphs, nreps=nreps,output_file=NULL, seed=-1) {
  set.seed(1)
  # ### Definition and initialization of parameters for graphfeatures
  graphProperties=list()
  # ### Definition and initialization of Graph Prefix
  graphid=list();graphreplicate=list(); graphname=list(); GraphPrefix=list(); analysis=list()
  for (g in 1:length(Graphs)){
    for (reps in 1:nreps) {
      ### Calculate the graph features for each simulated graph of all the synthetic networks
      print(paste("Calculating graph features on", Graphs[[g]]$name))
      #graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
      graphProperties[[reps]] <- calcGraphFeatures(Graphs[g])
    }
    graphname[[g]]=Graphs[[g]]$type
    graphid[[g]]=Graphs[[g]]$id
    graphreplicate[[g]]=c(1:nreps)
    GraphPrefix=cbind(graphname[[g]],graphid[[g]],graphreplicate[[g]])
    colnames(GraphPrefix)=c("GraphName","GraphID","GraphReplicate")
    analysis[[g]]=as.data.frame(cbind(GraphPrefix,graphProperties[[reps]]))
    row.names(analysis[[g]])=1:nreps
  }
  All_results=do.call(rbind,analysis)
  # write.csv(All_results, file=output_file)
  return( All_results)
}


#+++++++++++++++++++++++++++++++++++++++++++++++
#-------MAKING OF SPATIAL GRAPH---
#+++++++++++++++++++++++++++++++++++++++++++++++

# Helper functions for spatial graphs to convert between (row, column) pairs and the index in a list of Cells.
cellCoordsToIndex <- function(i, j, size) {
  return((i-1)*size+j)
}

indexToCellCoords <- function(idx, size) {
  j <- (idx-1) %% size + 1
  i <- (idx + size-1) %/% size
  return(c(i,j))
}





# Function to create SBM with specified parameters
MKSBMGraphs <- function(pm, n, k) {
  Graphs = list()
  i = 1
  print(paste("Creating SBM with", k, "communities and", n, "nodes"))
  # Calculate block sizes - equal distribution across communities
  block_sizes = rep(floor(n/k), k)
  # Adjust for any remainder
  block_sizes[1] = block_sizes[1] + (n - sum(block_sizes))
  
  Graphs[[i]] = sample_sbm(n, pref.matrix = pm, block.sizes = block_sizes)
  Graphs[[i]]$type = "sbm"
  Graphs[[i]]$id = paste0(k, "_communities")
  Graphs[[i]]$name = paste0("sbm_k", k)
  i <- i+1
  return(Graphs)
}


# Improved probability matrix generator for k communities
sbmfunc <- function(k) {
  # Ensure within-community probability is higher than between-community
  p_within <- runif(k, 0.1, 0.15)  # Higher within-community probability
  p_between <- runif(1, 0.01, 0.03) # Lower between-community probability
  
  mat = matrix(p_between, nrow = k, ncol = k)  # Fill matrix with between-community probability
  diag(mat) <- p_within  # Set diagonal to within-community probability
  
  return(mat)
}

# Function to simulate networks with different sizes and communities
sim.sbm.net <- function(n_nodes, k, n_instances=100) {
  print(paste("Simulating", n_instances, "networks with", k, "communities and", n_nodes, "nodes"))
  
  # Generate probability matrices
  pm_list = list()
  for (i in 1:n_instances) {
    pm_list[[i]] = sbmfunc(k)
  }
  
  # Generate networks
  net = lapply(pm_list, MKSBMGraphs, n=n_nodes, k=k)
  data = lapply(net, RunSimOnGraphFeatures, nreps=1)
  df = as.data.frame(do.call(rbind, data))
  
  # Add metadata
  df$n_nodes = n_nodes
  df$n_communities = k
  
  return(df)
}

# Parameters
node_sizes = 100#c(50, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000)
k_values = 4#2:6
n_instances = 4#100  # Number of instances per combination

# Create empty list to store results
all_results = list()

# Generate networks for all combinations
for(n in node_sizes) {
  for(k in k_values) {
    result = sim.sbm.net(n_nodes=n, k=k, n_instances=n_instances)
    all_results[[paste0("n", n, "_k", k)]] = result
  }
}

# Combine all results
final_df = do.call(rbind, all_results)

# Add timestamp and save results
timestamp = format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(final_df, file=paste0("sbm_graphs_simulation_results_", timestamp, ".rds"))





####################################################################################

# Function to create SBM with specified parameters
MKSBMGraphs <- function(pm, n, k) {
  print(paste("Creating SBM with", k, "communities and", n, "nodes"))
  
  # Calculate block sizes - equal distribution across communities
  block_sizes = rep(floor(n/k), k)
  # Adjust for any remainder
  block_sizes[1] = block_sizes[1] + (n - sum(block_sizes))
  
  # Create the graph
  g = sample_sbm(n, pref.matrix = pm, block.sizes = block_sizes)
  
  # Add community membership as vertex attribute
  V(g)$community = rep(1:k, times = block_sizes)
  
  # Add metadata
  g$type = "sbm"
  g$id = paste0(k, "_communities")
  g$name = paste0("sbm_k", k)
  
  return(g)
}

# Function to plot SBM with communities
plot_sbm <- function(g, title = "") {
  # Get community colors
  n_comm = length(unique(V(g)$community))
  comm_colors = rainbow(n_comm)
  
  # Create layout
  layout = layout_with_fr(g)
  
  # Create plot
  png(paste0("sbm_plot_", title, ".png"), width = 800, height = 800)
  plot(g,
       vertex.color = comm_colors[V(g)$community],
       vertex.size = 5,
       vertex.label = NA,
       layout = layout,
       main = title)
  dev.off()
  
  return(g)
}

# Improved probability matrix generator for k communities
sbmfunc <- function(k) {
  # Ensure within-community probability is higher than between-community
  p_within <- runif(k, 0.5, 0.9)  # Higher within-community probability
  p_between <- runif(1, 0.01, 0.04) # Lower between-community probability
  
  mat = matrix(p_between, nrow = k, ncol = k)  # Fill matrix with between-community probability
  diag(mat) <- p_within  # Set diagonal to within-community probability
  
  return(mat)
}

# Function to simulate networks with different sizes and communities
sim.sbm.net <- function(n_nodes, k, n_instances=100) {
  print(paste("Simulating", n_instances, "networks with", k, "communities and", n_nodes, "nodes"))
  
  # Lists to store graphs and features
  graphs_list = list()
  features_list = list()
  
  for (i in 1:n_instances) {
    # Generate probability matrix
    pm = sbmfunc(k)
    
    # Generate graph
    g = MKSBMGraphs(pm, n_nodes, k)
    
    # Save graph
    graphs_list[[i]] = g
    
    # Plot graph
    title = paste0("n", n_nodes, "_k", k, "_instance", i)
    plot_sbm(g, title)
    
    # Calculate features
    features = RunSimOnGraphFeatures(list(g), nreps=1)
    features_list[[i]] = features
  }
  
  # Combine features
  df = as.data.frame(do.call(rbind, features_list))
  
  # Add metadata
  df$n_nodes = n_nodes
  df$n_communities = k
  
  # Save graphs
  timestamp = format(Sys.time(), "%Y%m%d_%H%M")
  saveRDS(graphs_list, file=paste0("sbm_graphs_n", n_nodes, "_k", k, "_", timestamp, ".rds"))
  
  return(list(features = df, graphs = graphs_list))
}

# Parameters
node_sizes = 100
k_values = 4
n_instances = 3

# Create empty lists to store results
all_features = list()
all_graphs = list()

# Generate networks for all combinations
for(n in node_sizes) {
  for(k in k_values) {
    result = sim.sbm.net(n_nodes=n, k=k, n_instances=n_instances)
    all_features[[paste0("n", n, "_k", k)]] = result$features
    all_graphs[[paste0("n", n, "_k", k)]] = result$graphs
  }
}

# Combine all feature results
final_df = do.call(rbind, all_features)

# Save results
#timestamp = format(Sys.time(), "%Y%m%d_%H%M")
saveRDS(final_df, file=paste0("sbm_features_", ".rds"))
saveRDS(all_graphs, file=paste0("sbm_all_graphs_", ".rds"))


# Create summary plots
plot_summary <- function(graphs_list, title) {
  plots = list()
  for(i in 1:length(graphs_list)) {
    g = graphs_list[[i]]
    plots[[i]] = plot_sbm(g, paste0(title, "_", i))
  }
  return(plots)
}

# Generate summary plots for each combination
for(n in node_sizes) {
  for(k in k_values) {
    key = paste0("n", n, "_k", k)
    plot_summary(all_graphs[[key]], key)
  }
}

x=plot_summary(all_graphs[[key]], key)

plot(x[[2]])
#x=readRDS("sbm_graphs_simulation_results.rds")
x

# str(x)
# x$n_communities

x1=final_df%>%
  dplyr::filter(connected==1)
x1

nrow(x1)
colnames(x1)

#emp.data = emp.data[sample(1:nrow(emp.data)), ]
sim.data=x1%>%
  dplyr::rename("graph_name"="GraphName")%>%
  #dplyr::rename("NumOfNodes"="order")%>%
  dplyr::select(c(graph_name,order,edges,
                  mean_eccentr,mean_path_length,graph_energy,
                  modularity,diameter,betw_centr,transitivity,
                  spectral_radius,eigen_centr,deg_centr,
                  mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                  closeness_centr,deg_assort_coef))%>%
  mutate_if(is.character,factor)%>%
  clean_names()%>%
  as_tibble()

colSums(is.na(sim.data))


nrow(sim.data)

# emp.data=emp.data%>%
#   filter(order>10)
# 
# nrow(emp.data)
# colSums(is.na(emp.data))

# emp.data.rec=recipe(graph_name~., data = emp.data)%>%
#   step_impute_median(all_numeric_predictors())
# emp.data.prep=emp.data.rec%>%prep()
# emp.data.juice=emp.data.prep%>%juice()
# 
# saveRDS(emp.data.juice,"final.emp.data.rds")
# final.emp.data=readRDS("final.emp.data.rds")
# colSums(is.na(final.emp.data))

#extract_final.fitted.model=saveRDS(extract_final.fitted.model,"extract_final.fitted.model.rds")
#final.emp.data=readRDS("final.emp.data.rds")

####----Final preprocessed data frame----for empirical data----##
best.model.predictions <- predict_wf(sim.data) %>%
  bind_cols(sim.data)


best.model.predictions






###----data frame for predicted empirical network
# df.emp=data.frame(best.model.predictions$graph_name,best.model.predictions$.pred_class)
# colnames(df.emp)=c("target","Predicted_classes")  


x2=read.csv(file.path(dir_data, "Reduced-Animal-Social-network.csv"))
mean(x2$modularity)
range(x2$modularity)
summary(x2$modularity)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# New function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(igraph)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(janitor)
library(tidymodels)

generate_sbm_networks <- function(node_sizes, k_values, n_instances, 
                                  p_between_range=c(0.1, 0.4)) {
  # Input validation
  if(!is.numeric(node_sizes) || !is.numeric(k_values) || !is.numeric(n_instances)) {
    stop("All inputs must be numeric")
  }
  
  # Improved probability matrix generator for k communities
  sbmfunc <- function(k, p_between_range) {
    # Generate p_within for this specific matrix
    p_within <- runif(k, 0.5, 0.8)  # Higher within-community probability
    
    # Create symmetric matrix
    mat <- matrix(0, nrow = k, ncol = k)
    
    # Fill upper triangle with between-community probabilities
    for(i in 1:(k-1)) {
      for(j in (i+1):k) {
        mat[i,j] <- runif(1, p_between_range[1], p_between_range[2])
        mat[j,i] <- mat[i,j]  # Make symmetric
      }
    }
    
    # Set diagonal to within-community probabilities
    diag(mat) <- p_within
    
    return(list(matrix = mat, p_within = p_within))
  }
  
  # Function to create SBM with specified parameters
  MKSBMGraphs <- function(pm, n, k) {
    print(paste("Creating SBM with", k, "communities and", n, "nodes"))
    
    # Calculate block sizes - equal distribution across communities
    block_sizes = rep(floor(n/k), k)
    # Adjust for any remainder
    block_sizes[1] = block_sizes[1] + (n - sum(block_sizes))
    
    # Create the graph
    g = sample_sbm(n, pref.matrix = pm, block.sizes = block_sizes)
    
    # Add community membership as vertex attribute
    V(g)$community = rep(1:k, times = block_sizes)
    
    # Add metadata
    g$type = "sbm"
    g$id = paste0(k, "_communities")
    g$name = paste0("sbm_k", k)
    
    return(g)
  }
  
  # Function to plot SBM
  plot_sbm <- function(g, title = "") {
    # Get community colors
    n_comm = length(unique(V(g)$community))
    comm_colors = rainbow(n_comm)
    
    # Create layout
    layout = layout_with_fr(g)
    
    # Create plot
    p <- plot(g,
              vertex.color = comm_colors[V(g)$community],
              vertex.size = 5,
              vertex.label = NA,
              layout = layout,
              main = title)
    
    return(p)
  }
  
  # Create results list
  results <- list(
    features = list(),
    graphs = list(),
    probability_matrices = list(),
    p_within = list(),
    metadata = list(
      node_sizes = node_sizes,
      k_values = k_values,
      n_instances = n_instances,
      p_between_range = p_between_range,
      timestamp = format(Sys.time(), "%Y%m%d_%H%M")
    )
  )
  
  # Initialize empty dataframe for all predictions
  all_predictions <- data.frame()
  
  # Load the trained model
  
  # Generate networks for all combinations
  for(n in node_sizes) {
    for(k in k_values) {
      graphs_list <- list()
      features_list <- list()
      pm_list <- list()
      p_within_list <- list()
      
      for(i in 1:n_instances) {
        # Generate probability matrix
        pm_result <- sbmfunc(k, p_between_range)
        pm <- pm_result$matrix
        p_within <- pm_result$p_within
        
        pm_list[[i]] <- pm
        p_within_list[[i]] <- p_within
        
        # Verify matrix is symmetric
        if(!isSymmetric(pm)) {
          stop("Generated probability matrix is not symmetric")
        }
        
        # Generate graph
        g <- MKSBMGraphs(pm, n, k)
        graphs_list[[i]] <- g
        
        # Calculate features
        features <- RunSimOnGraphFeatures(list(g), nreps=1)
        features_list[[i]] <- features
        
        # Prepare data for prediction
        sim.data <- as.data.frame(features) %>%
          dplyr::filter(connected == 1) %>%
          dplyr::rename("graph_name" = "GraphName") %>%
          dplyr::select(c(graph_name, order, edges,
                          mean_eccentr, mean_path_length, graph_energy,
                          modularity, diameter, betw_centr, transitivity,
                          spectral_radius, eigen_centr, deg_centr,
                          mean_degree, minCut, FiedlerValue, Normalized_FiedlerValue,
                          closeness_centr, deg_assort_coef)) %>%
          mutate_if(is.character, factor) %>%
          clean_names() %>%
          as_tibble()
        
        # Make predictions
        if(nrow(sim.data) > 0) {
          predictions <- predict_wf(sim.data) %>%
            bind_cols(sim.data) %>%
            mutate(
              n_nodes = n,
              k_communities = k,
              instance = i,
              p_within = list(p_within)  # Store as list column
            )
          
          all_predictions <- bind_rows(all_predictions, predictions)
        }
      }
      
      key <- paste0("n", n, "_k", k)
      results$graphs[[key]] <- graphs_list
      results$features[[key]] <- do.call(rbind, features_list)
      results$probability_matrices[[key]] <- pm_list
      results$p_within[[key]] <- p_within_list
    }
  }
  
  # Store all predictions in results
  results$predictions <- all_predictions
  
  # Add feature access function
  results$get_features <- function(n, k) {
    key <- paste0("n", n, "_k", k)
    if(!key %in% names(results$features)) {
      stop("No features found for specified n and k values")
    }
    return(results$features[[key]])
  }
  
  # Add graph access function
  results$get_graphs <- function(n, k) {
    key <- paste0("n", n, "_k", k)
    if(!key %in% names(results$graphs)) {
      stop("No graphs found for specified n and k values")
    }
    return(results$graphs[[key]])
  }
  
  # Add probability matrix access function
  results$get_probability_matrices <- function(n, k) {
    key <- paste0("n", n, "_k", k)
    if(!key %in% names(results$probability_matrices)) {
      stop("No probability matrices found for specified n and k values")
    }
    return(results$probability_matrices[[key]])
  }
  
  # Add p_within access function
  results$get_p_within <- function(n, k) {
    key <- paste0("n", n, "_k", k)
    if(!key %in% names(results$p_within)) {
      stop("No p_within values found for specified n and k values")
    }
    return(results$p_within[[key]])
  }
  
  # Modify get_predictions to return filtered predictions
  results$get_predictions <- function(n, k) {
    results$predictions %>%
      filter(n_nodes == n, k_communities == k)
  }
  
  # Add summary function
  results$get_summary <- function(n, k) {
    predictions <- results$get_predictions(n, k)
    
    # Create summary
    summary <- list(
      total_networks = nrow(predictions),
      predicted_classes = table(predictions$.pred_class),
      features_summary = summary(predictions %>% select(-starts_with(".pred")))
    )
    
    return(summary)
  }
  
  # Add plotting function with predictions
  results$plot_graphs_with_predictions <- function(n, k) {
    graphs <- results$get_graphs(n, k)
    predictions <- results$get_predictions(n, k)
    
    plots <- list()
    for(i in seq_along(graphs)) {
      g <- graphs[[i]]
      pred <- predictions$.pred_class[predictions$instance == i][1]
      p_within <- unlist(predictions$p_within[predictions$instance == i][1])
      
      # title <- paste0("n=", n, ", k=", k, 
      #                 "\nPredicted: ", pred,
      #                 "\np_within: ", paste(round(p_within, 3), collapse=", "))
      plots[[i]] <- plot_sbm(g)#, title)
    }
    
    return(plots)
  }
  
  return(results)
}

# Example usage:
node_sizes <-100 #seq(100, 1000, by = 50)
k_values <- c(2, 3, 4, 5, 6)
n_instances <- 1

#---Generate networks for sbms with higher between community connection
results <- generate_sbm_networks(
  node_sizes = node_sizes,
  k_values = k_values,
  n_instances = n_instances,
  p_between_range = c(0.1, 0.4)
)

# Get all predictions
all_preds <- results$predictions
View(all_preds)

all_preds$p_within <- sapply(all_preds$p_within, function(x) paste(x, collapse = ","))
write.csv(all_preds, file.path(dir_outputs, "all_sbm_preds_higher_between_connectn_probs.csv"), row.names = FALSE)
mean(all_preds$modularity)#0.221



#---Generate networks for sbms with lower between community connection
results <- generate_sbm_networks(
  node_sizes = node_sizes,
  k_values = k_values,
  n_instances = n_instances,
  p_between_range = c(0.01, 0.04)
)

# Get all predictions
all_preds <- results$predictions
View(all_preds)

all_preds$p_within <- sapply(all_preds$p_within, function(x) paste(x, collapse = ","))
write.csv(all_preds, file.path(dir_outputs, "all_sbm_preds_lower_between_connectn_probs.csv"), row.names = FALSE)
mean(all_preds$modularity) #0.608

# Get predictions for specific parameters
preds_100_2 <- results$get_predictions(100, 2)

# Get summary for specific parameters
summary_100_2 <- results$get_summary(100, 2)

# Plot graphs with predictions
plots_100_2 <- results$plot_graphs_with_predictions(100, 6)

# Get p_within values
p_within_100_2 <- results$get_p_within(100, 2)

# Save all results
saveRDS(results, file = paste0("sbm_results_with_predictions_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"))





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Real world networks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##---Butterfly similarity network
edges_df1 <- read.table("butterfly-similarity.edges", header=TRUE)
# Increment node IDs by 1 in both columns
edges_df1$X0 <- edges_df1$X0 + 1
edges_df1$X4  <-edges_df1$X4  + 1

# Create igraph object using only the first two columns
g1 <- graph_from_edgelist(as.matrix(edges_df1[,1:2]), directed=FALSE)

##---Karate club network
edges_df2 <- read.table("soc-karate.edges", header=TRUE)
g2 <- graph_from_edgelist(as.matrix(edges_df2[,1:2]), directed=FALSE)

##---American football games network
edges_df3 <- read.table("misc-football.edges", header=TRUE)
g3 <- graph_from_edgelist(as.matrix(edges_df3[,1:2]), directed=FALSE)

##---Devil mating
edges_df4 <- read.table("Devil_mating.edges", header=TRUE)
# Increment node IDs by 1 in both columns
edges_df4$X0 <- edges_df4$X0 + 1
edges_df4$X1  <-edges_df4$X1  + 1

g4 <- graph_from_edgelist(as.matrix(edges_df4[,1:2]), directed=FALSE)


##---Devil non mating
edges_df5 <- read.table("Devil_nonmating.edges", header=TRUE)
# Increment node IDs by 1 in both columns
edges_df5$X0 <- edges_df5$X0 + 1
edges_df5$X1  <-edges_df5$X1  + 1

g5 <- graph_from_edgelist(as.matrix(edges_df5[,1:2]), directed=FALSE)

g=list(g1,g2,g3,g4,g5)
# 
# # Keep only the first two columns
# edges_df <- edges_df[, 1:2]
# 
# # Write the cleaned data
# write.table(edges_df, "misc-football.edges", row.names=FALSE, col.names=FALSE)




# Optional: Plot the graph
plot(g, vertex.size=3, vertex.label.cex=0.5)

# Calculate features

features = calcGraphFeatures(g)
                            
features$GraphNames <- c("Butterfly similarity network", 
                         "Karate club network", 
                         "American football games network", 
                         "Devil mating", 
                         "Devil non mating")

features <- features[, c("GraphNames", names(features)[names(features) != "GraphNames"])]



real.data=features%>%
  dplyr::rename("graph_name"="GraphNames")%>%
  #dplyr::rename("NumOfNodes"="order")%>%
  dplyr::select(c(order,edges,
                  mean_eccentr,mean_path_length,graph_energy,
                  modularity,diameter,betw_centr,transitivity,
                  spectral_radius,eigen_centr,deg_centr,
                  mean_degree,minCut,FiedlerValue,Normalized_FiedlerValue,
                  closeness_centr,deg_assort_coef))%>%
  mutate_if(is.character,factor)%>%
  clean_names()%>%
  as_tibble()

colSums(is.na(real.data))


nrow(real.data)



####----Final preprocessed data frame----for empirical data----##
best.model.predictions <- predict_wf(real.data) %>%
  bind_cols(real.data)


best.model.predictions


best.model.predictions$Graph_Names <- c("Butterfly similarity network", 
                                       "Karate club network", 
                                       "American football games network", 
                                       "Devil mating", 
                                       "Devil non mating")

best.model.predictions <- best.model.predictions[, c("Graph_Names", names(best.model.predictions)[names(best.model.predictions) != "Graph_Names"])]

View(best.model.predictions)

write.csv(best.model.predictions, file.path(dir_outputs, "real-world-datasets.csv")) 