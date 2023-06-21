## Collection of plots for paper

initialize <- TRUE
fig1 <- FALSE
fig1_init <- FALSE
fig2 <- TRUE
fig2_init <- TRUE
fig3 <- FALSE
fig_format <- "pdf"
limma <- FALSE
# pwr_calc <- tibble::tribble(~test, ~ES, ~d, ~sig.level, ~power, ~type, ~alternative)


# Initialize  ----
if(initialize){
  library("pairedGSEA")
  library("magrittr")
  library("dplyr")
  library("tidyr")
  library("purrr")
  library("ggplot2")
  library("patchwork")
  library("Cairo")
  theme_set(theme_classic(base_size = 19))
  
  # Both  DGE          DGS        Neither
  color_scheme <- c("gray", "lightblue", "darkred", "black")
  fill_scale <- function(reorder = NULL) {
    if(!is.null(reorder)) color_scheme <- color_scheme[reorder]
    scale_fill_manual(values = color_scheme, 
                      na.value = "white")}
  color_scale <- function(reorder = NULL) {
    if(!is.null(reorder)) color_scheme <- color_scheme[reorder]
    scale_color_manual(values = color_scheme, 
                       na.value = "white")}
  
  
  full_figure_theme <- function(){
    theme(plot.tag = element_text(face = "plain", size = 10))
  }
  
  utils::globalVariables(add = TRUE, c("color_scheme"))
  
  
  # Save metadata overview
  
  source("R/run_experiment.R")
  experiments <- combine_experiments(limma = limma)
  limma_suffix <- ifelse(limma, "_limma", "")
  experiments %>% select(-c("final_description", "group_nr", "id", "filename")) %>% 
    writexl::write_xlsx(paste0("figs/Experiments", limma_suffix, ".xlsx"))
  
  
}

# Run analysis ----

ora_all <- readRDS(paste0("results/ora", limma_suffix, "_all.RDS"))
gene_sets <- readRDS("results/gene_sets.RDS")
colnames(ora_all) <- colnames(ora_all) %>%
  stringr::str_replace("deseq", "expression") %>%
  stringr::str_replace("dexseq", "splicing")

quantile_overlap <- function(ora){
  
  print(ora$experiment[1])
  
  # Get gene sets for each analysis
  gs_expression <- ora %>% filter(padj_expression < 0.05) %>% pull(pathway)
  gs_splicing <- ora %>% filter(padj_splicing < 0.05) %>% pull(pathway)
  if(length(gs_expression) == 0 | length(gs_splicing) == 0) return(tibble("splicing" = NA, "expression" = NA))
  
  # Get relevant gene sets
  gs <- gene_sets[union(gs_expression, gs_splicing)]
  
  # Make a logical matrix of genes' presence in each gene set
  genes_tidy <- gs %>% 
    tibble::enframe(name = "gene_set", value = "gene") %>% 
    unnest(cols = c(gene))
  all_genes <- unique(genes_tidy$gene)
  
  ## Use pivot_wider to create the matrix
  gene_matrix <- genes_tidy %>%
    mutate(value = TRUE) %>%
    pivot_wider(names_from = gene_set, values_from = value, values_fill = FALSE)
  
  
  
  sample_expression_size <- min(length(gs_expression), length(gs_splicing))
  # Overlap coefficient for splicing vs expression and subset of expression vs expression
  oc_splicing <- proxy::simil(select(gene_matrix, all_of(gs_splicing)), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
  oc_expression <- proxy::simil(select(gene_matrix, all_of(sample(gs_expression, sample_expression_size))), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
  
  # Set diagonal to NA
  set_diag <- function(oc){
    oc_mat <- as.matrix(oc)
    name_combinations <- expand.grid(rownames(oc_mat), colnames(oc_mat))
    
    # Find matching indices based on name combinations
    matching_indices <- which(as.character(name_combinations$Var1) == as.character(name_combinations$Var2), arr.ind = TRUE)
    
    # Set the corresponding matrix elements to NA
    oc_mat[matching_indices] <- NA
    
    return(oc_mat)
  }
  oc_splicing <- set_diag(oc_splicing)
  oc_expression <- set_diag(oc_expression)
  
  # Compute 90th quantile
  quantiles_splicing <- apply(oc_splicing, 1, quantile, probs = 0.9, na.rm = TRUE)
  quantiles_expression <- apply(oc_expression, 1, quantile, probs = 0.9, na.rm = TRUE)
  
  return(tibble("splicing" = quantiles_splicing, "expression" = quantiles_expression))
}

system.time(
  mean_overlap <- map(unique(ora_all$experiment)[171:180], function(Experiment){
    ora_all %>% 
      filter(experiment == Experiment) %>% 
      group_by(experiment) %>% 
      summarise(results = quantile_overlap(.)) %>% 
      unnest_wider(results)
  }) %>% 
    bind_rows()
)

saveRDS(mean_overlap, file = "results/mean_overlap_171-180.RDS")

# "60_GSE47718_Effect of smoking on airway basal cells"
# "86_GSE153023_Folate high vs low"


# Plots ----

mean_overlap <- map(list.files("/home/projects/shd/pairedGSEA/results/", pattern = "mean_overlap*"), function(file){
  overlap <- readRDS(paste0("/home/projects/shd/pairedGSEA/results/", file))
}) %>% 
  bind_rows()

(plt_mean_overlap <- mean_overlap %>% 
  # filter(experiment == "77_GSE61220_TNF Treatment 12hrs") %>% 
  group_by(experiment) %>% 
  summarise(splicing = mean(splicing), expression = mean(expression)) %>% 
pivot_longer(cols = c("splicing", "expression"), 
             names_to = "type", values_to = "value") %>% 
  ggplot() +
  aes(x = value, color = type) +
  geom_density())

ggsave("figs/mean_overlap.pdf", plt_mean_overlap)




## Experiments ----

if(FALSE){
  ## Gene set OC ----
  overlap_coefficient <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    min_cardinality <- min(length(set1), length(set2))
    coefficient <- intersection / min_cardinality
    return(coefficient)
  }
  
  plot_GS_OC <- function(ora_all){
    
    
    ora_oc <- ora_all %>% 
      filter(
        # experiment == unique(ora_all$experiment)[4],
        experiment == paste0("77_GSE61220_TNF Treatment 12hrs", limma_suffix),
        padj_splicing < 0.05 | padj_expression < 0.05
      ) %>%
      group_by(experiment) %>% 
      mutate(
        OC_splicing = case_when(
          padj_splicing < 0.05 ~ map(
            overlapGenes_splicing,
            function(x) map_dbl(overlapGenes_expression[padj_expression < 0.05],
                                function(y) overlap_coefficient(x, y)))),
        OC_expression = case_when(
          padj_expression < 0.05 ~ map(
            overlapGenes_expression,
            function(x) map_dbl(overlapGenes_expression[padj_expression < 0.05],
                                function(y) overlap_coefficient(x, y))))
      )
    
    ora_oc <- ora_all %>% 
      filter(experiment == paste0("77_GSE61220_TNF Treatment 12hrs", limma_suffix),
             padj_splicing < 0.05) %>% 
      pull(pathway) %>% 
      map(
        function(x) map_dbl(ora_all$pathway[ora_all$padj_expression < 0.05],
                            function(y) overlap_coefficient(gene_sets[x], gene_sets[y])))
    
    quantile_overlap <- function(ora){
      
      print(ora$experiment[1])
      # Get relevant gene sets
      gs <- gene_sets[ora$pathway[ora$padj_splicing < 0.05 | ora$padj_expression < 0.05]]
      
      # Make a logical matrix of genes' presence in each gene set
      genes_tidy <- gs %>% 
        tibble::enframe(name = "gene_set", value = "gene") %>% 
        unnest(cols = c(gene))
      all_genes <- unique(genes_tidy$gene)
      
      ## Use pivot_wider to create the matrix
      gene_matrix <- genes_tidy %>%
        mutate(value = TRUE) %>%
        pivot_wider(names_from = gene_set, values_from = value, values_fill = FALSE)
      
      # Get gene sets for each analysis
      gs_expression <- ora %>% filter(padj_expression < 0.05) %>% pull(pathway)
      gs_splicing <- ora %>% filter(padj_splicing < 0.05) %>% pull(pathway)
      
      sample_expression_size <- min(length(gs_expression), length(gs_splicing))
      # Overlap coefficient for splicing vs expression and subset of expression vs expression
      oc_splicing <- proxy::simil(select(gene_matrix, all_of(gs_splicing)), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
      oc_expression <- proxy::simil(select(gene_matrix, all_of(sample(gs_expression, sample_expression_size))), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
      
      # Set diagonal to NA
      set_diag <- function(oc){
        oc_mat <- as.matrix(oc)
        name_combinations <- expand.grid(rownames(oc_mat), colnames(oc_mat))
        
        # Find matching indices based on name combinations
        matching_indices <- which(as.character(name_combinations$Var1) == as.character(name_combinations$Var2), arr.ind = TRUE)
        
        # Set the corresponding matrix elements to NA
        oc_mat[matching_indices] <- NA
        
        return(oc_mat)
      }
      oc_splicing <- set_diag(oc_splicing)
      oc_expression <- set_diag(oc_expression)
      
      # Compute 90th quantile
      quantiles_splicing <- apply(oc_splicing, 1, quantile, probs = 0.9, na.rm = TRUE)
      quantiles_expression <- apply(oc_expression, 1, quantile, probs = 0.9, na.rm = TRUE)
      
      return(tibble("splicing" = quantiles_splicing, "expression" = quantiles_expression))
    }
    
    system.time(
      mean_overlap <- map(unique(ora_all$experiment)[1], function(Experiment){
        ora_all %>% 
          filter(experiment == Experiment) %>% 
          group_by(experiment) %>% 
          summarise(results = quantile_overlap(.)) #%>% 
        # unnest_wider(results)
      }) %>% 
        bind_rows()
    )
    
    saveRDS(mean_overlap, file = "results/mean_overlap_11-20.RDS")
    
    mean_overlap <- map(list.files("results/", pattern = "mean_overlap*"), function(file){
      overlap <- readRDS(paste0("results/", file))
    }) %>% 
      bind_rows()
    
    # mean_overlap <- readRDS("results/mean_overlap_51-60.RDS")
    
    # mean_overlap <- ora_all %>% 
    #   group_by(experiment) %>% 
    #   filter(experiment %in% unique(ora_all$experiment)[1:2]) %>%
    #   summarise(results = quantile_overlap(.)) #%>%
    # unnest_wider(results)
    
    plt_mean_overlap <- mean_overlap %>% 
      # filter(experiment == "77_GSE61220_TNF Treatment 12hrs") %>% 
      group_by(experiment) %>% 
      summarise(splicing = mean(splicing), expression = mean(expression)) %>% 
      pivot_longer(cols = c("splicing", "expression"), 
                   names_to = "type", values_to = "value") %>% 
      ggplot() +
      aes(x = value, color = type) +
      geom_density()
    
    ggsave("figs/mean_overlap.pdf", plt_mean_overlap)
    
    or <- filter(ora_all, experiment == paste0("77_GSE61220_TNF Treatment 12hrs", limma_suffix))
    gs <- gene_sets[or$pathway[or$padj_splicing < 0.05 | or$padj_expression < 0.05]]
    
    genes_tidy <- gs %>% 
      tibble::enframe(name = "gene_set", value = "gene") %>% 
      unnest(cols = c(gene))
    
    all_genes <- unique(genes_tidy$gene)
    
    # Use pivot_wider to create the matrix
    gene_matrix <- genes_tidy %>%
      mutate(value = TRUE) %>%
      pivot_wider(names_from = gene_set, values_from = value, values_fill = FALSE)
    gs_expression <- or %>% filter(padj_expression < 0.05) %>% pull(pathway)
    gs_splicing <- or %>% filter(padj_splicing < 0.05) %>% pull(pathway)
    
    # similarity splicing
    system.time(
      sim_splicing <- proxy::simil(select(gene_matrix, all_of(gs_splicing)), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
    )
    simmat_splicing <- as.matrix(sim_splicing)
    name_combinations <- expand.grid(rownames(simmat_splicing), colnames(simmat_splicing))
    
    # Find matching indices based on name combinations
    matching_indices <- which(as.character(name_combinations$Var1) == as.character(name_combinations$Var2), arr.ind = TRUE)
    
    # Set the corresponding matrix elements to NA
    simmat_splicing[matching_indices] <- NA
    row_quantiles_splicing <- apply(simmat_splicing, 1, quantile, probs = 0.9, na.rm = TRUE)
    qplot(as.numeric(row_quantiles_splicing))
    
    # similarity expression
    system.time(
      sim_expression <- proxy::simil(select(gene_matrix, all_of(sample(gs_expression, 50))), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
    )
    simmat_expression <- as.matrix(sim_expression)
    name_combinations <- expand.grid(rownames(simmat_expression), colnames(simmat_expression))
    
    # Find matching indices based on name combinations
    matching_indices <- which(as.character(name_combinations$Var1) == as.character(name_combinations$Var2), arr.ind = TRUE)
    
    # Set the corresponding matrix elements to NA
    simmat_expression[matching_indices] <- NA
    row_quantiles_expression <- apply(simmat_expression, 1, quantile, probs = 0.9, na.rm = TRUE)
    qplot(as.numeric(row_quantiles_expression))
    
    # splicing + expression
    
    tibble(
      value = c(as.numeric(row_quantiles_splicing), as.numeric(row_quantiles_expression)),
      type = c(rep("splicing", length(row_quantiles_splicing)), rep("expression", length(row_quantiles_expression)))) %>% 
      ggplot() +
      aes(x = value, color = type) +
      geom_density()
    
    # similarity random
    num_gene_sets <- 50
    genes_per_set <- mean(colSums(select(gene_matrix, -gene)))
    # gs_random <- replicate(num_gene_sets, sample(gene_matrix$gene, genes_per_set), simplify = FALSE)
    
    # Add randomly generated gene sets to gene_matrix
    gene_matrix_r <- gene_matrix %>%
      tibble::add_column(!!!setNames(rep(list(FALSE), num_gene_sets), paste0("gene_set", 1:num_gene_sets)), .before = 1) %>% 
      mutate(across(starts_with("gene_set"), ~ if_else(row_number() %in% sample(nrow(gene_matrix), genes_per_set), TRUE, FALSE)))
    
    system.time(
      sim_random <- proxy::simil(select(gene_matrix_r, starts_with("gene_set")), select(gene_matrix, all_of(gs_expression)), method = "Simpson", by_rows = FALSE)
    )
    simmat_random <- as.matrix(sim_random, diag = 1)
    qplot(simmat_random)
    
    
    # Sina plot of found pathways
    plt_pathway_counts <- ora_oc %>% 
      select(OC_splicing, OC_expression) %>% 
      unnest(c(OC_splicing, OC_expression)) %>% 
      # filter(OC_expression == 1) %>% 
      # count()
      pivot_longer(cols = c(OC_splicing, OC_expression), names_to = "Analysis", values_to = "OC") %>% 
      # mutate(Analysis = factor(Analysis, levels = c("OC_splicing",
      #                                               "OC_expression")
      # ),
      # Analysis =  recode(Analysis,
      #                    OC_expression = paste0("Differential\nExpression\n", "(Median: ", median(ora_oc$OC_expression), ")"),
      #                    OC_splicing = paste0("Differential\nSplicing\n", "(Median: ", median(ora_oc$OC_splicing), ")"))
      # ) %>% 
      ggplot() +
      aes(x = OC, color = Analysis) +
      geom_density() +
      # scale_y_log10() +
      # geom_boxplot(width = 0.05) +
      # color_scale(reorder = c(2, 3, 1, 4)) + 
      labs(#x = "",
        #y = "",
        title = ora_oc$experiment[1]
      )# +
    # theme(legend.position = "none")
    
    return(plt_pathway_counts)
    
  }
}