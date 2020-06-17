#' Plot dimensionality reduction as hexagon bins with summary stats
#'
#' @param o a SingleCellExperiment class object
#' @param dimred character, name of the slot (called by \code{reducedDim}) containing the dimensionality reduction in o. Default is "PCA"
#' @param dims_use numeric, vector of 2 coordinates of the dimensionality reduction to be plotted. Default is 1, 2
#' @param alpha_by function used to set the alpha value of each bin. Default is \code{length}
#' @param stat_fun_numeric character, name of the function used to summarize the numeric value in each bin. Default is "median"
#' @param stat_fun_categorical character, name of the function used to summarize the factor value in each bin. If \code{prop}, 
#' it calculates the proportion of the most abundant factor against the others (useful for factors with 2 levels). If \code{maj},
#' it calculates which of the factor levels is the most abundant in each bin (useful for factors with > 2 levels).
#' @param nbins numeric, number of bins to group and plot the data. Default is 40
#' @param main_field character, name of the field to be plotted first. Defaults to "n_UMI". If NULL, gives an "unlabeled" density plot
#' @param nbins number of bins to group and plot the data. Default is 40
#' @param other_fields character, vector of column names from colData to be used as continuous or grouping variables.
#' @param plot.genes character, vector of gene names to be plotted
#' @param gene.counts character, assay used for gene plots
#' @param custom_cat_cols character, a vector of colours to be used in categorical variables. For the colour assignment to work, colours
#' must be' ordered in the same ordered as the levels of the factor that is being plotted. If NULL, defaults to the standard ggplot2 colour scheme
#' @param custom_cont_cols character, a vector of colours to be used in continuous variables. If NULL, defaults to viridis (option D)
#' @param subset list, colData field and levels to subset the object. Takes the name of the list element as colData field and 
#' its contents as levels of the field.
#' @param dark logical, use a dark theme? Default is FALSE, default theme is theme_classic()
#' @param labels character, add labels to the plot? Only used when viewing categorical data with "maj". Default is "none".
#' options include "sticker" for labels with an opaque rounded box and "text" to have text-only labels
#' @return a list of ggplot2 hex-bin plots.

hexDimRed <- function(o,
                      dimred = "PCA",
                      dims_use = c(1,2),
                      alpha_by = "length",
                      stat_fun_numeric = "median",
                      stat_fun_categorical = "prop",
                      nbins = 40,
                      main_field = "n_UMI",
                      other_fields = NULL,
                      plot.genes = NULL,
                      gene.counts = "logcounts",
                      custom_cat_cols = NULL,
                      custom_cont_cols = NULL,
                      subset = NULL,
                      dark = FALSE,
                      labels = "none"
                      ){

require(ggplot2)
require(ggplot.multistats)

# Sanity checks
if (class(o) != "SingleCellExperiment") stop("o must be an object of class SingleCellExperiment")

if (dimred %in% reducedDimNames(o) == FALSE) stop(paste0("Dimensional reduction ", dimred, " not found."))

if (length(dims_use) > 2) {
  print("More than 2 dimensions provided, will use the first 2")
  dims_use <- dims_use[1:2]
} 

if(!dark) theme =  theme = theme_classic()

if(dark) {

require(ggdark)

  theme = ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 14) + 
  theme(plot.title = element_text(family = "Fira Sans Condensed"),
        plot.background = element_rect(fill = "grey10"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey30", size = 0.2),
        panel.grid.minor = element_line(color = "grey30", size = 0.2),
        legend.background = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank())
}

if(!labels %in% c("none", "sticker", "text")) {
    print("labels argument invalid, setting to none.")
    labels = "none"
    }

# Subsetting the object before plotting

if(!is.null(subset)) {
    if(length(subset) > 1) stop("Can only subset using 1 field. Subset the object outside of the function for more complex subsetting")
    sampled_obj <- o[, which(colData(o)[,names(subset)] %in% unlist(subset))]
    attr(reducedDim(sampled_obj, "PCA"), "percentVar") <- attr(reducedDim(o, "PCA"), "percentVar")
    o <- sampled_obj
}

if(length(dims_use) <= 1) stop("Need 2 dimensions to plot")

if(!is.null(plot.genes) & length(intersect(plot.genes, rownames(o))) == 0) {
  plot.genes = NULL
  print("Selected genes were not found in the object, setting them to NULL")
} else if(!is.null(plot.genes) & !all(plot.genes %in% rownames(o))) {
  not.found <-  setdiff(plot.genes, rownames(o))
  plot.genes <- intersect(plot.genes, rownames(o))
  print(paste0("The following genes were not found in the object: ", paste(not.found, collapse = ", ")))
}

# Wrangling

tmp <- as.data.frame(cbind(colData(o), reducedDim(o, dimred)[,dims_use]))
coln_xy = colnames(tmp)[(ncol(tmp) - 1) : (ncol(tmp))]

colnames(tmp)[(ncol(tmp) - 1) : (ncol(tmp))] <- c("Dim1", "Dim2")

# Special case for PCA: adding % variance explained by each PC

if (dimred == "PCA" | dimred == "pca") {
  pctvar <- round(attr(reducedDim(o, "PCA"), "percentVar"), digits = 2)
}

# First plot: if main_field is null, colour as gray

if(is.null(main_field)) {
p1 <- ggplot(data = tmp, aes(x = Dim1, y = Dim2)) +
  stat_summaries_hex(
    aes(z = 1, alpha = stat(log10(n))),
    funs = c(n = alpha_by),
    bins = nbins
  ) +
  scale_fill_manual(values = "gray") +
  theme

} else {

# Otherwise use main_field for colouring

p1 <- ggplot(data = tmp, aes(x = Dim1, y = Dim2)) +
  stat_summaries_hex(
    aes(z = !!ensym(main_field), fill = stat(get(stat_fun_numeric)), alpha = stat(log10(n))),
    funs = c(stat_fun_numeric, n = alpha_by),
    bins = nbins
  ) 
     if(is.null(custom_cont_cols)) {
      p1 <- p1 + scale_fill_viridis(option = "D") + theme
    } else {
      p1 <- p1 + scale_fill_gradient(colours = custom_cont_cols) + theme
    }
}

# Add special axes/labels for PCA

if(is.null(main_field)) {
  if (dimred == "PCA" | dimred == "pca") {
  p1 <- p1 +
  labs(title = paste0(dimred, " by ", main_field), 
    x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
    y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
    alpha = paste0(alpha_by, " of observations in bin (log10)")
   ) } else {

   p1 <- p1 +
  labs(title = paste0(dimred, " by ", main_field), 
   x = coln_xy[1],
   y = coln_xy[2],
   alpha = paste0(alpha_by, " of observations in bin (log10)")
   )
  }
} else {
 if (dimred == "PCA" | dimred == "pca") {
  p1 <- p1 +
  labs(title = paste0(dimred, " by ", main_field), 
    x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
    y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
    fill = paste0(stat_fun_numeric, " ", main_field, " per bin"),
    alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) } else {

    p1 <- p1 +
  labs(title = paste0(dimred, " by ", main_field), 
   x = coln_xy[1],
   y = coln_xy[2],
   fill = paste0(stat_fun_numeric, " ", main_field, " per bin"),
   alpha = paste0(alpha_by, " of observations in bin (log10)")
   )
  }
}

# Other fields

if (!is.null(other_fields)) {

addedplotlist <- list()

for (f in other_fields) {

  if (!is.numeric(tmp[,f]) & !is.factor(tmp[,f]) & class(tmp[,f]) == "character") tmp[,f] <- factor(tmp[,f]) 

  if (is.numeric(tmp[,f])) {

    addedplotlist[[f]] <-  ggplot(data = tmp, aes(x = Dim1, y = Dim2)) +
      stat_summaries_hex(
        aes(z = !!ensym(f), fill = stat(get(stat_fun_numeric)), alpha = stat(log10(n))),
        funs = c(stat_fun_numeric, n = alpha_by),
        bins = nbins
      )
     
     if(is.null(custom_cont_cols)) {
      addedplotlist[[f]] <- addedplotlist[[f]] + scale_fill_viridis(option = "D") + theme
    } else {
      addedplotlist[[f]] <- addedplotlist[[f]] + scale_fill_gradient(colours = custom_cont_cols) + theme
    }
      

      if (dimred == "PCA") {

addedplotlist[[f]]  <- addedplotlist[[f]]  +
labs(title = paste0(dimred, " by ", f), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0(stat_fun_numeric, " ", f, " in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

addedplotlist[[f]]  <- addedplotlist[[f]]  +
labs(title = paste0(dimred, " by ", f), 
      x = coln_xy[1],
      y = coln_xy[2],
      fill = paste0(stat_fun_numeric, " ", f, " in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    )

      }

    } else if (is.factor(tmp[,f]) & stat_fun_categorical == "prop"){
      tmpp <- tmp[,c(f, "Dim1", "Dim2")]
      topl <- names(which.max(table(tmpp[,f])))
      tmpp$toplevel <- as.character(tmpp[,f])
      tmpp$toplevel[tmpp$toplevel == topl] = 1
      tmpp$toplevel[tmpp$toplevel != 1] = 0
      tmpp$toplevel <- as.numeric(tmpp$toplevel)

   addedplotlist[[f]] <- ggplot(data = tmpp, aes(x = Dim1, y = Dim2, group = 1)) +
      stat_summaries_hex(
      aes(z = toplevel, fill = stat(mean), alpha = stat(log10(n))),
      funs = c("mean", n = "length"),
      bins = nbins
       )      

      if(is.null(custom_cont_cols)) {
      addedplotlist[[f]] <- addedplotlist[[f]] + scale_fill_viridis(option = "D") + theme
    } else if(!is.null(custom_cont_cols)) {
      addedplotlist[[f]] <- addedplotlist[[f]] + scale_fill_gradient(colours = custom_cont_cols) + theme
    }

      if (dimred == "PCA") {

  addedplotlist[[f]]  <- addedplotlist[[f]]  +
  labs(title = paste0(dimred, " by ", f), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0("Proportion of ", topl, " in bin"),
        alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

  addedplotlist[[f]]  <- addedplotlist[[f]]  +
  labs(title = paste0(dimred, " by ", f), 
      x = coln_xy[1],
      y = coln_xy[2],
      fill = paste0("Proportion of ", topl, " in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
      )
        }

    } else if (is.factor(tmp[,f]) & stat_fun_categorical == "maj"){
        tmpm <- tmp[,c(f, "Dim1", "Dim2")]
   addedplotlist[[f]] <- ggplot(data = tmpm, aes(x = Dim1, y = Dim2, group = 1)) +
      stat_summaries_hex(
      aes(z = !!ensym(f), fill = factor(stat(most), levels = levels(tmpm[,f])), alpha = stat(log10(n))),
      funs = c(most = function(s) names(which.max(summary(s))), n = "length"),
      bins = nbins
       )     
          
if(labels != "none") {
        require(ggrepel)
        factor.labels = unique(as.character(tmpm[,f]))
        centers.x <- unlist(sapply(factor.labels, function(x) mean(tmpm[tmpm[,f] == x, "Dim1"])))
        centers.y <- unlist(sapply(factor.labels, function(x) mean(tmpm[tmpm[,f] == x, "Dim2"])))       
        labels.df <- as.data.frame(cbind(factor.labels, centers.x, centers.y), stringsAsFactors = FALSE)                       
        colnames(labels.df) <- c("label", "x", "y")
                                   
     if(labels == "sticker")  {                        
        addedplotlist[[f]] <- addedplotlist[[f]] + 
        geom_label_repel(data = labels.df, 
                         mapping = aes(x = as.numeric(x), y = as.numeric(y), label = label),
                         col = c("black", "white")[as.numeric(dark)+1],
                         fill = c("white", "black")[as.numeric(dark)+1],
                         inherit.aes = FALSE)
             
      } else if(labels == "text") {
        addedplotlist[[f]] <- addedplotlist[[f]] + 
        geom_text_repel(data = labels.df, 
                        mapping = aes(x = as.numeric(x), y = as.numeric(y), label = label), 
                        col = c("black", "white")[as.numeric(dark)+1],
                        inherit.aes = FALSE)
         }
      }   
                                   
      if(is.null(custom_cat_cols)) {
      addedplotlist[[f]] <- addedplotlist[[f]] + theme
    } else if(!is.null(custom_cat_cols)) {
      addedplotlist[[f]] <- addedplotlist[[f]] + scale_fill_manual(values = custom_cat_cols) + theme
    }

      if (dimred == "PCA") {

  addedplotlist[[f]]  <- addedplotlist[[f]]  +
  labs(title = paste0(dimred, " by ", f), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0(f, " most represented in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

  addedplotlist[[f]]  <- addedplotlist[[f]]  +
  labs(title = paste0(dimred, " by ", f), 
       x = coln_xy[1],
       y = coln_xy[2],
       fill = paste0(f, " most represented in bin"),
       alpha = paste0(alpha_by, " of observations in bin (log10)")
      )
     }
    }
  }
}

# Final assembly conditional on selected fields

  if(!is.null(main_field)) {
    finallist <- c(list("n_UMI" = p1))
  } else  {
    finallist <- c(list("cells_per_bin" = p1))
  }

  if(!is.null(other_fields)) finallist <- c(finallist, addedplotlist)

  if(!is.null(plot.genes)) finallist <- c(finallist, geneplotlist)

  return(finallist)

}
