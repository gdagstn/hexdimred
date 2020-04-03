
hexDimRed <- function(o,
                      dimred = "PCA",
                      dims_use = c(1,2),
                      alpha_by = "length",
                      stat_fun_numeric = "median",
                      stat_fun_categorical = "prop",
                      nbins = 40,
                      other_fields = NULL
                      ){

require(ggplot2)
require(ggplot.multistats)


if (class(o) != "SingleCellExperiment") stop("o must be an object of class SingleCellExperiment")
if (is.null(reducedDim(o, dimred))) stop(paste0("Dimensional reduction ", dimred, " not found."))
if (length(dims_use) > 2) {
  print("More than 2 dimensions provided, will use the first 2")
  dims_use <- dims_use[1:2]
} 

if(length(dims_use) <= 1) stop("Need 2 dimensions to plot")

tmp <- as.data.frame(cbind(colData(o), reducedDim(o, dimred)[,dims_use]))
coln_xy = colnames(tmp)[(ncol(tmp) - 1) : (ncol(tmp))]

colnames(tmp)[(ncol(tmp) - 1) : (ncol(tmp))] <- c("Dim1", "Dim2")

if (dimred == "PCA") {
  pctvar <- round(attr(reducedDim(o, "PCA"), "percentVar"), digits = 2)
}

p1 <- ggplot(data = tmp, aes(x = Dim1, y = Dim2)) +
  stat_summaries_hex(
    aes(z = n_UMI, fill = stat(get(stat_fun_numeric)), alpha = stat(log10(n))),
    funs = c(stat_fun_numeric, n = alpha_by),
    bins = nbins
  ) +
  scale_fill_viridis(option = "D") +
  theme_classic()

if (dimred == "PCA") {
p1 <- p1 +
labs(title = paste0(dimred, " by total UMI"), 
  x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
  y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
  fill = paste0(stat_fun_numeric, " total UMI per bin"),
  alpha = paste0(alpha_by, " of observations in bin (log10)")
  ) } else {

  p1 <- p1 +
labs(title = paste0(dimred, " by total UMI"), 
  x = coln_xy[1],
  y = coln_xy[2],
  fill = paste0(stat_fun_numeric, " total UMI per bin"),
  alpha = paste0(alpha_by, " of observations in bin (log10)")
  )
}



if (!is.null(other_fields)) {

addedplotlist <- list()

for (i in other_fields) {

  if (!is.numeric(tmp[,i]) & !is.factor(tmp[,i])) tmp[,i] <- factor(tmp[,i]) 

  if (is.numeric(tmp[,i])) {

    addedplotlist[[i]] <-  ggplot(data = tmp, aes(x = Dim1, y = Dim2)) +
      stat_summaries_hex(
        aes(z = !!ensym(i), fill = stat(get(stat_fun_numeric)), alpha = stat(log10(n))),
        funs = c(stat_fun_numeric, n = alpha_by),
        bins = nbins
      ) +
      scale_fill_viridis(option = "D") +
      theme_classic()

      if (dimred == "PCA") {

addedplotlist[[i]]  <- addedplotlist[[i]]  +
labs(title = paste0(dimred, " by ", i), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0(stat_fun_numeric, " ", i, " in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

addedplotlist[[i]]  <- addedplotlist[[i]]  +
labs(title = paste0(dimred, " by ", i), 
      x = coln_xy[1],
      y = coln_xy[2],
      fill = paste0(stat_fun_numeric, " ", i, " in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    )

      }

    } else if (is.factor(tmp[,i]) & stat_fun_categorical == "prop"){

      topl <- names(which.max(table(tmp[,i])))
      tmp$toplevel <- as.character(tmp[,i])
      tmp$toplevel[tmp$toplevel == topl] = 1
      tmp$toplevel[tmp$toplevel != 1] = 0
      tmp$toplevel <- as.numeric(tmp$toplevel)

   addedplotlist[[i]] <- ggplot(data = tmp, aes(x = Dim1, y = Dim2, group = 1)) +
      stat_summaries_hex(
      aes(z = toplevel, fill = stat(mean), alpha = stat(log10(n))),
      funs = c("mean", n = "length"),
      bins = nbins
       ) +
      scale_fill_viridis(option = "D") +
      theme_classic()

      if (dimred == "PCA") {

  addedplotlist[[i]]  <- addedplotlist[[i]]  +
  labs(title = paste0(dimred, " by ", i), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0("Proportion of ", topl, " in bin"),
        alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

  addedplotlist[[i]]  <- addedplotlist[[i]]  +
  labs(title = paste0(dimred, " by ", i), 
      x = coln_xy[1],
      y = coln_xy[2],
      fill = paste0("Proportion of ", topl, " in bin"),
        alpha = paste0(alpha_by, " of observations in bin (log10)")
      )
        }

    } else if (is.factor(tmp[,i]) & stat_fun_categorical == "maj"){

   addedplotlist[[i]] <- ggplot(data = tmp, aes(x = Dim1, y = Dim2, group = 1)) +
      stat_summaries_hex(
      aes(z = !!ensym(i), fill = stat(most), alpha = stat(log10(n))),
      funs = c(most = function(s) names(which.max(summary(s))), n = "length"),
      bins = nbins
       ) +
      theme_classic()

      if (dimred == "PCA") {

  addedplotlist[[i]]  <- addedplotlist[[i]]  +
  labs(title = paste0(dimred, " by ", i), 
      x = paste0(coln_xy[1], " - ", pctvar[dims_use[1]], "% variance"),
      y = paste0(coln_xy[2], " - ", pctvar[dims_use[2]], "% variance"),
      fill = paste0(i, " most represented in bin"),
      alpha = paste0(alpha_by, " of observations in bin (log10)")
    ) 

    } else {

  addedplotlist[[i]]  <- addedplotlist[[i]]  +
  labs(title = paste0(dimred, " by ", i), 
       x = coln_xy[1],
       y = coln_xy[2],
       fill = paste0(i, " most represented in bin"),
       alpha = paste0(alpha_by, " of observations in bin (log10)")
      )
        }

    }

 }

    finallist <- c(list("n_UMI" = p1), addedplotlist)
    names(finallist) <- c("n_UMI", other_fields)
  
    return(finallist)

} else if (is.null(other_fields)) {

    return(p1)
  }
}