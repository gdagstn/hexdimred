#Hexagonal binning with transparency for single cell data


`hexDimRed`


Plots a hexbin version of a dimensional reduction representing values from a series of columns of `colData`.
This is exactly what `schex`(https://github.com/SaskiaFreytag/schex) does, but with added graphical features afforded by `ggplot.multistats` (https://github.com/flying-sheep/ggplot.multistats) and a slightly more flexible approach.

1) for numeric variables, the `fill` aesthetic is given by a function that calculates summary statistics (e.g. mean, median, var) per bin. The function is specified by `stat_fun_numeric` and will be applied to each bin separately.

2) for categorical variables, the `fill` aesthetic can be given by either calculating the most represented factor level in each bin (`stat_fun_numeric` set to `"maj"`) or the proportion of a single reference level against all others (`"prop“`). This is achieved by setting all reference levels to 1 and all non-reference levels to 0, then calculating the mean. The reference level is chosen as the most abundant one in the column. In the future I may implement a choice for the reference level, although it could be tedious for the user who has to set each level independently for each column of `colData`.

3) the `alpha` aesthetic is set independently by `alpha_by` and can be the number (`length`) of observations, or any other operation performed on each bin.

A plot of `n_UMI` is always generated as the first in the list that this function returns.


Usage:

This is the traditional PCA plot as delivered by `scater`:

```r
p1 <- plotReducedDim(sce, dimred = "PCA", 
               colour_by = "individual"
              ) + labs(title = "PCA - by individual")
p2 <- plotReducedDim(sce, dimred = "PCA", 
               colour_by = "n_UMI"
              ) + labs(title = "PCA - by total UMI")
p3 <- plotReducedDim(sce, dimred = "PCA", 
               colour_by = "subset_mito_pct"
              ) + labs(title = "PCA - by % mitochondrial expression")
p4 <- plotReducedDim(sce, dimred = "PCA", 
               colour_by = "scDblFinder.class"
              ) + labs(title = "PCA - by doublet class")

gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
gridExtra::grid.arrange(grobs = list(p3, p4), ncol = 2)
```

![sc_pca_1](/figures/sc_pca1.png)
![sc_pca_2](/figures/sc_pca2.png)


It is hard to make out the relative contribution of each level because of how dots are superimposed.


Using `hexDimRed`:


```r
sce_hexdim <- hexDimRed(sce, 
                       stat_fun_categorical = "prop", 
                       other_fields = c("individual", 
                                        "subset_mito_pct", 
                                        "scDblFinder.class"
                                       )
                      )

sce_hexdim2 <- hexDimRed(sce, 
                        stat_fun_categorical = "maj", 
                        other_fields = c("individual", 
                                         "scDblFinder.class",
                                         "batch"
                                        )
                       )


gridExtra::grid.arrange(grobs = ct_hexdim[c(1,4)], ncol = 2)
gridExtra::grid.arrange(grobs = c(ct_hexdim[4], ct_hexdim2[2]), ncol = 2)
gridExtra::grid.arrange(grobs = ct_hexdim2[c(3,4)], ncol = 2)
```

![pca_1](/figures/pca1.png)
![pca_2](/figures/pca2.png)
![pca_3](/figures/pca3.png)

---
