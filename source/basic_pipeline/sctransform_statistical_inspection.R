# saveRDS(lung_tme_exp1, file='Q:\\BIML2020\\data\\mouse_lung_SCLC\\lung_tme_exp1_seuratset.rds')

# Inspect data
exprs_mat=lung_tme_exp1@assays$RNA@counts
gene_attr <- data.frame(mean = rowMeans(exprs_mat), 
                        detection_rate = rowMeans(exprs_mat > 0),
                        var = apply(exprs_mat, 1, var))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(exprs_mat)
cell_attr <- data.frame(n_umi = colSums(exprs_mat),
                        n_gene = colSums(exprs_mat > 0))
rownames(cell_attr) <- colnames(exprs_mat)

# We will now calculate some properties and visually inspect the data. 
# Our main interest is in the general trends not in individual outliers. 
# Neither genes nor cells that stand out are important at this step, 
# but we focus on the global trends.
# Derive gene and cell attributes from the UMI matrix.
ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color='red')+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

# add the expected detection rate under Poisson model
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

# Mean-detection-rate relationship
# In line with the previous plot, 
# we see a lower than expected detection rate in the medium expression range. 
# However, for the highly expressed genes, the rate is at or 
# very close to 1.0 suggesting that there is no zero-inflation in the counts
# for those genes and that zero-inflation is a result of overdispersion, rather than an independent systematic bias.

ggplot(cell_attr, aes(n_umi, n_gene)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

  
# Estimate model parameters and transform data
# The vst function estimates model parameters and performs the variance stabilizing transformation. 
# Here we use the log10 of the total UMI counts of a cell as variable for sequencing depth for each cell.
# After data transformation we plot the model parameters as a function of gene mean (geometric mean).

# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

set.seed(44)
vst_out <- sctransform::vst(exprs_mat, latent_var = c('log_umi'), return_gene_attr = TRUE, return_cell_attr = TRUE, show_progress = FALSE)
sctransform::plot_model_pars(vst_out)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

sctransform::plot_model(vst_out, exprs_mat, c('S100a8','S100a9'), plot_residual = TRUE)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggplot(vst_out$gene_attr, aes(residual_mean)) + geom_histogram(binwidth=0.01)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggplot(vst_out$gene_attr, aes(residual_variance)) + geom_histogram(binwidth=0.1) + geom_vline(xintercept=1, color='red') + xlim(0, 10) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) + geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
