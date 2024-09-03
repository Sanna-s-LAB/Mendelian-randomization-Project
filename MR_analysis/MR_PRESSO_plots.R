library(dplyr)
require(TwoSampleMR)
require(ieugwasr)
require(MRPRESSO)
library(officer)
library(ggplot2)
# Function for IVW before and after MR-PRESSO

mr_scatter_plot_IVW <- function(mr_results,mr_results2, dat, dat_kept, dat_removed)
{
  mr_results <- subset(mr_results, method == "Inverse variance weighted")
  mr_results2 <- subset(mr_results2, method == "Inverse variance weighted")
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  { 
    # DAT REMOVED
    dat_removed <- plyr::mutate(dat_removed)
    if(nrow(dat_removed) < 1 | sum(dat_removed$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    dat_removed <- subset(dat_removed, mr_keep)
    index <- dat_removed$beta.exposure < 0
    dat_removed$beta.exposure[index] <- dat_removed$beta.exposure[index] * -1
    dat_removed$beta.outcome[index] <- dat_removed$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == dat_removed$id.exposure[1] & id.outcome == dat_removed$id.outcome[1])
    
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(dat_removed$beta.exposure, dat_removed$beta.outcome, dat_removed$se.exposure, dat_removed$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(dat_removed$beta.exposure, dat_removed$beta.outcome, dat_removed$se.exposure, dat_removed$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    # DAT KEPT
    dat_kept <- plyr::mutate(dat_kept)
    if(nrow(dat_kept) < 1 | sum(dat_kept$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    dat_kept <- subset(dat_kept, mr_keep)
    index <- dat_kept$beta.exposure < 0
    dat_kept$beta.exposure[index] <- dat_kept$beta.exposure[index] * -1
    dat_kept$beta.outcome[index] <- dat_kept$beta.outcome[index] * -1
    mrres2 <- subset(mr_results2, id.exposure == dat_kept$id.exposure[1] & id.outcome == dat_kept$id.outcome[1])
    
    mrres2$a <- 0
    if("MR Egger" %in% mrres2$method)
    {
      temp <- mr_egger_regression(dat_kept$beta.exposure, dat_kept$beta.outcome, dat_kept$se.exposure, dat_kept$se.outcome, default_parameters())
      mrres2$a[mrres2$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres2$method)
    {
      temp <- mr_egger_regression_bootstrap(dat_kept$beta.exposure, dat_kept$beta.outcome, dat_kept$se.exposure, dat_kept$se.outcome, default_parameters())
      mrres2$a[mrres2$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    # data
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data = dat_kept, ggplot2::aes(x = beta.exposure, y = beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
      ggplot2::geom_abline(data = mrres2, ggplot2::aes(intercept = a, slope = b, colour = "MR-PRESSO (outlier removed)"), show.legend = TRUE) +
      ggplot2::scale_colour_manual(values = c("#a6cee3","#e31a1c", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour = "MR Test", x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position = "top", legend.direction = "vertical") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
      ggplot2::geom_errorbar(data = dat_removed, ggplot2::aes(x = beta.exposure, ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
      ggplot2::geom_errorbarh(data = dat_removed, ggplot2::aes(y = beta.outcome, xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
      ggplot2::geom_point(data = dat_removed, ggplot2::aes(x = beta.exposure, y = beta.outcome), color = "red")
  })
  mrres
}


# Function for scatter plot with different colour for outliers

mr_scatter_plot_col <- function(mr_results, dat, dat_kept, dat_removed)
{
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  { 
    # DAT REMOVED
    dat_removed <- plyr::mutate(dat_removed)
    if(nrow(dat_removed) < 1 | sum(dat_removed$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    dat_removed <- subset(dat_removed, mr_keep)
    index <- dat_removed$beta.exposure < 0
    dat_removed$beta.exposure[index] <- dat_removed$beta.exposure[index] * -1
    dat_removed$beta.outcome[index] <- dat_removed$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == dat_removed$id.exposure[1] & id.outcome == dat_removed$id.outcome[1])
    
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(dat_removed$beta.exposure, dat_removed$beta.outcome, dat_removed$se.exposure, dat_removed$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(dat_removed$beta.exposure, dat_removed$beta.outcome, dat_removed$se.exposure, dat_removed$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    # DAT KEPT
    dat_kept <- plyr::mutate(dat_kept)
    if(nrow(dat_kept) < 1 | sum(dat_kept$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    dat_kept <- subset(dat_kept, mr_keep)
    index <- dat_kept$beta.exposure < 0
    dat_kept$beta.exposure[index] <- dat_kept$beta.exposure[index] * -1
    dat_kept$beta.outcome[index] <- dat_kept$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == dat_kept$id.exposure[1] & id.outcome == dat_kept$id.outcome[1])
    
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(dat_kept$beta.exposure, dat_kept$beta.outcome, dat_kept$se.exposure, dat_kept$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(dat_kept$beta.exposure, dat_kept$beta.outcome, dat_kept$se.exposure, dat_kept$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data = dat_kept, ggplot2::aes(x = beta.exposure, y = beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
      ggplot2::scale_colour_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour = "MR Test", x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position = "top", legend.direction = "vertical") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
      ggplot2::geom_errorbar(data = dat_removed, ggplot2::aes(x = beta.exposure, ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
      ggplot2::geom_errorbarh(data = dat_removed, ggplot2::aes(y = beta.outcome, xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
      ggplot2::geom_point(data = dat_removed, ggplot2::aes(x = beta.exposure, y = beta.outcome), color = "red")
    
  })
  mrres
}

blank_plot <- function(message)
{
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + 
    ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + 
    ggplot2::labs(x=NULL,y=NULL) + 
    ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}