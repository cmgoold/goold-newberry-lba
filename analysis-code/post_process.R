############################################################################
# R post-processing script for Goold & Newberry (2021):
#     - Predicting individual shelter dog behaviour after adoption using 
#       longitudinal behavioural assessment: a hierarchical Bayesian approach

# Copyright(C) Conor Goold  2021
# goold.conor@gmail.com
############################################################################

path <- "~/Dropbox/PhD/PhD_NMBU/PaperIV/goold-newberry-lba/"
figure_path <- paste0(path, "paper/figures/")
source(paste0(path, "analysis-code/helper_functions.R"))
setwd(path)

get_needed_packages()

# load the posterior distribution
chain_1 <- rstan::read_stan_csv(paste0(path, "analysis-code/jhb_dog_samples1.csv"))
chain_2 <- rstan::read_stan_csv(paste0(path, "analysis_code/jhb_dog_samples2.csv"))

posterior_draws_1 <- posterior::as_draws_df(x = chain_1)
posterior_draws_2 <- posterior::as_draws_df(x = chain_2)

# chains combined
draws <- do.call("rbind.data.frame", list(posterior_draws_1, posterior_draws_2))

# ensure that analysis_script has been ran in this session so that 
#   we have access to d_s_stan, d_a_stan, X_s_stan, X_a_stan

# get the random effects for later plotting using a convenience function
random_effect_list <- get_marginal_samples()

#---------------------------------------------------------------------------------------------
# Key shelter results

# build data frame of primary shelter results
shelter_main_results_df <- data.frame(
  parameter = c("alpha", "beta", "alpha_miss", "beta_miss", "sigma", "kappa[1]")
)

shelter_main_results_df$posterior_mean <- apply(draws[,shelter_main_results_df$parameter], 2, mean)
shelter_main_results_df$hdi_2.5 <- apply(draws[,shelter_main_results_df$parameter], 2, hdi)[1,]
shelter_main_results_df$hdi_97.5 <- apply(draws[,shelter_main_results_df$parameter], 2, hdi)[2,]
print(shelter_main_results_df)

# build data frame of dog-level predictor variables
get_betas <- draws[, grep("Beta", colnames(draws))]
colnames(get_betas) <- c(colnames(X_s), paste0(colnames(X_s), ".miss"))

shelter_dl_predictors_df <- data.frame(
  parameter = c("length_of_stay", "weight", "age", "sex_male.female", 
                "neutered_no.yes", "neutered_no.onsite", "neutered_onsite.yes", 
                "source_gift.stray", "source_gift.return", "source_return.stray",
                paste0(
                  c("length_of_stay", "weight", "age", "sex_male.female", 
                  "neutered_no.yes", "neutered_no.onsite", "neutered_onsite.yes", 
                  "source_gift.stray", "source_gift.return", "source_return.stray"
                ), ".miss")), 
  
  posterior_mean = c(
    mean(get_betas$length_of_stay_Z), 
    mean(get_betas$latest_weight_Z), 
    mean(get_betas$estimated_age_at_departure_Z), 
    mean(get_betas$sex_MALE - get_betas$sex_MALE*-1),
    mean(get_betas$neutered_NO - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1),
    mean(get_betas$neutered_NO - get_betas$neutered_ONSITE),
    mean(get_betas$neutered_ONSITE - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1),
    mean(get_betas$source_GIFT - (get_betas$source_GIFT + get_betas$source_RETURN)*-1), 
    mean(get_betas$source_GIFT - get_betas$source_RETURN), 
    mean(get_betas$source_RETURN - (get_betas$source_GIFT + get_betas$source_RETURN)*-1),
    # missing
    mean(get_betas$length_of_stay_Z.miss), 
    mean(get_betas$latest_weight_Z.miss), 
    mean(get_betas$estimated_age_at_departure_Z.miss), 
    mean(get_betas$sex_MALE.miss - get_betas$sex_MALE.miss*-1),
    mean(get_betas$neutered_NO.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1),
    mean(get_betas$neutered_NO.miss - get_betas$neutered_ONSITE.miss),
    mean(get_betas$neutered_ONSITE.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1),
    mean(get_betas$source_GIFT.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1), 
    mean(get_betas$source_GIFT.miss - get_betas$source_RETURN.miss), 
    mean(get_betas$source_RETURN.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1)
  ),
  
  hdi_2.5 = c(
    hdi(get_betas$length_of_stay_Z)[1], 
    hdi(get_betas$latest_weight_Z)[1], 
    hdi(get_betas$estimated_age_at_departure_Z)[1], 
    hdi(get_betas$sex_MALE - get_betas$sex_MALE*-1)[1],
    hdi(get_betas$neutered_NO - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1)[1],
    hdi(get_betas$neutered_NO - get_betas$neutered_ONSITE)[1],
    hdi(get_betas$neutered_ONSITE - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1)[1],
    hdi(get_betas$source_GIFT - (get_betas$source_GIFT + get_betas$source_RETURN)*-1)[1], 
    hdi(get_betas$source_GIFT - get_betas$source_RETURN)[1], 
    hdi(get_betas$source_RETURN - (get_betas$source_GIFT + get_betas$source_RETURN)*-1)[1],
    # missing
    hdi(get_betas$length_of_stay_Z.miss)[1], 
    hdi(get_betas$latest_weight_Z.miss)[1], 
    hdi(get_betas$estimated_age_at_departure_Z.miss)[1], 
    hdi(get_betas$sex_MALE.miss - get_betas$sex_MALE.miss*-1)[1],
    hdi(get_betas$neutered_NO.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1)[1],
    hdi(get_betas$neutered_NO.miss - get_betas$neutered_ONSITE.miss)[1],
    hdi(get_betas$neutered_ONSITE.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1)[1],
    hdi(get_betas$source_GIFT.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1)[1], 
    hdi(get_betas$source_GIFT.miss - get_betas$source_RETURN.miss)[1], 
    hdi(get_betas$source_RETURN.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1)[1]
  ),
  
  hdi_97.5 = c(
    hdi(get_betas$length_of_stay_Z)[2], 
    hdi(get_betas$latest_weight_Z)[2], 
    hdi(get_betas$estimated_age_at_departure_Z)[2], 
    hdi(get_betas$sex_MALE - get_betas$sex_MALE*-1)[2],
    hdi(get_betas$neutered_NO - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1)[2],
    hdi(get_betas$neutered_NO - get_betas$neutered_ONSITE)[2],
    hdi(get_betas$neutered_ONSITE - (get_betas$neutered_NO + get_betas$neutered_ONSITE)*-1)[2],
    hdi(get_betas$source_GIFT - (get_betas$source_GIFT + get_betas$source_RETURN)*-1)[2], 
    hdi(get_betas$source_GIFT - get_betas$source_RETURN)[2], 
    hdi(get_betas$source_RETURN - (get_betas$source_GIFT + get_betas$source_RETURN)*-1)[2],
    # missing
    hdi(get_betas$length_of_stay_Z.miss)[2], 
    hdi(get_betas$latest_weight_Z.miss)[2], 
    hdi(get_betas$estimated_age_at_departure_Z.miss)[2], 
    hdi(get_betas$sex_MALE.miss - get_betas$sex_MALE.miss*-1)[2],
    hdi(get_betas$neutered_NO.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1)[2],
    hdi(get_betas$neutered_NO.miss - get_betas$neutered_ONSITE.miss)[2],
    hdi(get_betas$neutered_ONSITE.miss - (get_betas$neutered_NO.miss + get_betas$neutered_ONSITE.miss)*-1)[2],
    hdi(get_betas$source_GIFT.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1)[2], 
    hdi(get_betas$source_GIFT.miss - get_betas$source_RETURN.miss)[2], 
    hdi(get_betas$source_RETURN.miss - (get_betas$source_GIFT.miss + get_betas$source_RETURN.miss)*-1)[2]
  )
)

print(shelter_dl_predictors_df)

#---------------------------------------------------------------------------------------------
# FIGURE 1: plot code colour probabilities across time at shelter

# days to plot
min_max_days <- c(0, table_1$length_of_stay[[1]] + table_1$length_of_stay[[2]])
# grid of days
days_grid <- seq(min_max_days[1], min_max_days[2], length.out = 10)
# standardised grid of days (use values to standardise from analysis_script)
days_grid_Z <- ( days_grid - table_1$length_of_stay[[1]] ) / sd(d_s$day_since_arrival)

shelter_pred_list <- rep(list(list(list())), length(days_grid_Z))
shelter_missing_pred_list <- rep(list(list(list())), length(days_grid_Z))

for(i in 1:length(days_grid_Z)){#for each day
  message(paste0("Day ", i))
  for(j in 1:length(unique(d_s_stan$dog_id))){#for each dog
    
    # container
    out_jg <- rep(list(list()), length(unique(d_s_stan$context_id)))
    out_miss_jg <- rep(list(list()), length(unique(d_s_stan$context_id)))
    
    for(g in 1:length(unique(d_s_stan$context_id))){#for each context
      
      out_jg[[g]] <- draws$alpha + # mean
        random_effect_list$r_s_j_intercept[,j] + # dog
        random_effect_list$r_s_g_intercept[,g] + # context
        random_effect_list$r_s_jg_intercept[,make_cc_interaction(j,g)] + # dog x context
        (draws$beta + 
           random_effect_list$r_s_j_slope[,j] + # dog
           random_effect_list$r_s_g_slope[,g] + # context
           random_effect_list$r_s_jg_slope[,make_cc_interaction(j,g)]) * days_grid_Z[i]
      
      out_miss_jg[[g]] <- draws$alpha_miss + # mean
        random_effect_list$r_s_j_missing[,j] + # dog
        random_effect_list$r_s_g_missing[,g] + # context
        random_effect_list$r_s_jg_missing[,make_cc_interaction(j,g)] + # dog x context
        draws$beta_miss * days_grid_Z[i]
      
      }#g
    
    shelter_pred_list[[i]][[j]] <- matrix(unlist(out_jg), ncol=8)
    shelter_missing_pred_list[[i]][[j]] <- matrix(unlist(out_miss_jg), ncol=8)
    
    }
}

shelter_green_preds <- lapply(shelter_pred_list, 
                              function(x){
                                lapply(x, 
                                       function(z){
                                         apply(z, 2, function(y) get_ordinal_probs(y, draws$sigma)[[1]])
                                       })
                              })

shelter_amber_preds <- lapply(shelter_pred_list, 
                              function(x){
                                lapply(x, 
                                       function(z){
                                         apply(z, 2, function(y) get_ordinal_probs(y, draws$sigma)[[2]])
                                       })
                              })

shelter_red_preds <- lapply(shelter_pred_list, 
                              function(x){
                                lapply(x, 
                                       function(z){
                                         apply(z, 2, function(y) get_ordinal_probs(y, draws$sigma)[[3]])
                                       })
                              })

shelter_missing_preds <- lapply(shelter_missing_pred_list, 
                            function(x){
                              lapply(x, 
                                     function(z){
                                       apply(z, 2, inv_logit)
                                     })
                            })
set.seed(2020)
plot_dogs <- 241
n_each_samples <- 1

# green codes
png(filename = paste0(figure_path,"/figure_1_a.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, 1:10, 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after arrival", ylab = "probability",
     cex.lab = 1.5)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("a", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), n_each_samples, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_grid_Z, 
          unlist(lapply(shelter_green_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("darkolivegreen3", 0.6))
  }
}


polygon(
  x = c(days_grid_Z, rev(days_grid_Z)),
  y = c(
    unlist(lapply( shelter_green_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( shelter_green_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("green",0.5), border = scales::alpha("green",0.5)
)
lines(days_grid_Z, 
      unlist(lapply( shelter_green_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# amber codes
png(filename = paste0(figure_path,"/figure_1_b.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, 1:10, 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after arrival", ylab = "probability",
     cex.lab = 1.5)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("b", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), 1, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_grid_Z, 
          unlist(lapply(shelter_amber_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("orange", 0.6))
  }
}


polygon(
  x = c(days_grid_Z, rev(days_grid_Z)),
  y = c(
    unlist(lapply( shelter_amber_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( shelter_amber_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("darkorange2",0.5), border = scales::alpha("darkorange2",0.5)
)
lines(days_grid_Z, 
      unlist(lapply( shelter_amber_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# red codes
png(filename = paste0(figure_path,"/figure_1_c.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, 1:10, 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after arrival", ylab = "probability",
     cex.lab = 1.5, adj=0.5)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("c", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), 1, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_grid_Z, 
          unlist(lapply(shelter_red_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("tomato", 0.6))
  }
}


polygon(
  x = c(days_grid_Z, rev(days_grid_Z)),
  y = c(
    unlist(lapply( shelter_red_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( shelter_red_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("firebrick4",0.5), border = scales::alpha("firebrick4",0.5)
)
lines(days_grid_Z, 
      unlist(lapply( shelter_red_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# missing codes
png(filename = paste0(figure_path,"/figure_1_d.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, 1:10, 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after arrival", ylab = "probability",
     cex.lab = 1.5, adj=0.5, mar=c(5,4,4,2)+0.1)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("d", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), 1, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_grid_Z, 
          unlist(lapply(shelter_missing_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("steelblue", 0.4))
  }
}

polygon(
  x = c(days_grid_Z, rev(days_grid_Z)),
  y = c(
    unlist(lapply( shelter_missing_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( shelter_missing_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("steelblue",0.8), border = scales::alpha("steelblue",0.8)
)
lines(days_grid_Z, 
      unlist(lapply( shelter_missing_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# FIGURE 2: variation between the contexts

# make clean labels
clean_context_labels <- c("HND", "KNL", "OKNL", "FPL", "UFPL", "FOOD", "TOYS", "DOGS")
context_labels_numbers <- c(4, 5, 6, 3, 2, 8, 7, 1)
context_intercepts_ordered <- random_effect_list$r_s_g_intercept[, context_labels_numbers]
context_slopes_ordered <- random_effect_list$r_s_g_slope[, context_labels_numbers]
context_intercepts_missing_ordered <- random_effect_list$r_s_g_missing[,context_labels_numbers]
nudge <- 0.15

# Figure S1a
png(paste0(figure_path, "/figure_S1_a.png"), width=1100, height=1000, res=200)
plot(1:8, 
     apply(context_intercepts_ordered, 2, mean), 
     type="n",
     pch=16,
     bty="n", axes = F,
     ylim = c(-3,2), 
     xlab = "contexts", ylab="latent behaviour scale",
     cex.lab = 1.5
     )
mtext("a", font=2, at=1, line=0, cex=1.5)
axis(side = 1, at = 1:8, labels=FALSE, lwd=3)# labels = clean_context_labels, padj = 0, lwd = 3, cex.axis = 1)
text(1:8, par("usr")[3] - 0.4, labels = clean_context_labels, srt = 45, pos = 1, xpd = TRUE)
axis(side = 2,lwd = 3, cex.axis = 1.25)
abline(h=mean(draws$alpha),lty=2, lwd=2)
abline(h=mean(draws$beta),lty=2, lwd=2, col="darkgray")
abline(h=0, lwd=2, lty=3, col="red")
segments(x0 = 1:8-nudge, x1 = 1:8-nudge, 
         y0 = apply(draws$alpha + context_intercepts_ordered, 2, hdi)[1,],
         y1 = apply(draws$alpha + context_intercepts_ordered, 2, hdi)[2,], 
         lwd=3
)
segments(x0 = 1:8+nudge, x1 = 1:8+nudge, 
         y0 = apply(draws$beta + context_slopes_ordered, 2, hdi)[1,],
         y1 = apply(draws$beta + context_slopes_ordered, 2, hdi)[2,], 
         lwd=3, col="darkgray"
)
points(1:8-nudge, apply(draws$alpha + context_intercepts_ordered, 2, mean), pch=16, cex=1.3)
points(1:8+nudge, apply(draws$beta + context_slopes_ordered, 2, mean), pch=16, col="darkgray", cex=1.3)
legend("topright", 
       inset = c(-0.0, -0.03),
       bty="n",
       legend = c("intercept", "slope"), 
       pch = c(16, 16), 
       col = c("black", "darkgray"), 
       lwd=3, 
       cex = 1.1)
dev.off()

# figure S1b
png(paste0(figure_path, "/figure_S1_b.png"), width=1100, height=1000, res=200)
plot(1:8, 
     apply(context_intercepts_missing_ordered, 2, mean), 
     type="n",
     pch=16,
     bty="n", axes = F,
     ylim = c(0, 1), 
     xlab = "contexts", ylab="probability",
     cex.lab = 1.5
)
mtext("b", font=2, at=1, line=0, cex=1.5)
axis(side = 1, at = 1:8, labels=FALSE, lwd=3)
text(1:8, par("usr")[3] - 0.075, labels = clean_context_labels, srt = 45, pos = 1, xpd = TRUE)
axis(side = 2, lwd = 3, cex.axis = 1.25)
abline(h=inv_logit(mean(draws$alpha_miss)), lty=2, lwd=2, col="steelblue")
segments(x0 = 1:8, x1 = 1:8, 
         y0 = inv_logit(apply(draws$alpha_miss + context_intercepts_missing_ordered, 2, hdi)[1,]),
         y1 = inv_logit(apply(draws$alpha_miss + context_intercepts_missing_ordered, 2, hdi)[2,]), 
         lwd=3, col="steelblue"
)
points(1:8, inv_logit(apply(draws$alpha_miss + context_intercepts_missing_ordered, 2, mean)), 
       col="steelblue",
       pch=16, cex=1.3)
dev.off()

# Figure 2a
# png(paste0(figure_path, "/figure_2_a.png"), width=1000, height=1000, res=200)
# plot(days_grid_Z, rnorm(length(days_grid_Z), 0, 1), 
#      type="n", ylim=c(-5, 5), 
#      bty="n", axes = F, 
#      xlab="days after arrival", ylab="latent behaviour scale", 
#      cex.lab=1.5)
# axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
#      labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, lwd=3,  cex.axis=1.25)
# mtext("a", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)
# for(i in seq_along(clean_context_labels)){
#   lines(days_grid_Z, 
#         apply(
#           sapply(days_grid_Z, 
#                  function(x) {
#                    draws$alpha + random_effect_list$r_s_g_intercept[,context_labels_numbers[i]] + 
#                      (draws$beta + random_effect_list$r_s_g_slope[,context_labels_numbers[i]]) * x
#                  }
#           ), 
#           2, mean), 
#         lwd=3, 
#         col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
#         )
# }
# legend("topright", 
#        inset = c(0.3, 0),
#        bty="n",
#        legend = clean_context_labels[1:4], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
#        )
# legend("topright", 
#        inset = c(0, 0),
#        bty="n",
#        legend = clean_context_labels[5:8], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
# )
# dev.off()
# 
# # Figure 2b
# png(paste0(figure_path, "/figure_2_b.png"), width=1000, height=1000, res=200)
# plot(days_grid_Z, rnorm(length(days_grid_Z), 0, 1), 
#      type="n", ylim=c(0, 1.2), 
#      bty="n", axes = F, 
#      xlab="days after arrival", ylab="probability missing", 
#      cex.lab=1.5)
# axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
#      labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, at = c(0, 0.5, 1), lwd=3,  cex.axis=1.25)
# mtext("b", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)
# 
# for(i in seq_along(clean_context_labels)){
#   lines(days_grid_Z, 
#         apply(
#           sapply(days_grid_Z, 
#                  function(x) {
#                    inv_logit(draws$alpha_miss + random_effect_list$r_s_g_missing[,context_labels_numbers[i]] + 
#                      draws$beta_miss * x
#                    )
#                  }
#           ), 
#           2, mean), 
#         lwd=3, 
#         col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
#   )
# }
# legend("topright", 
#        inset = c(0.3, -0),
#        bty="n",
#        legend = clean_context_labels[1:4], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
# )
# legend("topright", 
#        inset = c(0, -0),
#        bty="n",
#        legend = clean_context_labels[5:8], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
# )
# dev.off()

#---------------------------------------------------------------------------------------------
# Key post-adoption results

# build data frame of primary shelter results
adoption_main_results_df <- data.frame(
  parameter = c("delta", "gamma", "delta_miss", "gamma_miss", "epsilon", "kappa[2]")
)

adoption_main_results_df$posterior_mean <- apply(draws[,adoption_main_results_df$parameter], 2, mean)
adoption_main_results_df$hdi_2.5 <- apply(draws[,adoption_main_results_df$parameter], 2, hdi)[1,]
adoption_main_results_df$hdi_97.5 <- apply(draws[,adoption_main_results_df$parameter], 2, hdi)[2,]
print(adoption_main_results_df)

# build data frame of dog-level predictor variables
get_gammas <- draws[, grep("Gamma", colnames(draws))]
colnames(get_gammas) <- c(colnames(X_a), paste0(colnames(X_a), ".miss"))

adoption_dl_predictors_df <- data.frame(
  parameter = c("length_of_stay", "weight", "age", "sex_male.female", 
                "neutered_no.yes", "neutered_no.onsite", "neutered_onsite.yes", 
                "source_gift.stray", "source_gift.return", "source_return.stray",
                "number_of_surveys_two.one",
                paste0(
                  c("length_of_stay", "weight", "age", "sex_male.female", 
                    "neutered_no.yes", "neutered_no.onsite", "neutered_onsite.yes", 
                    "source_gift.stray", "source_gift.return", "source_return.stray", 
                    "number_of_surveys_two.one"
                  ), ".miss")), 
  
  posterior_mean = c(
    mean(get_gammas$length_of_stay_Z), 
    mean(get_gammas$latest_weight_Z), 
    mean(get_gammas$estimated_age_at_departure_Z), 
    mean(get_gammas$sex_MALE - get_gammas$sex_MALE*-1),
    mean(get_gammas$neutered_NO - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1),
    mean(get_gammas$neutered_NO - get_gammas$neutered_ONSITE),
    mean(get_gammas$neutered_ONSITE - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1),
    mean(get_gammas$source_GIFT - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1), 
    mean(get_gammas$source_GIFT - get_gammas$source_RETURN), 
    mean(get_gammas$source_RETURN - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1),
    mean(get_gammas$total_number_of_surveys_TWO - get_gammas$total_number_of_surveys_TWO*-1),
    # missing
    mean(get_gammas$length_of_stay_Z.miss), 
    mean(get_gammas$latest_weight_Z.miss), 
    mean(get_gammas$estimated_age_at_departure_Z.miss), 
    mean(get_gammas$sex_MALE.miss - get_gammas$sex_MALE.miss*-1),
    mean(get_gammas$neutered_NO.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1),
    mean(get_gammas$neutered_NO.miss - get_gammas$neutered_ONSITE.miss),
    mean(get_gammas$neutered_ONSITE.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1),
    mean(get_gammas$source_GIFT.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1), 
    mean(get_gammas$source_GIFT.miss - get_gammas$source_RETURN.miss), 
    mean(get_gammas$source_RETURN.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1),
    mean(get_gammas$total_number_of_surveys_TWO.miss - get_gammas$total_number_of_surveys_TWO.miss*-1)
  ),
  
  hdi_2.5 = c(
    hdi(get_gammas$length_of_stay_Z)[1], 
    hdi(get_gammas$latest_weight_Z)[1], 
    hdi(get_gammas$estimated_age_at_departure_Z)[1], 
    hdi(get_gammas$sex_MALE - get_gammas$sex_MALE*-1)[1],
    hdi(get_gammas$neutered_NO - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1)[1],
    hdi(get_gammas$neutered_NO - get_gammas$neutered_ONSITE)[1],
    hdi(get_gammas$neutered_ONSITE - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1)[1],
    hdi(get_gammas$source_GIFT - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1)[1], 
    hdi(get_gammas$source_GIFT - get_gammas$source_RETURN)[1], 
    hdi(get_gammas$source_RETURN - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1)[1],
    hdi(get_gammas$total_number_of_surveys_TWO - get_gammas$total_number_of_surveys_TWO*-1)[1],
    # missing
    hdi(get_gammas$length_of_stay_Z.miss)[1], 
    hdi(get_gammas$latest_weight_Z.miss)[1], 
    hdi(get_gammas$estimated_age_at_departure_Z.miss)[1], 
    hdi(get_gammas$sex_MALE.miss - get_gammas$sex_MALE.miss*-1)[1],
    hdi(get_gammas$neutered_NO.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1)[1],
    hdi(get_gammas$neutered_NO.miss - get_gammas$neutered_ONSITE.miss)[1],
    hdi(get_gammas$neutered_ONSITE.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1)[1],
    hdi(get_gammas$source_GIFT.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1)[1], 
    hdi(get_gammas$source_GIFT.miss - get_gammas$source_RETURN.miss)[1], 
    hdi(get_gammas$source_RETURN.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1)[1],
    hdi(get_gammas$total_number_of_surveys_TWO.miss - get_gammas$total_number_of_surveys_TWO.miss*-1)[1]
  ),
  
  hdi_97.5 = c(
    hdi(get_gammas$length_of_stay_Z)[2], 
    hdi(get_gammas$latest_weight_Z)[2], 
    hdi(get_gammas$estimated_age_at_departure_Z)[2], 
    hdi(get_gammas$sex_MALE - get_gammas$sex_MALE*-1)[2],
    hdi(get_gammas$neutered_NO - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1)[2],
    hdi(get_gammas$neutered_NO - get_gammas$neutered_ONSITE)[2],
    hdi(get_gammas$neutered_ONSITE - (get_gammas$neutered_NO + get_gammas$neutered_ONSITE)*-1)[2],
    hdi(get_gammas$source_GIFT - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1)[2], 
    hdi(get_gammas$source_GIFT - get_gammas$source_RETURN)[2], 
    hdi(get_gammas$source_RETURN - (get_gammas$source_GIFT + get_gammas$source_RETURN)*-1)[2],
    hdi(get_gammas$total_number_of_surveys_TWO - get_gammas$total_number_of_surveys_TWO*-1)[2],
    # missing
    hdi(get_gammas$length_of_stay_Z.miss)[2], 
    hdi(get_gammas$latest_weight_Z.miss)[2], 
    hdi(get_gammas$estimated_age_at_departure_Z.miss)[2], 
    hdi(get_gammas$sex_MALE.miss - get_gammas$sex_MALE.miss*-1)[2],
    hdi(get_gammas$neutered_NO.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1)[2],
    hdi(get_gammas$neutered_NO.miss - get_gammas$neutered_ONSITE.miss)[2],
    hdi(get_gammas$neutered_ONSITE.miss - (get_gammas$neutered_NO.miss + get_gammas$neutered_ONSITE.miss)*-1)[2],
    hdi(get_gammas$source_GIFT.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1)[2], 
    hdi(get_gammas$source_GIFT.miss - get_gammas$source_RETURN.miss)[2], 
    hdi(get_gammas$source_RETURN.miss - (get_gammas$source_GIFT.miss + get_gammas$source_RETURN.miss)*-1)[2],
    hdi(get_gammas$total_number_of_surveys_TWO.miss - get_gammas$total_number_of_surveys_TWO.miss*-1)[2]
  )
)

print(adoption_dl_predictors_df)

#---------------------------------------------------------------------------------------------
# FIGURE 2: plot code colour probabilities across days after adoption

# days to plot
min_max_a_days <- c(0, mean(d_a_stan$days_after_adoption) + sd(dog_means_days_after_adoption))
# grid of days
days_a_grid <- seq(min_max_a_days[1], min_max_a_days[2], length.out = 10)
# standardised grid of days (use values to standardise from analysis_script)
days_a_grid_Z <- ( days_a_grid - mean(dog_means_days_after_adoption))/ sd(dog_means_days_after_adoption)

adoption_pred_list <- rep(list(list(list())), length(days_a_grid_Z))
adoption_missing_pred_list <- rep(list(list(list())), length(days_a_grid_Z))

for(i in 1:length(days_a_grid_Z)){#for each day
  message(paste0("Day ", i))
  for(j in 1:length(unique(d_a_stan$dog_id))){#for each dog
    
    # container
    out_jg <- rep(list(list()), length(unique(d_a_stan$context_id)))
    out_miss_jg <- rep(list(list()), length(unique(d_a_stan$context_id)))
    
    for(g in 1:length(unique(d_a_stan$context_id))){#for each context
      
      out_jg[[g]] <- draws$delta + # mean
        random_effect_list$r_a_j_intercept[,j] + # dog
        random_effect_list$r_a_g_intercept[,g] + # context
        random_effect_list$r_a_jg_intercept[,make_cc_interaction(j,g)] + # dog x context
        (draws$gamma + 
           random_effect_list$r_a_j_slope[,j] + # dog
           random_effect_list$r_a_g_slope[,g] + # context
           random_effect_list$r_a_jg_slope[,make_cc_interaction(j,g)]) * days_a_grid_Z[i]
      
      out_miss_jg[[g]] <- draws$delta_miss + # mean
        random_effect_list$r_a_j_missing[,j] + # dog
        random_effect_list$r_a_g_missing[,g] + # context
        random_effect_list$r_a_jg_missing[,make_cc_interaction(j,g)] + # dog x context
        draws$gamma_miss * days_a_grid_Z[i]
      
    }#g
    
    adoption_pred_list[[i]][[j]] <- matrix(unlist(out_jg), ncol=8)
    adoption_missing_pred_list[[i]][[j]] <- matrix(unlist(out_miss_jg), ncol=8)
    
  }
}

adoption_green_preds <- lapply(adoption_pred_list, 
                              function(x){
                                lapply(x, 
                                       function(z){
                                         apply(z, 2, function(y) get_ordinal_probs(y, draws$epsilon)[[1]])
                                       })
                              })

adoption_amber_preds <- lapply(adoption_pred_list, 
                              function(x){
                                lapply(x, 
                                       function(z){
                                         apply(z, 2, function(y) get_ordinal_probs(y, draws$epsilon)[[2]])
                                       })
                              })

adoption_red_preds <- lapply(adoption_pred_list, 
                            function(x){
                              lapply(x, 
                                     function(z){
                                       apply(z, 2, function(y) get_ordinal_probs(y, draws$epsilon)[[3]])
                                     })
                            })

adoption_missing_preds <- lapply(adoption_missing_pred_list, 
                                function(x){
                                  lapply(x, 
                                         function(z){
                                           apply(z, 2, inv_logit)
                                         })
                                })

set.seed(2020)
plot_dogs <- 241
n_each_samples <- 1

# Figure 2a
png(filename = paste0(figure_path,"/figure_2_a.png"), width=1000, height=1000, res=200)
plot(days_a_grid_Z, 1:length(days_a_grid_Z), 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after adoption", ylab = "probability",
     cex.lab = 1.5)
axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
     labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("a", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), n_each_samples, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_a_grid_Z, 
          unlist(lapply(adoption_green_preds, 
                        function(x) x[[i]][row_samples[j], context_sample])),
          lwd=0.5, col=scales::alpha("darkolivegreen3", 0.6))
  }
}


polygon(
  x = c(days_a_grid_Z, rev(days_a_grid_Z)),
  y = c(
    unlist(lapply( adoption_green_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( adoption_green_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("green",0.5), border = scales::alpha("green",0.5)
)
lines(days_a_grid_Z, 
      unlist(lapply( adoption_green_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# amber codes
png(filename = paste0(figure_path,"/figure_2_b.png"), width=1000, height=1000, res=200)
plot(days_a_grid_Z, 1:length(days_a_grid_Z), 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after adoption", ylab = "probability",
     cex.lab = 1.5)
axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
     labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("b", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), n_each_samples, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_a_grid_Z, 
          unlist(lapply(adoption_amber_preds, 
                        function(x) x[[i]][row_samples[j], context_sample])),
          lwd=0.5, col=scales::alpha("orange", 0.6))
  }
}


polygon(
  x = c(days_a_grid_Z, rev(days_a_grid_Z)),
  y = c(
    unlist(lapply( adoption_amber_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( adoption_amber_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("darkorange2",0.5), border = scales::alpha("darkorange2",0.5)
)
lines(days_a_grid_Z, 
      unlist(lapply( adoption_amber_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# red codes
png(filename = paste0(figure_path,"/figure_2_c.png"), width=1000, height=1000, res=200)
plot(days_a_grid_Z, 1:length(days_a_grid_Z), 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after adoption", ylab = "probability",
     cex.lab = 1.5, adj=0.5)
axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
     labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("c", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), n_each_samples, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_a_grid_Z, 
          unlist(lapply(adoption_red_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("tomato", 0.6))
  }
}


polygon(
  x = c(days_a_grid_Z, rev(days_a_grid_Z)),
  y = c(
    unlist(lapply( adoption_red_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( adoption_red_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("firebrick4",0.5), border = scales::alpha("firebrick4",0.5)
)
lines(days_a_grid_Z, 
      unlist(lapply( adoption_red_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# missing codes
png(filename = paste0(figure_path,"/figure_2_d.png"), width=1000, height=1000, res=200)
plot(days_a_grid_Z, 1:length(days_a_grid_Z), 
     type="n", ylim=c(0,1), bty="n", axes = F,
     xlab = "days after adoption", ylab = "probability",
     cex.lab = 1.5, adj=0.5, mar=c(5,4,4,2)+0.1)
axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
     labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = seq(0, 1, 0.5), lwd=3,  cex.axis=1.25)
mtext("d", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)

for(i in 1:plot_dogs){
  row_samples <- sample(1:nrow(draws), n_each_samples, replace = FALSE)
  context_sample <- sample(1:8, 1)
  for(j in seq_along(row_samples)){
    lines(days_a_grid_Z, 
          unlist(lapply(adoption_missing_preds, 
                        function(x) x[[i]][row_samples, context_sample])),
          lwd=0.5, col=scales::alpha("steelblue", 0.4))
  }
}

polygon(
  x = c(days_a_grid_Z, rev(days_a_grid_Z)),
  y = c(
    unlist(lapply( adoption_missing_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[1])),
    rev(unlist(lapply( adoption_missing_preds, function(x) hdi(unlist(lapply(x, function(z) mean(colMeans(z)))))[2])))
  ),
  col = scales::alpha("steelblue",0.8), border = scales::alpha("steelblue",0.8)
)
lines(days_a_grid_Z, 
      unlist(lapply(adoption_missing_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()


#---------------------------------------------------------------------------------------------------------------
# FIGURE S2 & 4: variation between the contexts post-adoption

# make clean labels
clean_context_labels <- c("HND", "HOUSE", "OTSD", "FPL", "UFPL", "FOOD", "TOYS", "DOGS")
context_labels_numbers <- c(4, 5, 6, 3, 2, 8, 7, 1)
context_intercepts_ordered <- random_effect_list$r_a_g_intercept[, context_labels_numbers]
context_slopes_ordered <- random_effect_list$r_a_g_slope[, context_labels_numbers]
context_intercepts_missing_ordered <- random_effect_list$r_a_g_missing[,context_labels_numbers]
nudge <- 0.15

# Figure S2
png(paste0(figure_path, "/figure_S2_a.png"), width=1100, height=1000, res=200)
plot(1:8, 
     apply(context_intercepts_ordered, 2, mean), 
     type="n",
     pch=16,
     bty="n", axes = F,
     ylim = c(-2,3), 
     xlab = "contexts", ylab="latent behaviour scale",
     cex.lab = 1.5
)
mtext("a", font=2, at=1, line=0, cex=1.5)
axis(side = 1, at = 1:8, labels=FALSE, lwd=3)# labels = clean_context_labels, padj = 0, lwd = 3, cex.axis = 1)
text(1:8, par("usr")[3] - 0.4, labels = clean_context_labels, srt = 45, pos = 1, xpd = TRUE)
axis(side = 2,lwd = 3, cex.axis = 1.25)
abline(h=mean(draws$delta),lty=2, lwd=2)
abline(h=mean(draws$gamma),lty=2, lwd=2, col="darkgray")
abline(h=0, lwd=2, lty=3, col="red")
segments(x0 = 1:8-nudge, x1 = 1:8-nudge, 
         y0 = apply(draws$delta + context_intercepts_ordered, 2, hdi)[1,],
         y1 = apply(draws$delta + context_intercepts_ordered, 2, hdi)[2,], 
         lwd=3
)
segments(x0 = 1:8+nudge, x1 = 1:8+nudge, 
         y0 = apply(draws$gamma + context_slopes_ordered, 2, hdi)[1,],
         y1 = apply(draws$gamma + context_slopes_ordered, 2, hdi)[2,], 
         lwd=3, col="darkgray"
)
points(1:8-nudge, apply(draws$delta + context_intercepts_ordered, 2, mean), pch=16, cex=1.3)
points(1:8+nudge, apply(draws$gamma + context_slopes_ordered, 2, mean), pch=16, col="darkgray", cex=1.3)
legend("topright", 
       inset = c(-0.0, -0.03),
       bty="n",
       legend = c("intercept", "slope"), 
       pch = c(16, 16), 
       col = c("black", "darkgray"), 
       lwd=3, 
       cex = 1.1)
dev.off()

# figure S2b
png(paste0(figure_path, "/figure_S2_b.png"), width=1100, height=1000, res=200)
plot(1:8, 
     apply(context_intercepts_missing_ordered, 2, mean), 
     type="n",
     pch=16,
     bty="n", axes = F,
     ylim = c(0, 0.1), 
     xlab = "contexts", ylab="probability",
     cex.lab = 1.5
)
mtext("b", font=2, at=1, line=0, cex=1.5)
axis(side = 1, at = 1:8, labels=FALSE, lwd=3)
text(1:8, par("usr")[3] - 0.01, labels = clean_context_labels, srt = 45, pos = 1, xpd = TRUE)
axis(side = 2, lwd = 3, cex.axis = 1.25)
abline(h=inv_logit(mean(draws$delta_miss)), lty=2, lwd=2, col="steelblue")
segments(x0 = 1:8, x1 = 1:8, 
         y0 = inv_logit(apply(draws$delta_miss + context_intercepts_missing_ordered, 2, hdi)[1,]),
         y1 = inv_logit(apply(draws$delta_miss + context_intercepts_missing_ordered, 2, hdi)[2,]), 
         lwd=3, col="steelblue"
)
points(1:8, inv_logit(apply(draws$delta_miss + context_intercepts_missing_ordered, 2, mean)), 
       col="steelblue",
       pch=16, cex=1.3)
dev.off()

# # Figure 2a OLD
# png(paste0(figure_path, "/figure_2_a.png"), width=1000, height=1000, res=200)
# plot(days_grid_Z, rnorm(length(days_grid_Z), 0, 1), 
#      type="n", ylim=c(-5, 5), 
#      bty="n", axes = F, 
#      xlab="days after arrival", ylab="latent behaviour scale", 
#      cex.lab=1.5)
# axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
#      labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, lwd=3,  cex.axis=1.25)
# mtext("a", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)
# for(i in seq_along(clean_context_labels)){
#   lines(days_grid_Z, 
#         apply(
#           sapply(days_grid_Z, 
#                  function(x) {
#                    draws$alpha + random_effect_list$r_s_g_intercept[,context_labels_numbers[i]] + 
#                      (draws$beta + random_effect_list$r_s_g_slope[,context_labels_numbers[i]]) * x
#                  }
#           ), 
#           2, mean), 
#         lwd=3, 
#         col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
#   )
# }
# legend("topright", 
#        inset = c(0.3, 0),
#        bty="n",
#        legend = clean_context_labels[1:4], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
# )
# legend("topright", 
#        inset = c(0, 0),
#        bty="n",
#        legend = clean_context_labels[5:8], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
# )
# dev.off()
# 
# # Figure 4a
# png(paste0(figure_path, "/figure_4_a.png"), width=1000, height=1000, res=200)
# plot(days_a_grid_Z, rnorm(length(days_a_grid_Z), 0, 1), 
#      type="n", ylim=c(-5, 5), 
#      bty="n", axes = F, 
#      xlab="days after adoption", ylab="latent behaviour scale", 
#      cex.lab=1.5)
# axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
#      labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, lwd=3,  cex.axis=1.25)
# mtext("a", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)
# for(i in seq_along(clean_context_labels)){
#   lines(days_a_grid_Z, 
#         apply(
#           sapply(days_a_grid_Z, 
#                  function(x) {
#                    draws$delta + random_effect_list$r_a_g_intercept[,context_labels_numbers[i]] + 
#                      (draws$gamma + random_effect_list$r_a_g_slope[,context_labels_numbers[i]]) * x
#                  }
#           ), 
#           2, mean), 
#         lwd=3, 
#         col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
#   )
# }
# legend("topright", 
#        inset = c(0.3, 0),
#        bty="n",
#        legend = clean_context_labels[1:4], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
# )
# legend("topright", 
#        inset = c(0, 0),
#        bty="n",
#        legend = clean_context_labels[5:8], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
# )
# dev.off()

# # Figure 4b
# png(paste0(figure_path, "/figure_4_b.png"), width=1000, height=1000, res=200)
# plot(days_a_grid_Z, rnorm(length(days_a_grid_Z), 0, 1), 
#      type="n", ylim=c(0, 0.3), 
#      bty="n", axes = F, 
#      xlab="days after arrival", ylab="probability missing", 
#      cex.lab=1.5)
# axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
#      labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, at = c(0, 0.5, 1), lwd=3,  cex.axis=1.25)
# mtext("b", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)
# 
# for(i in seq_along(clean_context_labels)){
#   lines(days_a_grid_Z, 
#         apply(
#           sapply(days_a_grid_Z, 
#                  function(x) {
#                    inv_logit(draws$delta_miss + random_effect_list$r_a_g_missing[,context_labels_numbers[i]] + 
#                                draws$gamma_miss * x
#                    )
#                  }
#           ), 
#           2, mean), 
#         lwd=3, 
#         col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
#   )
# }
# legend("topright", 
#        inset = c(0.3, -0),
#        bty="n",
#        legend = clean_context_labels[1:4], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
# )
# legend("topright", 
#        inset = c(0, -0),
#        bty="n",
#        legend = clean_context_labels[5:8], 
#        lwd = 3, 
#        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
# )
# dev.off()

#-----------------------------------------------------------------------------------
# Figure 3
colours_ <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")

png(paste0(figure_path, "/figure_3_a.png"), width=1000, height=1000, res=200)
plot(density(
  draws$`sigma_j[1]`,
  adjust = 0.1
  ), 
  type="n",
  bty="n", axes=F,
  xlab = "repeatability", 
  ylab = "density",
  main = "",
  xlim=c(0,1), ylim=c(0, 25),
  cex.lab=1.5)
axis(side = 1, lwd = 3, cex.axis=1.25)
axis(side = 2, lwd=3,  cex.axis=1.25)
mtext("a", line=0, at=0, cex = 1.5, font = 2)
shelter_denom <- draws$`sigma_j[1]`^2 + draws$`sigma_g[1]`^2 + draws$`sigma_jg[1]`^2 + draws$sigma^2
adoption_denom <- draws$`sigma_j[4]`^2 + draws$`sigma_g[4]`^2 + draws$`sigma_jg[4]`^2 + draws$epsilon^2
lines(density(draws$`sigma_j[1]`^2/shelter_denom, adjust = 0.5), col=colours_[1], lwd=3)
lines(density(draws$`sigma_j[4]`^2/adoption_denom, adjust = 0.5), col=colours_[1], lwd=3, lty=2)
lines(density(draws$`sigma_g[1]`^2/shelter_denom, adjust = 0.5), col=colours_[2], lwd=3)
lines(density(draws$`sigma_g[4]`^2/adoption_denom, adjust = 0.5), col=colours_[2], lwd=3, lty=2)
lines(density(draws$`sigma_jg[1]`^2/shelter_denom, adjust = 0.5), col=colours_[3], lwd=3)
lines(density(draws$`sigma_jg[4]`^2/adoption_denom, adjust = 0.5), col=colours_[3], lwd=3, lty=2)

legend(
  "topright", 
  inset=c(-0.02, -0.05),
  bty="n",
  legend = c("dogs", "contexts", "dogs x contexts"),
  lwd=c(2,2,2),
  seg.len = 1.5,
  col = colours_,
  cex = 1.15
)
legend(
  "topright",
  inset=c(0.02, 0.2),
  bty="n",
  legend = c("shelter", "post-adoption", ""),
  lwd=c(2,2,0), lty=c(1,2,0),
  seg.len = 1.5,
  cex = 1.15
)
#text(x = 0.7, y = 5, "solid line = shelter", cex = 0.8)
dev.off()

#--------------------------------------------------------------------------------------
# correlations between the random effects

Rho_j_raw <- draws[ ,grep("\\bRho_j\\b", colnames(draws))]
Rho_g_raw <- draws[ ,grep("\\bRho_g\\b", colnames(draws))]
Rho_jg_raw <- draws[ ,grep("\\bRho_jg\\b", colnames(draws))]

Rho_res <- data.frame(
  parameter_name = c(colnames(Rho_j_raw), colnames(Rho_g_raw), colnames(Rho_jg_raw)),
  re_type = rep(c("dogs","contexts","dogs_x_contexts"), each=36),
  mu = apply(cbind(Rho_j_raw, Rho_g_raw, Rho_jg_raw), 2, mean),
  hid_2.5 = apply(cbind(Rho_j_raw, Rho_g_raw, Rho_jg_raw), 2, hdi)[1,],
  hid_97.5 = apply(cbind(Rho_j_raw, Rho_g_raw, Rho_jg_raw), 2, hdi)[2,], 
  row.names = NULL
)

# which have HDIs surpassing zero
Rho_res[ which(Rho_res$hid_2.5 > 0 | Rho_res$hid_97.5 < 0), ]

#------------------------------------------------------------------------------------
# FIGURE 3b: correlations

# 3b (personality)
png(paste0(figure_path, "/figure_3_b.png"), width=1000, height=1000, res=200)
plot(density(
  Rho_jg_raw$`Rho_jg[4,1]`, 
  adjust = 0.1
  ), 
  type="n",
  bty="n", axes=F,
  xlab = "correlation", 
  ylab = "density",
  main = "",
  xlim=c(-1,1), ylim=c(0, 6),
  cex.lab=1.5)
axis(side = 1, at = c(-1, 0, 1), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, lwd=3,  cex.axis=1.25)
mtext("b", line=0, at=-1, cex = 1.5, font = 2)
abline(v=0, lty=2, col=scales::alpha("black", 0.5))
lines(density(Rho_j_raw$`Rho_j[4,1]`, adjust = 0.5), col=colours_[1], lwd=3)
lines(density(Rho_g_raw$`Rho_g[4,1]`, adjust = 0.5), col=colours_[2], lwd=3)
lines(density(Rho_jg_raw$`Rho_jg[4,1]`, adjust = 0.5), col=colours_[3], lwd=3)
legend(
  "topright", 
  inset=c(-0.02, -0.05),
  bty="n",
  legend = c("dogs", "contexts", "dogs x contexts"),
  lwd=c(3,3,3),
  seg.len = 1,
  col = colours_,
  cex = 1.15
)
dev.off()

# Figure 5b
# png(paste0(figure_path, "/figure_5_b.png"), width=1000, height=1000, res=200)
# plot(density(
#   Rho_j_raw$`Rho_j[2,5]`, 
#   adjust = 0.1
#   ), 
#   type="n",
#   bty="n", axes=F,
#   xlab = "correlation", 
#   ylab = "density",
#   main = "",
#   xlim=c(-1,1), ylim=c(0, 6),
#   cex.lab=1.5)
# axis(side = 1, at = c(-1, 0, 1), 
#      lwd = 3, cex.axis=1.25)
# axis(side = 2, at=c(0:6), lwd=3,  cex.axis=1.25)
# mtext("b", line=0, at=-1, cex = 1.5, font = 2)
# colours_ <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")
# abline(v=0, lty=2, col=scales::alpha("black", 0.5))
# lines(density(Rho_j_raw$`Rho_j[2,5]`, adjust = 0.5), col=colours_[1], lwd=3)
# lines(density(Rho_g_raw$`Rho_g[2,5]`, adjust = 0.5), col=colours_[2], lwd=3)
# lines(density(Rho_jg_raw$`Rho_jg[2,5]`, adjust = 0.5), col=colours_[3], lwd=3)
# legend(
#   "topright", 
#   inset=c(-0.02, -0.05),
#   bty="n",
#   legend = c("dogs", "contexts", "dogs x contexts"),
#   lwd=c(3,3,3),
#   seg.len = 1,
#   col = colours_,
#   cex = 1.15
# )
# dev.off()

#------------------------------------------------------------------------------------------------------
# VARIATION BETWEEN DOGS
# png("~/Documents/sds1.png", width=1000, height=1000, res=200)
# plot(density(
#   draws$`sigma_j[1]`,
#   adjust = 0.1
#   ), 
#   type="n",
#   bty="n", axes=F,
#   xlab = "standard deviation", 
#   ylab = "density",
#   main = "",
#   xlim=c(0,2), ylim=c(0, 15),
#   cex.lab=1.5)
# axis(side = 1, lwd = 3, cex.axis=1.25)
# axis(side = 2, lwd=3,  cex.axis=1.25)
# lines(density(draws$`sigma_j[1]`, adjust = 0.5), col=colours_[1], lwd=3)
# lines(density(draws$`sigma_j[4]`, adjust = 0.5), col=colours_[1], lwd=3, lty=2)
# lines(density(draws$`sigma_g[1]`, adjust = 0.5), col=colours_[2], lwd=3)
# lines(density(draws$`sigma_g[4]`, adjust = 0.5), col=colours_[2], lwd=3, lty=2)
# lines(density(draws$`sigma_jg[1]`, adjust = 0.5), col=colours_[3], lwd=3)
# lines(density(draws$`sigma_jg[4]`, adjust = 0.5), col=colours_[3], lwd=3, lty=2)
# legend(
#   "topright", 
#   inset=c(-0.02, -0.05),
#   bty="n",
#   legend = c("dogs", "contexts", "dogs x contexts"),
#   lwd=c(3,3,3),
#   seg.len = 1,
#   col = colours_,
#   cex = 1.15
# )
# legend(
#   "topleft", 
#   inset=c(0, -0.05),
#   bty="n",
#   legend = c("shelter", "post-adoption"),
#   lwd=c(3,3), lty=c(1,2),
#   seg.len = 2,
#   cex = 1.15
# )
# #text(x = 1.6, y = 3, "solid line = shelter", cex = 0.8)
# dev.off()


#------------------------------------------------------------------------------------------------------
# ANALYSIS OF SHELTER-ADOPTION BEHAVIOUR DIFFERENCES

old_context_labels <- c("DOGS", "UFPL", "FPL", "HND", "KNL/HOUSE", "OKNL/OTSD", "TOYS", "FOOD")
sorted_context_numbers <- c(8, 5, 4, 1, 2, 3, 7, 6)
clean_context_labels <- c("HND", "KNL/HOUSE", "OKNL/OSTD", "FPL", "UFPL", "FOOD", "TOYS", "DOGS")
context_labels_numbers <- c(4, 5, 6, 3, 2, 8, 7, 1)

# Figure4
n_dogs <- 4
sample_dogs <- c(9, 89, 112, 155) #sample(1:241, n_dogs, replace=F)
which_contexts <- c(8,1,7,6) #sample(1:8, n_dogs, replace = F)
# op <- par(mfrow = c(1,3),
#           oma = c(5,4,0,0) + 0.1,
#           mar = c(0,0,1,1) + 0.1)

for(i in seq_along(sample_dogs)){  
    png( paste0(figure_path, "/figure_4_", letters[i], ".png"), width = 2000, height = 500, res=200)
    par(mfrow=c(1,3))
    label_size <- 1.5
    dogs_shelter_data <- d_s_stan[d_s_stan$dog_id %in% sample_dogs[i] & d_s_stan$context_id %in% which_contexts[i], 
                                  c("day_since_arrival_Z", "day_since_arrival", "behaviour_code_colour")]
    
    dogs_adoption_data <- d_a_stan[d_a_stan$dog_id %in% sample_dogs[i] & d_a_stan$context_id %in% which_contexts[i], 
                                   c("days_after_adoption_Z", "days_after_adoption", "behaviour_code_colour")]
    
    total_length <- nrow(dogs_shelter_data) + nrow(dogs_adoption_data)
    
    dogs_shelter_preds <- sapply(dogs_shelter_data$day_since_arrival_Z, 
                                 function(x){
                                   draws$alpha + 
                                     random_effect_list$r_s_j_intercept[,sample_dogs[i]] + 
                                     random_effect_list$r_s_g_intercept[,which_contexts[i]] + 
                                     random_effect_list$r_s_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[i])] + 
                                     (
                                       draws$beta + 
                                         random_effect_list$r_s_j_slope[,sample_dogs[i]] + 
                                         random_effect_list$r_s_g_slope[,which_contexts[i]] + 
                                         random_effect_list$r_s_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[i])]
                                     ) * x
                                 })
    
    dogs_adoption_preds <- sapply(dogs_adoption_data$days_after_adoption_Z, 
                                  function(x){
                                    draws$delta + 
                                      random_effect_list$r_a_j_intercept[,sample_dogs[i]] + 
                                      random_effect_list$r_a_g_intercept[,which_contexts[i]] + 
                                      random_effect_list$r_a_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[i])] + 
                                      (
                                        draws$gamma + 
                                          random_effect_list$r_a_j_slope[,sample_dogs[i]] + 
                                          random_effect_list$r_a_g_slope[,which_contexts[i]] + 
                                          random_effect_list$r_a_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[i])]
                                      ) * x
                                  })
    
    
    dog_predict_df <- data.frame(
      dog_id = rep(sample_dogs[i], total_length),
      which_context = rep(which_contexts[i], total_length),
      days = c(dogs_shelter_data$day_since_arrival, max(dogs_shelter_data$day_since_arrival) + dogs_adoption_data$days_after_adoption),
      codes = c(ifelse(is.na(dogs_shelter_data$behaviour_code_colour), 0, dogs_shelter_data$behaviour_code_colour), 
                ifelse(is.na(dogs_adoption_data$behaviour_code_colour), 0, dogs_adoption_data$behaviour_code_colour)),
      time = rep(c("shelter","adoption"), c(nrow(dogs_shelter_data), nrow(dogs_adoption_data))),
      mu = c(apply(dogs_shelter_preds, 2, mean), apply(dogs_adoption_preds, 2, mean)),
      hdi_2.5 = c(apply(dogs_shelter_preds, 2, hdi)[1,], apply(dogs_adoption_preds, 2, hdi)[1,]),
      hdi_97.5 = c(apply(dogs_shelter_preds, 2, hdi)[2,], apply(dogs_adoption_preds, 2, hdi)[2,])
    )
    
    dog_predict_df$code_colour <- ifelse(
      dog_predict_df$codes %in% 0, scales::alpha("black", 0.1), 
      ifelse(
        dog_predict_df$codes %in% 1, "darkolivegreen4", 
        ifelse(
          dog_predict_df$codes %in% 2, "darkorange2",
          "firebrick2"
        )
      )
    )
    
    plot(dog_predict_df$days, dog_predict_df$codes, type="n", ylim=c(-5,5), 
         axes=F, bty="n",
         xlab = "observation days", ylab = "behaviour", 
         cex.lab = label_size)
    abline(h=c(1, 2, 3), lty=3, 
           col=scales::alpha("black", 0.5))
    polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5, 
                min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5,
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1, 
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1), 
            y = c(-10, 10, 10, -10), 
            col = scales::alpha("black", 0.1), border = NA
    )
    polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1, 
                min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1,
                max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5, 
                max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5), 
            y = c(-10, 10, 10, -10), 
            col = scales::alpha("black", 0.1), border = NA
    )
    axis(side = 2, cex.axis=1, las=1,
         at = c(-5, 1,2,3,5))
    axis(side=1, cex.axis=1,
         at = c(0, 
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"]),
                min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"]),
                max(dog_predict_df$days))
    )
    
    mtext(letters[i], cex=1, adj=0, font=2)
    mtext(paste0(sample_dogs[i], " ", clean_context_labels[which(context_labels_numbers %in% which_contexts[i])]),
          cex = 1, adj = 0, at = 3)
    
    # shelter preds
    sample_lines <- sample(1:nrow(draws), 50, replace=T)
    for(j in seq_along(sample_lines)){
      lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"],
            dogs_shelter_preds[sample_lines[j], ],
            col = scales::alpha("slateblue", 0.2))
    } 
    for(j in seq_along(sample_lines)){
      lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"],
            dogs_adoption_preds[sample_lines[j], ],
            col = scales::alpha("slateblue", 0.2))
    }
    
    points(dog_predict_df$days, ifelse(is.na(dog_predict_df$codes), 0, dog_predict_df$codes), 
           pch=ifelse(dog_predict_df$codes > 0, 16, 1), 
           cex = ifelse(dog_predict_df$codes > 0, 1, 0),
           col = scales::alpha(dog_predict_df$code_colour, 1))
    
    # shelter code probs plot
    shelter_code_probs <- table(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "shelter", "codes"])/
      length(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "shelter", "codes"])
    
    which_codes <- as.numeric(names(shelter_code_probs))
    
    if(length(which_codes) != 3){
      missing_codes <- numeric()
      for(j in 1:3){
        if(length(which_codes[which_codes %in% j]) == 0) {
          missing_codes <- c(missing_codes, j)
          #which_codes <- sort(c(which_codes, i), decreasing = F)
        }
      }
      
      temp_codes <- rep(NA, 3)
      for(j in 1:3){
        if(j %in% which_codes) temp_codes[j] <- shelter_code_probs[names(shelter_code_probs)==j] else temp_codes[j] <- 0
      }
      
      shelter_code_probs <- temp_codes
      names(shelter_code_probs) <- 1:3
    }
    
    # ordinal code probs
    shelter_pred_codes <- do.call("rbind.data.frame", apply(dogs_shelter_preds, 2, function(x) get_ordinal_probs(x, draws$sigma)))
    
    the_colours <- c("darkolivegreen4","darkorange2","firebrick2")
    
    plot(density(rowMeans(dogs_shelter_preds)), lwd=2, col=scales::alpha("slateblue",0.8), 
         xlim=c(-5,5), ylim=c(0,1.6), bty="n", axes=F, main="", type="n",
         xlab = "behaviour", 
         ylab = "probability & density", 
         cex.lab = label_size)
    for(l in 1:3) if(shelter_code_probs[l] > 0) lines(c(l,l), c(0, shelter_code_probs[l]), lwd=10, col=scales::alpha(the_colours[l],0.6) )
    points(1:3, colMeans(shelter_pred_codes))
    segments(x0 = 1:3, x1 = 1:3, 
             y0 = apply(shelter_pred_codes, 2, hdi)[1,], 
             y1 = apply(shelter_pred_codes, 2, hdi)[2,]
    )
    lines(density(rowMeans(dogs_shelter_preds)), lwd=2, col=scales::alpha("slateblue",0.8))
    axis(side = 1, at = c(-5, 1, 2, 3, 5))
    axis(side = 2, at = c(0, 1, 1.6), labels = c(0, 1, ""), las=1)
    
    mtext(
      text = paste0("n = ", sum(!is.na(dogs_shelter_data$behaviour_code_colour)), 
                    "; missing = ", round(sum(is.na(dogs_shelter_data$behaviour_code_colour))/nrow(dogs_shelter_data),2)*100, 
                    "%"
      ),
      side = 3, at = -2, cex = 1
    )
    
    # adoption code probs plot
    adoption_code_probs <- table(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "adoption", "codes"])/
      length(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "adoption", "codes"])
    
    which_codes <- as.numeric(names(adoption_code_probs))
    
    if(length(which_codes) != 3){
      missing_codes <- numeric()
      for(j in 1:3){
        if(length(which_codes[which_codes %in% j]) == 0) {
          missing_codes <- c(missing_codes, j)
        }
      }
      
      temp_codes <- rep(NA, 3)
      for(j in 1:3){
        if(j %in% which_codes) temp_codes[j] <- adoption_code_probs[names(adoption_code_probs)==j] else temp_codes[j] <- 0
      }
      
      adoption_code_probs <- temp_codes
      names(adoption_code_probs) <- 1:3
    }
    
    # ordinal code probs
    adoption_pred_codes <- do.call("rbind.data.frame", 
                                   apply(dogs_adoption_preds, 2, function(x) get_ordinal_probs(x, draws$epsilon)))
    
    plot(density(rowMeans(dogs_adoption_preds)), lwd=2,
         xlim=c(-5,5), ylim=c(0,1.6), bty="n", axes=F, main="", type="n",
         xlab = "behaviour", 
         ylab = "probability & density", 
         cex.lab = label_size)
    for(l in 1:3) if(adoption_code_probs[l] > 0) lines(c(l,l), c(0, adoption_code_probs[l]), lwd=10, col=scales::alpha(the_colours[l],0.6) )
    points(1:3, colMeans(adoption_pred_codes))
    segments(x0 = 1:3, x1 = 1:3, 
             y0 = apply(adoption_pred_codes, 2, hdi)[1,], 
             y1 = apply(adoption_pred_codes, 2, hdi)[2,]
    )
    lines(density(rowMeans(dogs_adoption_preds)), lwd=2, col=scales::alpha("slateblue",0.8))
    axis(side = 1, at = c(-5, 1, 2, 3, 5))
    axis(side = 2, at = c(0, 1, 1.6), labels = c(0, 1, ""), las=1)
    
    mtext(
      text = paste0("n = ", sum(!is.na(dogs_adoption_data$behaviour_code_colour)), 
                    "; missing = ", round(sum(is.na(dogs_adoption_data$behaviour_code_colour))/nrow(dogs_adoption_data),2)*100, 
                    "%"
      ),
      side = 3, at = -2, cex = 1
    )
    dev.off()
}#i

par(mfrow=c(1,1))

# # Figure 4
# png(paste0(figure_path, "/figure_4.png"), width=1400, height=1400, res=200)
# set.seed(2020)
# n_dogs <- 8
# sample_dogs <- sample(1:length(unique(d_s_stan$dog_id)), n_dogs, replace=F)
# which_contexts <- context_labels_numbers #1:8 #sample(1:8, n_dogs, replace=T)
# op <- par(mfrow = c(n_dogs/2,2),
#           oma = c(5,4,0,0) + 0.1,
#           mar = c(2,0,1,1) + 0.1)
# 
# for(i in seq_along(sample_dogs)){
#   dogs_shelter_data <- d_s_stan[d_s_stan$dog_id %in% sample_dogs[i] & d_s_stan$context_id %in% which_contexts[i], 
#                                 c("day_since_arrival_Z", "day_since_arrival", "behaviour_code_colour")]
#   
#   dogs_adoption_data <- d_a_stan[d_a_stan$dog_id %in% sample_dogs[i] & d_a_stan$context_id %in% which_contexts[i], 
#                                  c("days_after_adoption_Z", "days_after_adoption", "behaviour_code_colour")]
#   
#   total_length <- nrow(dogs_shelter_data) + nrow(dogs_adoption_data)
#   
#   dogs_shelter_preds <- sapply(dogs_shelter_data$day_since_arrival_Z, 
#                                function(x){
#                                  draws$alpha + 
#                                    random_effect_list$r_s_j_intercept[,sample_dogs[i]] + 
#                                    random_effect_list$r_s_g_intercept[,which_contexts[i]] + 
#                                    random_effect_list$r_s_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[i])] + 
#                                    (
#                                      draws$beta + 
#                                        random_effect_list$r_s_j_slope[,sample_dogs[i]] + 
#                                        random_effect_list$r_s_g_slope[,which_contexts[i]] + 
#                                        random_effect_list$r_s_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[i])]
#                                    ) * x
#                                })
#   
#   dogs_adoption_preds <- sapply(dogs_adoption_data$days_after_adoption_Z, 
#                                 function(x){
#                                   draws$delta + 
#                                     random_effect_list$r_a_j_intercept[,sample_dogs[i]] + 
#                                     random_effect_list$r_a_g_intercept[,which_contexts[i]] + 
#                                     random_effect_list$r_a_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[i])] + 
#                                     (
#                                       draws$gamma + 
#                                         random_effect_list$r_a_j_slope[,sample_dogs[i]] + 
#                                         random_effect_list$r_a_g_slope[,which_contexts[i]] + 
#                                         random_effect_list$r_a_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[i])]
#                                     ) * x
#                                 })
#   
#   
#   dog_predict_df <- data.frame(
#     dog_id = rep(sample_dogs[i], total_length),
#     which_context = rep(which_contexts[i], total_length),
#     days = c(dogs_shelter_data$day_since_arrival, max(dogs_shelter_data$day_since_arrival) + dogs_adoption_data$days_after_adoption),
#     codes = c(ifelse(is.na(dogs_shelter_data$behaviour_code_colour), 0, dogs_shelter_data$behaviour_code_colour), 
#               ifelse(is.na(dogs_adoption_data$behaviour_code_colour), 0, dogs_adoption_data$behaviour_code_colour)),
#     time = rep(c("shelter","adoption"), c(nrow(dogs_shelter_data), nrow(dogs_adoption_data))),
#     mu = c(apply(dogs_shelter_preds, 2, mean), apply(dogs_adoption_preds, 2, mean)),
#     hdi_2.5 = c(apply(dogs_shelter_preds, 2, hdi)[1,], apply(dogs_adoption_preds, 2, hdi)[1,]),
#     hdi_97.5 = c(apply(dogs_shelter_preds, 2, hdi)[2,], apply(dogs_adoption_preds, 2, hdi)[2,])
#   )
#   
#   dog_predict_df$code_colour <- ifelse(
#     dog_predict_df$codes %in% 0, scales::alpha("black", 0.1), 
#     ifelse(
#       dog_predict_df$codes %in% 1, "darkolivegreen4", 
#       ifelse(
#         dog_predict_df$codes %in% 2, "darkorange2",
#         "red"
#       )
#     )
#   )
#   
#   plot(dog_predict_df$days, dog_predict_df$codes, type="n", ylim=c(-5,5), axes=F, bty="n",
#        #xlab = "observation days", ylab = "behaviour", 
#        cex.lab = 1)
#   abline(h=c(1, 2, 3), lty=3, 
#          col=scales::alpha("black", 0.5))
#   polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5, 
#               min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5,
#               max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1, 
#               max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1), 
#           y = c(-10, 10, 10, -10), 
#           col = scales::alpha("black", 0.1), border = NA
#   )
#   polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1, 
#               min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1,
#               max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5, 
#               max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5), 
#           y = c(-10, 10, 10, -10), 
#           col = scales::alpha("black", 0.1), border = NA
#   )
#   if(i %in% seq(1,20,2)){
#     axis(side = 2, cex.axis=1, las=1,
#          at = c(-5, 1,2,3,5))
#   }
#   axis(side=1, cex.axis=1,
#        at = c(0, 
#               max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"]),
#               min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"]),
#               max(dog_predict_df$days))
#   )
#   
#   if(i == 0){
#     mtext("shelter", line = 0, at = mean(dog_predict_df[dog_predict_df$time %in% "shelter", "days"]), cex = 0.8)
#     mtext("adoption", line=0, at = mean(dog_predict_df[dog_predict_df$time %in% "adoption", "days"]), cex = 0.8)
#   }
#   
#    mtext(paste0(sample_dogs[i], " ", clean_context_labels[which(context_labels_numbers %in% which_contexts[i])]),
#          cex = 0.7, adj = 0)
#   
#   # shelter preds
#   sample_lines <- sample(1:nrow(draws), 50, replace=T)
#   for(i in seq_along(sample_lines)){
#     lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"],
#           dogs_shelter_preds[sample_lines[i], ],
#           col = scales::alpha("slateblue", 0.2))
#   } 
#   for(i in seq_along(sample_lines)){
#     lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"],
#           dogs_adoption_preds[sample_lines[i], ],
#           col = scales::alpha("slateblue", 0.2))
#   }
#    
#   # lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "shelter", "mu"], 
#   #       lwd=2)
#   # lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "shelter", "hdi_2.5"], 
#   #       lwd=2, lty=2)
#   # lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "shelter", "hdi_97.5"], 
#   #       lwd=2, lty=2)
#   # adoption preds
#   # lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "adoption", "mu"], 
#   #       lwd=2)
#   # lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "adoption", "hdi_2.5"], 
#   #       lwd=2, lty=2)
#   # lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"], 
#   #       dog_predict_df[dog_predict_df$time %in% "adoption", "hdi_97.5"], 
#   #       lwd=2, lty=2)
#   
#   points(dog_predict_df$days, ifelse(is.na(dog_predict_df$codes), 0, dog_predict_df$codes), 
#          pch=ifelse(dog_predict_df$codes > 0, 16, 1), 
#          cex = ifelse(dog_predict_df$codes > 0, 1, 0),
#          col = scales::alpha(dog_predict_df$code_colour, 1))
#   
# }
# title(
#   xlab = "observation days (shelter & rehoming)",
#   ylab = "behaviour (ordinal & latent scales)",
#   outer = TRUE, line = 2, 
#   cex.lab = 2
#   )
# par(op)
# par(mfrow=c(1,1))
# dev.off()


#-------------------------------------------------------------------------------------------------------
# raw data - for each context, count up fp, fn, tp, tn for chosen code
# this is a CRUDE estimate!
# Should not be trusted because it does not take into account measurement error, differing amounts of data etc.
diagnostics_vec <- c("fp", "fn", "tp", "tn")
diagnostics_green <- diagnostics_amber <- diagnostics_red <- array(NA, dim=c(4, 8), dimnames = list(c(diagnostics_vec), c(old_context_labels)))
for(k in 1:length(unique(d_s_stan$context_id))){
  for(d in 1:length(diagnostics_vec)){
    diagnostics_green[d, k] <- sum( 
      get_diagnostic(d1 = d_s_stan[d_s_stan$context_id %in% k, ], 
                     d2 = d_a_stan[d_a_stan$context_id %in% k, ], 
                     code_colour = 1, 
                     type = diagnostics_vec[d]
      ) 
    )
    diagnostics_amber[d, k] <- sum( 
      get_diagnostic(d1 = d_s_stan[d_s_stan$context_id %in% k, ], 
                     d2 = d_a_stan[d_a_stan$context_id %in% k, ], 
                     code_colour = 2, 
                     type = diagnostics_vec[d]
      ) 
    )
    diagnostics_red[d, k] <- sum( 
      get_diagnostic(d1 = d_s_stan[d_s_stan$context_id %in% k, ], 
                     d2 = d_a_stan[d_a_stan$context_id %in% k, ], 
                     code_colour = 3, 
                     type = diagnostics_vec[d]
      ) 
    )
  }
}

print( ppvs_green <- apply(diagnostics_green, 2, function(x) ppv(x["tp"], x["fp"])) )
print( npvs_green <- apply(diagnostics_green, 2, function(x) npv(x["tn"], x["fn"])) )
print( ppvs_amber <- apply(diagnostics_amber, 2, function(x) ppv(x["tp"], x["fp"])) )
print( npvs_amber <- apply(diagnostics_amber, 2, function(x) npv(x["tn"], x["fn"])) )
print( ppvs_red <- apply(diagnostics_red, 2, function(x) ppv(x["tp"], x["fp"])) )
print( npvs_red <- apply(diagnostics_red, 2, function(x) npv(x["tn"], x["fn"])) )


# now make predictions from the model posterior distributions - better method!
green_diff_es <- amber_diff_es <- red_diff_es <- array(NA, 
                                                       dim=c( length(unique(d_s_stan$dog_id)),
                                                              length(unique(d_s_stan$context_id))*4))

# which dogs have different probabilities between shelter and post-adoption
for(i in 1:length(unique(d_s_stan$dog_id))){
  
  # get the day after arrival/adoption for each dog desired for predictions
  # set to zero to hold at overall means (~ 30 days after arrival, ~30 days post-adoption)
  x <- mean( d_s_stan[d_s_stan$dog_id %in% i, "day_since_arrival_Z"] )
  x_joint_1 <- d_a_stan[d_a_stan$dog_id %in% i, "days_after_adoption_Z"][[1]]
  x_joint_2 <- d_a_stan[d_a_stan$dog_id %in% i, "days_after_adoption_Z"][[2]]
  
  for(k in 1:length(unique(d_s_stan$context_id))){
    
    dog_shelter_pred <- draws$alpha + 
      random_effect_list$r_s_j_intercept[,i] + 
      random_effect_list$r_s_g_intercept[,k] + 
      random_effect_list$r_s_jg_intercept[, make_cc_interaction(i,k)] + 
      (draws$beta + 
         random_effect_list$r_s_j_slope[,i] + 
         random_effect_list$r_s_g_slope[,k] + 
         random_effect_list$r_s_jg_slope[, make_cc_interaction(i,k)]) * x
    
    dog_adoption_pred_1 <- draws$delta + 
      random_effect_list$r_a_j_intercept[,i] + 
      random_effect_list$r_a_g_intercept[,k] + 
      random_effect_list$r_a_jg_intercept[, make_cc_interaction(i,k)] +
      (draws$gamma + 
         random_effect_list$r_a_j_slope[,i] + 
         random_effect_list$r_a_g_slope[,k] + 
         random_effect_list$r_a_jg_slope[, make_cc_interaction(i,k)]) * x_joint_1
    
    dog_adoption_pred_2 <- draws$delta + 
      random_effect_list$r_a_j_intercept[,i] + 
      random_effect_list$r_a_g_intercept[,k] + 
      random_effect_list$r_a_jg_intercept[, make_cc_interaction(i,k)] +
      (draws$gamma + 
         random_effect_list$r_a_j_slope[,i] + 
         random_effect_list$r_a_g_slope[,k] + 
         random_effect_list$r_a_jg_slope[, make_cc_interaction(i,k)]) * x_joint_2
    
    dog_green_diffs1 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[1]] - 
                          get_ordinal_probs(dog_adoption_pred_1, draws$epsilon)[[1]]
    )
    dog_amber_diffs1 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[2]] - 
                          get_ordinal_probs(dog_adoption_pred_1, draws$epsilon)[[2]]
    )
    dog_red_diffs1 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[3]] - 
                        get_ordinal_probs(dog_adoption_pred_1, draws$epsilon)[[3]]
    )
    dog_green_diffs2 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[1]] - 
                           get_ordinal_probs(dog_adoption_pred_2, draws$epsilon)[[1]]
    )
    dog_amber_diffs2 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[2]] - 
                           get_ordinal_probs(dog_adoption_pred_2, draws$epsilon)[[2]]
    )
    dog_red_diffs2 <- (get_ordinal_probs(dog_shelter_pred, draws$sigma)[[3]] - 
                         get_ordinal_probs(dog_adoption_pred_2, draws$epsilon)[[3]]
    )
    
    green_diff_es[i,k] <- ifelse(hdi(dog_green_diffs1)[1] > 0, mean(dog_green_diffs1), NA)
    green_diff_es[i,k+8] <- ifelse(hdi(dog_green_diffs1)[2] < 0, mean(dog_green_diffs1), NA)
    amber_diff_es[i,k] <- ifelse(hdi(dog_amber_diffs1)[1] > 0, mean(dog_amber_diffs1), NA)
    amber_diff_es[i,k+8] <- ifelse(hdi(dog_amber_diffs1)[2] < 0, mean(dog_amber_diffs1), NA)
    red_diff_es[i,k] <- ifelse(hdi(dog_red_diffs1)[1] > 0, mean(dog_red_diffs1), NA)
    red_diff_es[i,k+8] <- ifelse(hdi(dog_red_diffs1)[2] < 0, mean(dog_red_diffs1), NA)
    
    green_diff_es[i,k+8*2] <- ifelse(hdi(dog_green_diffs2)[1] > 0, mean(dog_green_diffs2), NA)
    green_diff_es[i,k+8*3] <- ifelse(hdi(dog_green_diffs2)[2] < 0, mean(dog_green_diffs2), NA)
    amber_diff_es[i,k+8*2] <- ifelse(hdi(dog_amber_diffs2)[1] > 0, mean(dog_amber_diffs2), NA)
    amber_diff_es[i,k+8*3] <- ifelse(hdi(dog_amber_diffs2)[2] < 0, mean(dog_amber_diffs2), NA)
    red_diff_es[i,k + 8*2] <- ifelse(hdi(dog_red_diffs2)[1] > 0, mean(dog_red_diffs2), NA)
    red_diff_es[i,k+8*3] <- ifelse(hdi(dog_red_diffs2)[2] < 0, mean(dog_red_diffs2), NA)
  }
}

print(
  colour_diffs_df <- data.frame(
    colour = rep(c("green","amber","red"), each=8*4),
    higher_lower_s = rep(c("higher","lower"), each=8, times=3*2),
    adoption_call = rep(c("first","second"), each = 8*2, times = 3),
    contexts_num = rep(sorted_context_numbers, 4*3),
    context_label = rep(old_context_labels, 4*3),
    num_dogs = c(
      apply(green_diff_es[,1:8], 2, function(x) sum( !is.na(x) )),
      -1*apply(green_diff_es[,9:16], 2, function(x) sum(!is.na(x))),
      apply(green_diff_es[,17:24], 2, function(x) sum( !is.na(x) )),
      -1*apply(green_diff_es[,25:32], 2, function(x) sum(!is.na(x))),
      apply(amber_diff_es[,1:8], 2, function(x) sum(!is.na(x))),
      -1*apply(amber_diff_es[,9:16], 2, function(x) sum(!is.na(x))),
      apply(amber_diff_es[,17:24], 2, function(x) sum( !is.na(x) )),
      -1*apply(amber_diff_es[,25:32], 2, function(x) sum(!is.na(x))),
      apply(red_diff_es[,1:8], 2, function(x) sum(!is.na(x))),
      -1*apply(red_diff_es[,9:16], 2, function(x) sum(!is.na(x))),
      apply(red_diff_es[,17:24], 2, function(x) sum( !is.na(x) )),
      -1*apply(red_diff_es[,25:32], 2, function(x) sum(!is.na(x)))
    ),
    mu = c(
      apply(green_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(green_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(green_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(green_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(amber_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(amber_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(amber_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(amber_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(red_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(red_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(red_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA),
      apply(red_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) mean(x, na.rm = T) else NA)
    ),
    
    hdi_2.5 = c(
      apply(green_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(green_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(green_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(green_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(amber_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(amber_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(amber_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(amber_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(red_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(red_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(red_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA),
      apply(red_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[1] else NA)
    ),
    
    hdi_97.5 = c(
      apply(green_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(green_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(green_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(green_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(amber_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(amber_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(amber_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(amber_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(red_diff_es[,1:8], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(red_diff_es[,9:16], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(red_diff_es[,17:24], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA),
      apply(red_diff_es[,25:32], 2, function(x) if(sum(!is.na(x))>1) hdi(x)[2] else NA)
    )
  ))

apply(red_diff_es, 2, function(x) which(!is.na(x))) # which dogs differ

colour_diffs_df_ordered <- colour_diffs_df[order(colour_diffs_df$contexts_num), ]

# figure 5 - depending on days chosen
png(paste0(figure_path, "/figure_5.png"), width=1400, height=1200, res=200)
barplot_axis_names <- as.vector(sapply(clean_context_labels, function(x) c(x, " ")))
op <- par(mfrow = c(3, 2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

green1_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                              colour_diffs_df_ordered$adoption_call %in% "first", "num_dogs"], nrow=2)
green1_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                        colour_diffs_df_ordered$adoption_call %in% "first", "mu"]), na.rm = T), 2)
green1_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                           colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[1], 2)
green1_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                             colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[2], 2)

green2_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                              colour_diffs_df_ordered$adoption_call %in% "second", "num_dogs"], nrow=2)
green2_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                        colour_diffs_df_ordered$adoption_call %in% "second", "mu"]), na.rm = T), 2)
green2_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                           colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[1], 2)
green2_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "green"&
                                                             colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[2], 2)

barplot( green1_mat/241, beside = T, ylim=c(-1,1), col = "darkolivegreen2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
gbp1 <- barplot( green1_mat/241, beside = T, ylim=c(-1,1), col = "darkolivegreen2", axes=F, add = T)
mtext(at=1, text="a", font=2)
axis(side=2, at=c(-1, 0, 1), lwd=1, las=2)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = gbp1, y=as.vector(green1_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(green1_mat)))
# mtext(text=paste0("effect size: ", green1_mu_es, " [", green1_hdi_es2.5, ", ", green1_hdi_es_97.5, "]"),
#       side = 3, line = -1.5, at = 7, cex = 0.8)
mtext(text=bquote(expr = delta == .(green1_mu_es) ~ "["*.(green1_hdi_es2.5)*","*.(green1_hdi_es_97.5)*"]"),
      side = 3, line = -1.5, at = 7, cex = 0.8)

barplot( green2_mat/241, beside = T, ylim=c(-1,1), col = "darkolivegreen2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
gbp2 <- barplot( green2_mat/241, beside = T, ylim=c(-1,1), col = "darkolivegreen2", axes=F, add = T)
mtext(at=1, text="b", font=2)
#axis(side=2, at=c(-1, 0, 1), lwd=1)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = gbp2, y=as.vector(green2_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(green2_mat)))
mtext(text=bquote(expr = delta == .(green2_mu_es) ~ "["*.(green2_hdi_es2.5)*","*.(green2_hdi_es_97.5)*"]"),
      side = 3, line = -1.5, at = 7, cex = 0.8)

amber1_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                              colour_diffs_df_ordered$adoption_call %in% "first", "num_dogs"], nrow=2)
amber1_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                        colour_diffs_df_ordered$adoption_call %in% "first", "mu"]), na.rm = T), 2)
amber1_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                           colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[1], 2)
amber1_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                               colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[2], 2)

amber2_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                               colour_diffs_df_ordered$adoption_call %in% "second", "num_dogs"], nrow=2)
amber2_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                         colour_diffs_df_ordered$adoption_call %in% "second", "mu"]), na.rm = T), 2)
amber2_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                            colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[1], 2)
amber2_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "amber"&
                                                              colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[2], 2)

barplot( amber1_mat/241, beside = T, ylim=c(-1,1), col = "darkorange2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
abp1 <- barplot( amber1_mat/241, beside = T, ylim=c(-1,1), col = "darkorange2", axes=F, add = T)
mtext(at=1, text="c", font=2)
axis(side=2, at=c(-1, 0, 1), lwd=1, las=2)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = abp1, y=as.vector(amber1_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(amber1_mat)))
mtext(text=bquote(expr = delta == "NA"),
      side = 3, line = -1.5, at = 4, cex = 0.8)

barplot( amber2_mat/241, beside = T, ylim=c(-1,1), col = "darkorange2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
abp2 <- barplot( amber2_mat/241, beside = T, ylim=c(-1,1), col = "darkorange2", axes=F, add = T)
mtext(at=1, text="d", font=2)
#axis(side=2, at=c(-1, 0, 1), lwd=1)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = abp2, y=as.vector(amber2_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(amber2_mat)))
mtext(text=bquote(expr = delta == "NA"),
      side = 3, line = -1.5, at = 4, cex = 0.8)


red1_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                               colour_diffs_df_ordered$adoption_call %in% "first", "num_dogs"], nrow=2)
red1_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                         colour_diffs_df_ordered$adoption_call %in% "first", "mu"]), na.rm = T), 2)
red1_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                            colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[1], 2)
red1_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                              colour_diffs_df_ordered$adoption_call %in% "first", "mu"]))[2], 2)

red2_mat <- matrix(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                               colour_diffs_df_ordered$adoption_call %in% "second", "num_dogs"], nrow=2)
red2_mu_es <- round(mean(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                         colour_diffs_df_ordered$adoption_call %in% "second", "mu"]), na.rm = T), 2)
red2_hdi_es2.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                            colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[1], 2)
red2_hdi_es_97.5 <- round(hdi(abs(colour_diffs_df_ordered[colour_diffs_df_ordered$colour %in% "red"&
                                                              colour_diffs_df_ordered$adoption_call %in% "second", "mu"]))[2], 2)

barplot( red1_mat/241, beside = T, ylim=c(-1,1), col = "firebrick2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
rbp1 <- barplot( red1_mat/241, beside = T, ylim=c(-1,1), col = "firebrick2", axes=F, add = T)
mtext(at=1, text="e", font=2)
axis(side=2, at=c(-1, 0, 1), lwd=1, las=2)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = rbp1, y=as.vector(red1_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(red1_mat)))
mtext(text=bquote(expr = delta == .(red1_mu_es) ~ "["*.(red1_hdi_es2.5)*","*.(red1_hdi_es_97.5)*"]"),
      side = 3, line = -1.5, at = 7, cex = 0.8)
text(as.vector(rbp1), par("usr")[3] + 0.42, labels = barplot_axis_names, cex = 1, srt = 45, pos = 1, xpd = TRUE)

barplot( red2_mat/241, beside = T, ylim=c(-1,1), col = "firebrick2", axes=F)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = scales::alpha("lightgray", 0.6), border = NA)
rbp2 <- barplot( red2_mat/241, beside = T, ylim=c(-1,1), col = "firebrick2", axes=F, add = T)
mtext(at=1, text="f", font=2)
#axis(side=2, at=c(-1, 0, 1), lwd=1)
abline(h=0, col=scales::alpha("black",0.8), lty=1)
text(x = rbp2, y=as.vector(red2_mat/241) + 0.2*rep(c(1,-1), 8), labels = abs(as.vector(red2_mat)))
mtext(text=bquote(expr = delta == .(red2_mu_es) ~ "["*.(red2_hdi_es2.5)*","*.(red2_hdi_es_97.5)*"]"),
      side = 3, line = -1.5, at = 7, cex = 0.8)
text(as.vector(rbp2), par("usr")[3] + 0.42, labels = barplot_axis_names, cex = 1, srt = 45, pos = 1, xpd = TRUE)

title(ylab = "proportion of dogs", xlab = "contexts", outer = T, line = 3, cex.lab = 1.7)
par(op)
par(mfrow=c(1,1))
dev.off()



avg_int <- apply(random_effect_list$r_s_j_intercept, 2, mean)
avg_slp <- apply(random_effect_list$r_s_j_slope, 2, mean)
avg_int_missing_a <- apply(random_effect_list$r_a_j_missing, 2, mean)
plot(avg_int, avg_int_missing_a, col=scales::alpha("black",0.3))
plot(avg_slp, avg_int_missing_a, col=scales::alpha("black",0.3))

#---------------------------------------------------------------------------------
#PLOT ALL DOGS' POSTERIOR PREDICTIONS

sample_dogs <- 1:241 
which_contexts <- 1:8 
# op <- par(mfrow = c(1,3),
#           oma = c(5,4,0,0) + 0.1,
#           mar = c(0,0,1,1) + 0.1)

for(i in seq_along(sample_dogs)){  
  for(k in seq_along(which_contexts)){
    png(paste0("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/analysis_code/all_dogs_plots/dog_", 
               sample_dogs[i], "_context_", which_contexts[k]), 
        width = 2000, height=800, res=200)
    par(mfrow=c(1,3))
   
    label_size <- 1.5
    dogs_shelter_data <- d_s_stan[d_s_stan$dog_id %in% sample_dogs[i] & d_s_stan$context_id %in% which_contexts[k], 
                                  c("day_since_arrival_Z", "day_since_arrival", "behaviour_code_colour")]
    
    dogs_adoption_data <- d_a_stan[d_a_stan$dog_id %in% sample_dogs[i] & d_a_stan$context_id %in% which_contexts[k], 
                                   c("days_after_adoption_Z", "days_after_adoption", "behaviour_code_colour")]
    
    total_length <- nrow(dogs_shelter_data) + nrow(dogs_adoption_data)
    
    dogs_shelter_preds <- sapply(dogs_shelter_data$day_since_arrival_Z, 
                                 function(x){
                                   draws$alpha + 
                                     random_effect_list$r_s_j_intercept[,sample_dogs[i]] + 
                                     random_effect_list$r_s_g_intercept[,which_contexts[k]] + 
                                     random_effect_list$r_s_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[k])] + 
                                     (
                                       draws$beta + 
                                         random_effect_list$r_s_j_slope[,sample_dogs[i]] + 
                                         random_effect_list$r_s_g_slope[,which_contexts[k]] + 
                                         random_effect_list$r_s_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[k])]
                                     ) * x
                                 })
    
    dogs_adoption_preds <- sapply(dogs_adoption_data$days_after_adoption_Z, 
                                  function(x){
                                    draws$delta + 
                                      random_effect_list$r_a_j_intercept[,sample_dogs[i]] + 
                                      random_effect_list$r_a_g_intercept[,which_contexts[k]] + 
                                      random_effect_list$r_a_jg_intercept[, make_cc_interaction(sample_dogs[i],which_contexts[k])] + 
                                      (
                                        draws$gamma + 
                                          random_effect_list$r_a_j_slope[,sample_dogs[i]] + 
                                          random_effect_list$r_a_g_slope[,which_contexts[k]] + 
                                          random_effect_list$r_a_jg_slope[, make_cc_interaction(sample_dogs[i],which_contexts[k])]
                                      ) * x
                                  })
    
    
    dog_predict_df <- data.frame(
      dog_id = rep(sample_dogs[i], total_length),
      which_context = rep(which_contexts[i], total_length),
      days = c(dogs_shelter_data$day_since_arrival, max(dogs_shelter_data$day_since_arrival) + dogs_adoption_data$days_after_adoption),
      codes = c(ifelse(is.na(dogs_shelter_data$behaviour_code_colour), 0, dogs_shelter_data$behaviour_code_colour), 
                ifelse(is.na(dogs_adoption_data$behaviour_code_colour), 0, dogs_adoption_data$behaviour_code_colour)),
      time = rep(c("shelter","adoption"), c(nrow(dogs_shelter_data), nrow(dogs_adoption_data))),
      mu = c(apply(dogs_shelter_preds, 2, mean), apply(dogs_adoption_preds, 2, mean)),
      hdi_2.5 = c(apply(dogs_shelter_preds, 2, hdi)[1,], apply(dogs_adoption_preds, 2, hdi)[1,]),
      hdi_97.5 = c(apply(dogs_shelter_preds, 2, hdi)[2,], apply(dogs_adoption_preds, 2, hdi)[2,])
    )
    
    dog_predict_df$code_colour <- ifelse(
      dog_predict_df$codes %in% 0, scales::alpha("black", 0.1), 
      ifelse(
        dog_predict_df$codes %in% 1, "darkolivegreen3", 
        ifelse(
          dog_predict_df$codes %in% 2, "darkorange2",
          "firebrick2"
        )
      )
    )
    
    plot(dog_predict_df$days, dog_predict_df$codes, type="n", ylim=c(-5,5), 
         axes=F, bty="n",
         xlab = "observation days", ylab = "behaviour (latent & ordinal scales)", 
         cex.lab = label_size)
    abline(h=c(1, 2, 3), lty=3, 
           col=scales::alpha("black", 0.5))
    polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5, 
                min(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])-5,
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1, 
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"])+1), 
            y = c(-10, 10, 10, -10), 
            col = scales::alpha("black", 0.1), border = NA
    )
    polygon(x=c(min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1, 
                min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])-1,
                max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5, 
                max(dog_predict_df[dog_predict_df$time %in% "adoption", "days"])+5), 
            y = c(-10, 10, 10, -10), 
            col = scales::alpha("black", 0.1), border = NA
    )
    axis(side = 2, cex.axis=1, las=1,
         at = c(-5, 1,2,3,5))
    axis(side=1, cex.axis=1,
         at = c(0, 
                max(dog_predict_df[dog_predict_df$time %in% "shelter", "days"]),
                min(dog_predict_df[dog_predict_df$time %in% "adoption", "days"]),
                max(dog_predict_df$days))
    )
    
    mtext(paste0(sample_dogs[i], " ", clean_context_labels[which(context_labels_numbers %in% which_contexts[k])]),
          cex = 0.7, adj = 0)
    
    # shelter preds
    sample_lines <- sample(1:nrow(draws), 50, replace=T)
    for(j in seq_along(sample_lines)){
      lines(dog_predict_df[dog_predict_df$time %in% "shelter", "days"],
            dogs_shelter_preds[sample_lines[j], ],
            col = scales::alpha("slateblue", 0.2))
    } 
    for(j in seq_along(sample_lines)){
      lines(dog_predict_df[dog_predict_df$time %in% "adoption", "days"],
            dogs_adoption_preds[sample_lines[j], ],
            col = scales::alpha("slateblue", 0.2))
    }
    
    points(dog_predict_df$days, ifelse(is.na(dog_predict_df$codes), 0, dog_predict_df$codes), 
           pch=ifelse(dog_predict_df$codes > 0, 16, 1), 
           cex = ifelse(dog_predict_df$codes > 0, 1, 0),
           col = scales::alpha(dog_predict_df$code_colour, 1))
    
    # shelter code probs plot
    shelter_code_probs <- table(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "shelter", "codes"])/
      length(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "shelter", "codes"])
    
    which_codes <- as.numeric(names(shelter_code_probs))
    
    if(length(which_codes) != 3){
      missing_codes <- numeric()
      for(j in 1:3){
        if(length(which_codes[which_codes %in% j]) == 0) {
          missing_codes <- c(missing_codes, j)
          #which_codes <- sort(c(which_codes, i), decreasing = F)
        }
      }
      
      temp_codes <- rep(NA, 3)
      for(j in 1:3){
        if(j %in% which_codes) temp_codes[j] <- shelter_code_probs[names(shelter_code_probs)==j] else temp_codes[j] <- 0
      }
      
      shelter_code_probs <- temp_codes
      names(shelter_code_probs) <- 1:3
    }
    
    # ordinal code probs
    shelter_pred_codes <- do.call("rbind.data.frame", apply(dogs_shelter_preds, 2, function(x) get_ordinal_probs(x, draws$sigma)))
    
    the_colours <- c("darkolivegreen3","darkorange2","firebrick2")
    
    plot(density(rowMeans(dogs_shelter_preds)), lwd=2, col=scales::alpha("slateblue",0.8), 
         xlim=c(-5,5), ylim=c(0,1.6), bty="n", axes=F, main="", type="n",
         xlab = "behaviour (latent & ordinal scales)", 
         ylab = "probability & density", 
         cex.lab = label_size)
    for(l in 1:3) if(shelter_code_probs[l] > 0) lines(c(l,l), c(0, shelter_code_probs[l]), lwd=10, col=the_colours[l] )
    points(1:3, colMeans(shelter_pred_codes))
    segments(x0 = 1:3, x1 = 1:3, 
             y0 = apply(shelter_pred_codes, 2, hdi)[1,], 
             y1 = apply(shelter_pred_codes, 2, hdi)[2,]
    )
    lines(density(rowMeans(dogs_shelter_preds)), lwd=2, col=scales::alpha("slateblue",0.8))
    axis(side = 1, at = c(-5, 1, 2, 3, 5))
    axis(side = 2, at = c(0, 1, 1.6), labels = c(0, 1, ""), las=1)
    
    mtext(
      text = paste0("n = ", sum(!is.na(dogs_shelter_data$behaviour_code_colour)), 
                    "; missing = ", round(sum(is.na(dogs_shelter_data$behaviour_code_colour))/nrow(dogs_shelter_data),2)*100, 
                    "%"
      ),
      side = 3, at = -2, cex = 0.8
    )
    
    # adoption code probs plot
    adoption_code_probs <- table(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "adoption", "codes"])/
      length(dog_predict_df[dog_predict_df$codes %in% 1:3 & dog_predict_df$time %in% "adoption", "codes"])
    
    which_codes <- as.numeric(names(adoption_code_probs))
    
    if(length(which_codes) != 3){
      missing_codes <- numeric()
      for(j in 1:3){
        if(length(which_codes[which_codes %in% j]) == 0) {
          missing_codes <- c(missing_codes, j)
        }
      }
      
      temp_codes <- rep(NA, 3)
      for(j in 1:3){
        if(j %in% which_codes) temp_codes[j] <- adoption_code_probs[names(adoption_code_probs)==j] else temp_codes[j] <- 0
      }
      
      adoption_code_probs <- temp_codes
      names(adoption_code_probs) <- 1:3
    }
    
    # ordinal code probs
    adoption_pred_codes <- do.call("rbind.data.frame", apply(dogs_adoption_preds, 2, function(x) get_ordinal_probs(x, draws$epsilon)))
    
    plot(density(rowMeans(dogs_adoption_preds)), lwd=2,
         xlim=c(-5,5), ylim=c(0,1.6), bty="n", axes=F, main="", type="n",
         xlab = "behaviour (latent & ordinal scales)", 
         ylab = "probability & density", 
         cex.lab = label_size)
    for(l in 1:3) if(adoption_code_probs[l] > 0) lines(c(l,l), c(0, adoption_code_probs[l]), lwd=10, col=the_colours[l] )
    points(1:3, colMeans(adoption_pred_codes))
    segments(x0 = 1:3, x1 = 1:3, 
             y0 = apply(adoption_pred_codes, 2, hdi)[1,], 
             y1 = apply(adoption_pred_codes, 2, hdi)[2,]
    )
    lines(density(rowMeans(dogs_adoption_preds)), lwd=2, col=scales::alpha("slateblue",0.8))
    axis(side = 1, at = c(-5, 1, 2, 3, 5))
    axis(side = 2, at = c(0, 1, 1.6), labels = c(0, 1, ""), las=1)
    
    mtext(
      text = paste0("n = ", sum(!is.na(dogs_adoption_data$behaviour_code_colour)), 
                    "; missing = ", round(sum(is.na(dogs_adoption_data$behaviour_code_colour))/nrow(dogs_adoption_data),2)*100, 
                    "%"
      ),
      side = 3, at = -2, cex = 0.8
    )
    dev.off()
  }#k
}#i

par(mfrow=c(1,1))




