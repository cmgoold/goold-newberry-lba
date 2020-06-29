############################################################################
# Post-processing results for Goold & Newberry (2020):
#     - Longitudinal behavioural assessment of shelter dogs predicts 
#       behaviour post-adoption

# Copyright Conor Goold (2020)
# c.goold@leeds.ac.uk
############################################################################

setwd("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba")
source("analysis_code/helper_functions.R")
figure_path <- "~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/paper/figures"

get_needed_packages()

# load the posterior distribution
chain_1 <- rstan::read_stan_csv("analysis_code/jhb_dog_samples1.csv")
chain_2 <- rstan::read_stan_csv("analysis_code/jhb_dog_samples1.csv")

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
  col = scales::alpha("darkorange4",0.5), border = scales::alpha("darkorange4",0.5)
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
png(paste0(figure_path, "/figure_2_a.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, rnorm(length(days_grid_Z), 0, 1), 
     type="n", ylim=c(-5, 5), 
     bty="n", axes = F, 
     xlab="days after arrival", ylab="latent behaviour scale", 
     cex.lab=1.5)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, lwd=3,  cex.axis=1.25)
mtext("a", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)
for(i in seq_along(clean_context_labels)){
  lines(days_grid_Z, 
        apply(
          sapply(days_grid_Z, 
                 function(x) {
                   draws$alpha + random_effect_list$r_s_g_intercept[,context_labels_numbers[i]] + 
                     (draws$beta + random_effect_list$r_s_g_slope[,context_labels_numbers[i]]) * x
                 }
          ), 
          2, mean), 
        lwd=3, 
        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
        )
}
legend("topright", 
       inset = c(0.3, 0),
       bty="n",
       legend = clean_context_labels[1:4], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
       )
legend("topright", 
       inset = c(0, 0),
       bty="n",
       legend = clean_context_labels[5:8], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
)
dev.off()

# Figure 2b
png(paste0(figure_path, "/figure_2_b.png"), width=1000, height=1000, res=200)
plot(days_grid_Z, rnorm(length(days_grid_Z), 0, 1), 
     type="n", ylim=c(0, 1.2), 
     bty="n", axes = F, 
     xlab="days after arrival", ylab="probability missing", 
     cex.lab=1.5)
axis(side = 1, at = c(min(days_grid_Z), median(days_grid_Z), max(days_grid_Z)), 
     labels = round(c(min(days_grid), median(days_grid), max(days_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, at = c(0, 0.5, 1), lwd=3,  cex.axis=1.25)
mtext("b", line=0, at=min(days_grid_Z), cex = 1.5, font = 2)

for(i in seq_along(clean_context_labels)){
  lines(days_grid_Z, 
        apply(
          sapply(days_grid_Z, 
                 function(x) {
                   inv_logit(draws$alpha_miss + random_effect_list$r_s_g_missing[,context_labels_numbers[i]] + 
                     draws$beta_miss * x
                   )
                 }
          ), 
          2, mean), 
        lwd=3, 
        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
  )
}
legend("topright", 
       inset = c(0.3, -0),
       bty="n",
       legend = clean_context_labels[1:4], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
)
legend("topright", 
       inset = c(0, -0),
       bty="n",
       legend = clean_context_labels[5:8], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
)
dev.off()

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
# FIGURE 3: plot code colour probabilities across days after adoption

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

# Figure 3a
png(filename = paste0(figure_path,"/figure_3_a.png"), width=1000, height=1000, res=200)
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
png(filename = paste0(figure_path,"/figure_3_b.png"), width=1000, height=1000, res=200)
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
  col = scales::alpha("darkorange4",0.5), border = scales::alpha("darkorange4",0.5)
)
lines(days_a_grid_Z, 
      unlist(lapply( adoption_amber_preds, function(x) mean(unlist(lapply(x, function(z) mean(colMeans(z))))))), col="black", lwd=2)
dev.off()

# red codes
png(filename = paste0(figure_path,"/figure_3_c.png"), width=1000, height=1000, res=200)
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
png(filename = paste0(figure_path,"/figure_3_d.png"), width=1000, height=1000, res=200)
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
clean_context_labels <- c("HND", "HOUSE", "OSTD", "FPL", "UFPL", "FOOD", "TOYS", "DOGS")
context_labels_numbers <- c(4, 5, 6, 3, 2, 8, 7, 1)
context_intercepts_ordered <- random_effect_list$r_a_g_intercept[, context_labels_numbers]
context_slopes_ordered <- random_effect_list$r_a_g_slope[, context_labels_numbers]
context_intercepts_missing_ordered <- random_effect_list$r_a_g_missing[,context_labels_numbers]
nudge <- 0.15

# Figure S1
png(paste0(figure_path, "/figure_S2.png"), width=1100, height=1000, res=200)
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

# Figure 4
png(paste0(figure_path, "/figure_4.png"), width=1000, height=1000, res=200)
plot(days_a_grid_Z, rnorm(length(days_a_grid_Z), 0, 1), 
     type="n", ylim=c(-5, 5), 
     bty="n", axes = F, 
     xlab="days after adoption", ylab="latent behaviour scale", 
     cex.lab=1.5)
axis(side = 1, at = c(min(days_a_grid_Z), median(days_a_grid_Z), max(days_a_grid_Z)), 
     labels = round(c(min(days_a_grid), median(days_a_grid), max(days_a_grid))), 
     lwd = 3, cex.axis=1.25)
axis(side = 2, lwd=3,  cex.axis=1.25)
mtext("a", line=0, at=min(days_a_grid_Z), cex = 1.5, font = 2)
for(i in seq_along(clean_context_labels)){
  lines(days_a_grid_Z, 
        apply(
          sapply(days_a_grid_Z, 
                 function(x) {
                   draws$delta + random_effect_list$r_a_g_intercept[,context_labels_numbers[i]] + 
                     (draws$gamma + random_effect_list$r_a_g_slope[,context_labels_numbers[i]]) * x
                 }
          ), 
          2, mean), 
        lwd=3, 
        col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[i]
  )
}
legend("topright", 
       inset = c(0.3, 0),
       bty="n",
       legend = clean_context_labels[1:4], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4]
)
legend("topright", 
       inset = c(0, 0),
       bty="n",
       legend = clean_context_labels[5:8], 
       lwd = 3, 
       col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[5:8]
)
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
# FIGURE 5: correlations

# 5a (personality)
png(paste0(figure_path, "/figure_5_a.png"), width=1000, height=1000, res=200)
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
mtext("a", line=0, at=-1, cex = 1.5, font = 2)
colours_ <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")
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
png(paste0(figure_path, "/figure_5_b.png"), width=1000, height=1000, res=200)
plot(density(
  Rho_j_raw$`Rho_j[2,5]`, 
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
axis(side = 2, at=c(0:6), lwd=3,  cex.axis=1.25)
mtext("b", line=0, at=-1, cex = 1.5, font = 2)
colours_ <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")
abline(v=0, lty=2, col=scales::alpha("black", 0.5))
lines(density(Rho_j_raw$`Rho_j[2,5]`, adjust = 0.5), col=colours_[1], lwd=3)
lines(density(Rho_g_raw$`Rho_g[2,5]`, adjust = 0.5), col=colours_[2], lwd=3)
lines(density(Rho_jg_raw$`Rho_jg[2,5]`, adjust = 0.5), col=colours_[3], lwd=3)
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

#------------------------------------------------------------------------------------------------------
# ANALYSIS OF SHELTER-ADOPTION BEHAVIOUR DIFFERENCES

sample_dogs <- sample(1:length(unique(d_s_stan$dog_id)), 1, replace=F)
which_context <- 1

diff_preds <- rep(list(list()), length(sample_dogs))

for(i in seq_along(sample_dogs)){
  res_df <- matrix(NA, nrow = nrow(draws), ncol=length(days_grid))
  for(j in 1:length(days_grid)){
    shelter_behaviour <- shelter_green_preds[[j]][[i]][, which_context]
    adoption_behaviour <- adoption_green_preds[[j]][[i]][, which_context]
    res_df[,j] <- shelter_behaviour - adoption_behaviour
  }
  diff_preds[[i]] <- res_df
}

plot(1:10, apply(diff_preds[[1]], 2, mean), type="l", ylim=c(-1, 1))
lines(apply(diff_preds[[1]], 2, hdi)[1,], lty=2)
lines(apply(diff_preds[[1]], 2, hdi)[2,], lty=2)


for(i in seq_along(sample_dogs)){
  res_df <- matrix(NA, nrow = nrow(draws), ncol=length(days_grid))
  for(j in 1:length(days_grid)){
    shelter_behaviour <- shelter_pred_list[[j]][[i]][, which_context]
    adoption_behaviour <- adoption_pred_list[[j]][[i]][, which_context]
    res_df[,j] <- shelter_behaviour - adoption_behaviour
  }
  diff_preds[[i]] <- res_df
}

plot(days_a_grid_Z, apply(diff_preds[[1]], 2, mean), type="l", ylim=c(-5, 5))
lines(days_a_grid_Z, apply(diff_preds[[1]], 2, hdi)[1,], lty=2)
lines(days_a_grid_Z, apply(diff_preds[[1]], 2, hdi)[2,], lty=2)
points(d_a_stan[d_a_stan$dog_id %in% sample_dogs & d_a_stan$context_id %in% which_context, "days_after_adoption_Z"], 
       d_a_stan[d_a_stan$dog_id %in% sample_dogs & d_a_stan$context_id %in% which_context, "behaviour_code_colour"]
       )

