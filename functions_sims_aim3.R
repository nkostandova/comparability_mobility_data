single_run <- function(pop_mat = population_matrix, 
                       InC = initial_conditions,
                       params = parameters
) {
  M0 = InC$M
  lower_S0 = InC$S.lower_S0
  upper_S0 = InC$S.upper_S0
  mean_S0 = InC$S.mean_S0
  S0= (1 - runif(length(lower_S0), min = lower_S0, max = upper_S0))*InC$N0
  I0 = InC$I0
  V1R0 = InC$V1R0
  V1C0 = InC$V1C0
  V20 = InC$V20
  incid_deaths0 = InC$incid_deaths0
 # R0 = InC$R0  ## will need to update after a random draw
  N0 = InC$N0
  
  sd.dw_val = params$sd.dw.val
  t = params$t
  mob_data_prob = params$mob_data_prob
  mcv1_cov_dat = params$mcv1_cov_dat
  mcv2_cov_dat = params$mcv2_cov_dat
  sia_cov = params$sia_cov
  births_rate = params$births_rate
  wane_mat = params$wane_mat
  beta_values = params$beta_values
  alpha = params$alpha
  deaths_rate = params$deaths_rate
  measles_cfr = params$measles_cfr
  MR1_time = params$MR1_time
  MR2_time = params$MR2_time
  p_eff_V1r = params$p_eff_V1r
  p_eff_V1c = params$p_eff_V1c
  intro_mat = params$intro_mat
  
  mat_matrix = matrix(NA,dim(pop_mat)[1],t)
  suscep_matrix <- matrix(NA,dim(pop_mat)[1],t)
  infect_matrix <- matrix(NA,dim(pop_mat)[1],t)
  V1R_matrix = matrix(NA,dim(pop_mat)[1],t)
  V1C_matrix = matrix(NA,dim(pop_mat)[1],t)
  V2_matrix = matrix(NA,dim(pop_mat)[1],t)
  R_matrix = matrix(NA,dim(pop_mat)[1],t)
  incid_deaths_matrix = matrix(NA,dim(pop_mat)[1],t)
  N_matrix = matrix(NA, dim(pop_mat)[1],t)
  
  rownames(mat_matrix) <- rownames(pop_mat)
  rownames(suscep_matrix) <- rownames(pop_mat)
  rownames(infect_matrix) <- rownames(pop_mat)
  rownames(V1R_matrix) <- rownames(pop_mat)
  rownames(V1C_matrix) <- rownames(pop_mat)
  rownames(V2_matrix) <- rownames(pop_mat)
  rownames(incid_deaths_matrix) <- rownames(pop_mat)
  rownames(R_matrix) <- rownames(pop_mat)
  rownames(N_matrix) = rownames(pop_mat)
  
  mat_matrix[,1] = M0
  # for susceptibles - randomly draw from uniform distribution above
  suscep_matrix[,1] = S0
  infect_matrix[,1] = I0 + intro_mat[,1]
  V1R_matrix[,1] = V1R0
  V1C_matrix[,1] = V1C0
  V2_matrix[,1] = V20
  incid_deaths_matrix[,1] = incid_deaths0
  R_matrix[,1] = N0 - (M0 + S0 +I0 + V1R0+ V1C0+ V20+ V1C0) 
  N_matrix[,1] = N0
  
  
  sd.dw = sd.dw_val
  
  nt0 = rbind(mat_matrix[,1],  suscep_matrix[,1], infect_matrix[,1],
              V1R_matrix[,1], V1C_matrix[,1], V2_matrix[,1], 
              R_matrix[,1], incid_deaths_matrix[,1])
  nt0 = t(nt0)
  nt0 = as.data.frame(nt0)
  nt0$district = rownames(nt0)
  nt0$time = 1
  colnames(nt0) = c("M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
  rownames(nt0) = NULL
  all_out = nt0
  
  
  for (it in 2:(t-1)) {
    
    ### move people around before doing transitions
    # maternal compartment
    M_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      M_move[dd,] = rmultinom(n = 1, size = mat_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    M_move = colSums(M_move)
    mat_matrix[,it-1] = M_move 
    
    # susceptible compartment
    S_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      S_move[dd,] = rmultinom(n = 1, size = suscep_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    S_move = colSums(S_move)
    suscep_matrix[,it-1] = S_move 
    
    # infected compartment
    I_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      I_move[dd,] = rmultinom(n = 1, size = infect_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    I_move = colSums(I_move)
    infect_matrix[,it-1] = I_move 
    
    infect_matrix[, it-1] = infect_matrix[, it-1] + intro_mat[, it-1]
    
    # first vax didn't work
    MR1_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MR1_move[dd,] = rmultinom(n = 1, size = V1R_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MR1_move = colSums(MR1_move)
    V1R_matrix[,it-1] = MR1_move 
    
    
    # first campaign vax that didn't work
    MC1_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MC1_move[dd,] = rmultinom(n = 1, size = V1C_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MC1_move = colSums(MC1_move)
    V1C_matrix[,it-1] = MC1_move   
    
    # second vax
    MR2_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MR2_move[dd,] = rmultinom(n = 1, size = V2_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MR2_move = colSums(MR2_move)
    V2_matrix[,it-1] = MR2_move 
    
    # recovered compartment
    R_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      R_move[dd,] = rmultinom(n = 1, size = R_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    R_move = colSums(R_move)
    R_matrix[,it-1] = R_move 
    
    # update N 
    N_matrix[, it-1] = mat_matrix[, it-1] + suscep_matrix[,it-1] + infect_matrix[,it-1] + V1R_matrix[,it-1] +  V1C_matrix[,it-1] +V2_matrix[,it-1] +  R_matrix[,it-1]
    
    
    dw = rtruncnorm(Ncomp, a =0, mean = 1, sd = sd.dw)
    
    # transitions
    wane_prob = 1 - exp(-wane_mat * delta.t) # probability of moving from maternal compartment to susceptible
    foi_prob <- 1 - exp( - beta_values[season_index[it]]*infect_matrix[, it-1]^alpha/ N_matrix[, it-1] *delta.t * dw) # probability of moving from S to I
    death_prob = 1 - exp(-deaths_rate/1000 / 26 * delta.t)
    death_prob_I = 1 - exp(-deaths_rate/100 / 26 * delta.t - measles_cfr * delta.t)
    MR1_prob = 1 - exp(-mcv1_cov_dat[, it-1] / MR1_time *  delta.t)
    MR2_prob = 1 - exp(-mcv2_cov_dat[, it-1] / MR2_time * delta.t)
    MR1C_prob = 1 - exp(-sia_cov[, it-1]  * delta.t)
    
    # replace any negative values with zero just in case
    wane_prob = replace(wane_prob, wane_prob<0, 0)
    foi_prob = replace(foi_prob, foi_prob<0, 0)
    death_prob = replace(death_prob, death_prob<0, 0)
    death_prob_I = replace(death_prob_I, death_prob_I<0, 0)
    MR1_prob = replace(MR1_prob, MR1_prob<0, 0)
    MR2_prob = replace(MR2_prob, MR2_prob<0, 0)
    MR1C_prob = replace(MR1C_prob, MR1C_prob<0, 0)
    
    ## set up a transition matrix A
    nt = rbind(mat_matrix[,it-1],  suscep_matrix[,it-1], infect_matrix[,it-1],
               V1R_matrix[,it-1], V1C_matrix[,it-1], V2_matrix[,it-1], 
               R_matrix[,it-1], incid_deaths_matrix[,it-1])
    
    # create a list of empty transition matrices for each district
    TM_list = lapply(1:nrow(mat_matrix), matrix, data=0, nrow=dim(nt)[1], ncol=dim(nt)[1])
    for (k in 1:(nrow(mat_matrix))){
      Tm_a = TM_list[[k]]
      
      # transitions out of M compartment
      Tm_a[1, 2] = wane_prob
      Tm_a[1, 8] = death_prob[k,]
      Tm_a[1, 1] = max(0, 1 - wane_prob - death_prob[k,])      
      
      # transitions out of S compartment
      Tm_a[2, 3] = foi_prob[k]
      Tm_a[2, 4] = MR1_prob[k] * (1 - p_eff_V1r)
      Tm_a[2, 5] = MR1C_prob[k] * (1 - p_eff_V1c)
      Tm_a[2, 7] = MR1_prob[k] * p_eff_V1r + MR1C_prob[k] * p_eff_V1c
      Tm_a[2, 8] = death_prob[k,]
      Tm_a[2, 2] = max(0, 1 - Tm_a[2, 3] - Tm_a[2, 4] - Tm_a[2, 5] - Tm_a[2, 7] - Tm_a[2, 8])
      
      
      # transitions out of I compartment
      Tm_a[3, 7] = 1
      Tm_a[3, 8] = death_prob_I[k,]
      Tm_a[3, 3] = max(0, 1 - Tm_a[3,7] - Tm_a[3,8])
      
      # transitions out of V1R compartment
      Tm_a[4, 3] = foi_prob[[k]]
      Tm_a[4, 6] = MR2_prob[k] + MR1C_prob[k] * (p_eff_V1c)
      Tm_a[4, 8] = death_prob[k,]
      Tm_a[4, 4] = max(0, 1 - Tm_a[4, 3] - Tm_a[4, 6] - Tm_a[4, 8])
      
      # transitions out of V1C compartment
      Tm_a[5, 3] = foi_prob[[k]]
      Tm_a[5, 6] = MR2_prob[k]
      Tm_a[5, 8] = death_prob[k,]
      Tm_a[5, 5] = max(0, 1 - Tm_a[5, 3] - Tm_a[5, 6] - Tm_a[5, 8])
      
      # transitions out of V2 compartment
      Tm_a[6, 7] = 1
      Tm_a[6, 8] = death_prob[k,]
      Tm_a[6, 6] = max(0, 1 - Tm_a[6, 7] - Tm_a[6, 8])
      
      # transitions out of R compartment
      Tm_a[7, 8] = death_prob[k,]
      Tm_a[7, 7] = max(0, 1 - Tm_a[7,8])
      
      # transitions out of incid_deaths
      Tm_a[8, 8] = 1
      
      TM_list[[k]] = (Tm_a)  
    }
    
    ### do multinomial draws
    ## set up list to store matrices, one for each district
    TmN_i = list()
    for (kk in 1:nrow(pop_mat)){
      TmN_i[[kk]] = matrix(NA, nrow = nrow(nt), ncol = nrow(nt))
      for (jj in 1:nrow(nt)){
        TmN_i[[kk]][jj,] = rmultinom(1, size = nt[jj,kk], prob = TM_list[[kk]][jj,])
      }
      TmN_i[[kk]] = colSums(TmN_i[[kk]])
    }
    
    out = do.call("rbind", TmN_i)
    out = as.data.frame(out)
    out$district = rownames(pop_mat)
    out$time = it
    colnames(out) = c("M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
    
    # update M compartment with new births
    new_births = rep(NA, nrow(out))
    for (k in 1:nrow(out)){
      new_births[k] = rbinom(1, N_matrix[k, it-1], births_rate$cbr_province[k]/1000/26) 
    }
    new_births = as.data.frame(new_births)
    # update M compartment with new births
    out = out %>% mutate(M = M + new_births$new_births) 
    
    # update N
    N_matrix[, it] = out$M + out$S + out$I + out$V1R + out$V1C + out$V2 + out$R
    
    mat_matrix[,it] = out$M
    suscep_matrix[,it] = out$S 
    infect_matrix[,it] = out$I
    V1R_matrix[,it] = out$V1R
    V1C_matrix[,it] = out$V1C
    V2_matrix[,it] = out$V2
    incid_deaths_matrix[,it] = out$incid_deaths
    R_matrix[,it] = out$R
    
    all_out = rbind(all_out, out)
  }    
  
  return(all_out)
}



## function for running multiple runs
multi_run <- function(num_runs = n_runs, p_mat = pop_mat, ICs = initial_conditions, pars = parameters){
  
  dat <-matrix(nrow = 0, ncol = 11)
  colnames(dat) = c("run_index", "M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
  
  for(n in 1:num_runs){ 
    ## run the simulation one time
    ## update so that it picks the first 
    pars_mob = pars
    pars_mob$mob_data_prob = pars_mob$mob_data_prob[[n]]

    single.sim <- single_run(pop_mat = p_mat, InC = ICs, params = pars_mob)
    
    #add on a value for the run_index (simulation number)
    run_index = rep(n, nrow(single.sim))
    single.sim <- cbind(run_index, single.sim)
    dat <- rbind(single.sim, dat)
    print(n)
  }
  
  return(dat)  
}



plots_post_run <- function(runs_output = run_sc_dd_scale_choma,
                           time_max = 36){
  
  # get total cases
  sum_nat = runs_output %>% group_by(run_index, time) %>% summarize(totalI = sum(I))
  sum_nat = sum_nat %>% ungroup() %>% group_by(time) %>% summarize(meanI = mean(totalI),
                                                                   llI = quantile(totalI, prob= c(0.025)),
                                                                   ulI = quantile(totalI, prob = c(0.975)))
  
  # plot mean number of infections nationwide 
  totalI_plot = ggplot(sum_nat , aes(x=time, y=meanI)) +
    geom_line() +
    xlab("Time") + 
    ylab("Measles infections") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 8),strip.text=element_text(size=6))  +
    #  scale_x_continuous(breaks = seq(79,234,26), 
    #                       labels = c("79" = "2019", "105" = "2020", "131" = "2021", "157" = "2022",
    #                                 "183" = "2023","209" = "2024")) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
  
  # look at introduction from 2020 onwards (until 2022)
  prop_intro = runs_output %>% ungroup() %>% filter(time<=time_max) %>% group_by(district, time) %>% summarize(num_nonzeroI = sum(I > 0))
  prop_intro = prop_intro %>% mutate(prop_nonzeroI = num_nonzeroI / n_sims)
  prop_intro_nonzero = prop_intro  %>%  ungroup() %>% group_by(district) %>% mutate(any_intros = sum(prop_nonzeroI>0)) %>% filter(any_intros !=0)

  
  # make a heatmap of introductions
  heatmap_p = ggplot(prop_intro_nonzero %>% filter(time <=time_max), aes(x = time, y = district, fill = prop_nonzeroI)) +
    geom_tile()+
    scale_x_continuous(breaks = c(2, 11, 19, 28, 37), 
                       labels = c("2" = "Sept 2020", "11" = "Jan 2021", "19" = "Apr 2021", "28" = "Sept 2021",  "37" = "Jan 2022"))+
    guides(fill=guide_legend(title="Proportion scenarios with at least 1 incident case"))
  
  
  
  # Check when introduction first happens for each district
  
  
  
  # look at when outbreaks are introduced
  # add cumulative infection
  runs_output = runs_output %>% group_by(run_index, district) %>% mutate(cumI = cumsum(I))
  totalI_index = runs_output %>% ungroup() %>% group_by(district, run_index) %>% arrange(time) %>% dplyr::slice(n())
  # make a table with % of scenarios that have case for each district
  totalI_index = totalI_index %>% ungroup() %>% group_by(district) %>% summarize(prop_non0 = sum(cumI>0)/length(unique(run_index)))
  # keep only nonzero ones
  totalI_index = totalI_index %>% filter(prop_non0>0) %>% arrange(desc(prop_non0))
  totalI_index = totalI_index %>% mutate(district_tit = str_to_title(district))
  totalI_index_dd_scale_choma = totalI_index
  #totalI_index %>% flextable()
  # restrict to those with at least 10% chance of introduction
  prop_w_intro_plot = ggplot(totalI_index_dd_scale_choma %>% filter(prop_non0>0.1)) +
    geom_bar( aes(x=(reorder(district_tit, desc(prop_non0))), y=prop_non0), stat="identity", fill="#F39B7FB2", alpha=0.7) +
    theme_bw() + xlab("District") + theme(text = element_text(size = 14)) +ylab("Proportion simulations with non-zero cases") + ggtitle("Proportion of simulations with at least one measles case")+
    theme(axis.text.x = element_text(angle = 45))
  
  outbreaks_only = runs_output %>% ungroup() %>% filter(cumI>=1)
  # filter to those with at least 10% probability of outbreak
  outbreaks_only = outbreaks_only %>% filter(district %in% prop_w_intro_plot$data$district)
  # only keep the first date in each simulation
  outbreaks_only = outbreaks_only %>% group_by(run_index, district) %>% arrange(time) %>% slice(1)
  
  # median time to outbreak
  outbreaks_only = outbreaks_only %>% ungroup() %>% group_by(district) %>% summarize(median_time = median(time),
                                                                                     ll_time = quantile(time, probs = c(0.025)), 
                                                                                     ul_time = quantile(time, probs = c(0.975)))
  
  # filter out - those that have median time > 36 weeks
  outbreaks_only = outbreaks_only %>% filter(median_time<=time_max)
  
  outbreaks_only_bar = ggplot(outbreaks_only) +
    geom_bar( aes(x=reorder(district, median_time), y=median_time), stat="identity", fill="#00A087B2", alpha=0.7) +
    geom_pointrange( aes(x=reorder(district, median_time), y=median_time, ymin=ll_time, ymax = ul_time), colour="#8491B4B2", alpha=0.9, size=1)+
    theme_bw() + xlab("District") + theme(text = element_text(size = 14)) +ylab("Median time (2-weeks)") + ggtitle("Median time to outbreak")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  
  # make a map
  # check how many districts are in shp_zam
  outbreaks_only = outbreaks_only %>% mutate(in_map = ifelse(district %in% shp_zam$district_lower, 1, 0))
  outbreaks_only = outbreaks_only %>% dplyr::select(-c(in_map)) # all in map alhumdulilah
  
  # scale time by 2
  outbreaks_only$median_weeks = outbreaks_only$median_time * 2
  outbreaks_only$ll_weeks = outbreaks_only$ll_time * 2
  outbreaks_only$ul_weeks = outbreaks_only$ul_time * 2
  
  outbreak_time = merge(shp_zam, outbreaks_only, by.x = "district_lower", by.y = "district", all.x = TRUE)
  outbreak_time@data$id = rownames(outbreak_time@data)
  new_zm <- fortify(outbreak_time, region = "id")
  newdf_zm <- merge(new_zm, outbreak_time@data, by = "id")
  
  # centroids for labels
  centroids_nonNA = centroids_SAdf %>% mutate(lower_district = tolower(ADM2))
  centroids_nonNA = centroids_nonNA %>% mutate(intros = ifelse(lower_district %in% outbreaks_only$district, 1, 0))
  centroids_nonNA = centroids_nonNA %>% filter(intros==1)
  
  outbreak_time_map <- ggplot() +
    #geom_line(aes(x = long, y = lat), color = "gray")+
    geom_polygon(data = newdf_zm, aes(fill = median_weeks,
                                      x = long, 
                                      y = lat, 
                                      group = group),
                 color = "grey",size = 0.2) + 
    #lims(fill = c(0, 20))+
    #geom_text_repel(data = centroids_nonNA, aes(x, y, label = ADM2))+
    #scale_fill_material("pink", reverse = TRUE, na.value="grey95")+
    scale_fill_viridis_c(option = "mako", direction = 1, na.value = "grey95")+
    #geom_point(data = ndola_point, aes(x = x, y = y), color = "red") +
    #geom_point(data = choma_point, aes(x = x, y = y), color = "red") +
    labs(fill = "Weeks to outbreak")+ #scale_y_continuous(trans='log2')+
    coord_map() + ggtitle("Median time to outbreak") + theme(plot.title = element_text(hjust =0.5))+ theme_void()#
  

  # epi curve for each district with an outbreak
  # limit to first 36 weeks (ending on 2022)
  
  ts_36wk = runs_output %>% ungroup() %>% filter(time <=time_max) %>% group_by(district, time) %>% summarize(meanI = mean(I), medianI = median(I),
                                                                                                                 llI = quantile(I, probs = c(0.025)),
                                                                                                                 ulI = quantile(I, probs = c(0.975)),
                                                                                                                 median_cumI = median(cumI),
                                                                                                                 mean_cumI = mean(cumI),
                                                                                                                 ll_cumI = quantile(cumI, probs = c(0.025)),
                                                                                                                 ul_cumI = quantile(cumI, probs = c(0.975)))
  # epi curve for each district where there was at least 1 case (median)
  ts_36wk_nonzer = ts_36wk %>% ungroup() %>% group_by(district) %>% filter(max(meanI)>=1) # districts
  
  epi_curves = ggplot(ts_36wk_nonzer)+
    geom_line(aes(x = time, y = meanI, color = district), show.legend = FALSE) + 
    geom_ribbon(aes(x = time, ymin = llI, ymax = ulI, fill = district), alpha = 0.4, show.legend = FALSE)+
    theme_minimal()+ ylab("Number of infections")+
    scale_x_continuous(breaks = c(2, 11, 19, 28, 37), 
                       labels = c("2" = "Sept 2020", "11" = "Jan 2021", "19" = "Apr 2021", "28" = "Sept 2021",  "37" = "Jan 2022")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~district, scales = "free_y")
  
  out = list(totalI_plot, heatmap_p, prop_w_intro_plot, outbreaks_only_bar, outbreak_time_map, epi_curves, totalI_index, outbreaks_only)
  
}



### create single_run version for district outbreak response

single_run_outbreak <- function(pop_mat = population_matrix, 
                       InC = initial_conditions,
                       params = parameters
) {
  M0 = InC$M
  lower_S0 = InC$S.lower_S0
  upper_S0 = InC$S.upper_S0
  mean_S0 = InC$S.mean_S0
  S0= (1 - runif(length(lower_S0), min = lower_S0, max = upper_S0))*InC$N0
  I0 = InC$I0
  V1R0 = InC$V1R0
  V1C0 = InC$V1C0
  V20 = InC$V20
  incid_deaths0 = InC$incid_deaths0
  # R0 = InC$R0  ## will need to update after a random draw
  N0 = InC$N0
  
  sd.dw_val = params$sd.dw.val
  t = params$t
  mob_data_prob = params$mob_data_prob
  mcv1_cov_dat = params$mcv1_cov_dat
  mcv2_cov_dat = params$mcv2_cov_dat
  sia_cov = params$sia_cov
  births_rate = params$births_rate
  wane_mat = params$wane_mat
  beta_values = params$beta_values
  alpha = params$alpha
  deaths_rate = params$deaths_rate
  measles_cfr = params$measles_cfr
  MR1_time = params$MR1_time
  MR2_time = params$MR2_time
  p_eff_V1r = params$p_eff_V1r
  p_eff_V1c = params$p_eff_V1c
  intro_mat = params$intro_mat
  
  mat_matrix = matrix(NA,dim(pop_mat)[1],t)
  suscep_matrix <- matrix(NA,dim(pop_mat)[1],t)
  infect_matrix <- matrix(NA,dim(pop_mat)[1],t)
  V1R_matrix = matrix(NA,dim(pop_mat)[1],t)
  V1C_matrix = matrix(NA,dim(pop_mat)[1],t)
  V2_matrix = matrix(NA,dim(pop_mat)[1],t)
  R_matrix = matrix(NA,dim(pop_mat)[1],t)
  incid_deaths_matrix = matrix(NA,dim(pop_mat)[1],t)
  N_matrix = matrix(NA, dim(pop_mat)[1],t)
  
  rownames(mat_matrix) <- rownames(pop_mat)
  rownames(suscep_matrix) <- rownames(pop_mat)
  rownames(infect_matrix) <- rownames(pop_mat)
  rownames(V1R_matrix) <- rownames(pop_mat)
  rownames(V1C_matrix) <- rownames(pop_mat)
  rownames(V2_matrix) <- rownames(pop_mat)
  rownames(incid_deaths_matrix) <- rownames(pop_mat)
  rownames(R_matrix) <- rownames(pop_mat)
  rownames(N_matrix) = rownames(pop_mat)
  
  # separate matrix for SIA triggered by outbreak
  sia_outbreak = matrix(0, nrow = nrow(sia_cov), ncol = ncol(sia_cov))
  rownames(sia_outbreak) = rownames(sia_cov)
  
  mat_matrix[,1] = M0
  # for susceptibles - randomly draw from uniform distribution above
  suscep_matrix[,1] = S0
  infect_matrix[,1] = I0 + intro_mat[,1]
  V1R_matrix[,1] = V1R0
  V1C_matrix[,1] = V1C0
  V2_matrix[,1] = V20
  incid_deaths_matrix[,1] = incid_deaths0
  R_matrix[,1] = N0 - (M0 + S0 +I0 + V1R0+ V1C0+ V20+ V1C0) 
  N_matrix[,1] = N0
  

  
  sd.dw = sd.dw_val
  
  nt0 = rbind(mat_matrix[,1],  suscep_matrix[,1], infect_matrix[,1],
              V1R_matrix[,1], V1C_matrix[,1], V2_matrix[,1], 
              R_matrix[,1], incid_deaths_matrix[,1])
  nt0 = t(nt0)
  nt0 = as.data.frame(nt0)
  nt0$district = rownames(nt0)
  nt0$time = 1
  colnames(nt0) = c("M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
  rownames(nt0) = NULL
  all_out = nt0
  
  
  for (it in 2:(t-1)) {
    
    ### move people around before doing transitions
    # maternal compartment
    M_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      M_move[dd,] = rmultinom(n = 1, size = mat_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    M_move = colSums(M_move)
    mat_matrix[,it-1] = M_move 
    
    # susceptible compartment
    S_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      S_move[dd,] = rmultinom(n = 1, size = suscep_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    S_move = colSums(S_move)
    suscep_matrix[,it-1] = S_move 
    
    # infected compartment
    I_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      I_move[dd,] = rmultinom(n = 1, size = infect_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    I_move = colSums(I_move)
    infect_matrix[,it-1] = I_move 
    
    infect_matrix[, it-1] = infect_matrix[, it-1] + intro_mat[, it-1]
    
    # first vax didn't work
    MR1_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MR1_move[dd,] = rmultinom(n = 1, size = V1R_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MR1_move = colSums(MR1_move)
    V1R_matrix[,it-1] = MR1_move 
    
    
    # first campaign vax that didn't work
    MC1_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MC1_move[dd,] = rmultinom(n = 1, size = V1C_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MC1_move = colSums(MC1_move)
    V1C_matrix[,it-1] = MC1_move   
    
    # second vax
    MR2_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      MR2_move[dd,] = rmultinom(n = 1, size = V2_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    MR2_move = colSums(MR2_move)
    V2_matrix[,it-1] = MR2_move 
    
    # recovered compartment
    R_move = matrix(NA, nrow = nrow(mob_data_prob), ncol = ncol(mob_data_prob))
    for (dd in 1:nrow(mat_matrix)){
      R_move[dd,] = rmultinom(n = 1, size = R_matrix[dd, it-1], prob = mob_data_prob[dd,])
    }
    # add up columns to get new number of people in each district
    R_move = colSums(R_move)
    R_matrix[,it-1] = R_move 
    
    # update N 
    N_matrix[, it-1] = mat_matrix[, it-1] + suscep_matrix[,it-1] + infect_matrix[,it-1] + V1R_matrix[,it-1] +  V1C_matrix[,it-1] +V2_matrix[,it-1] +  R_matrix[,it-1]
    
    
    dw = rtruncnorm(Ncomp, a =0, mean = 1, sd = sd.dw)
    
    # transitions
    wane_prob = 1 - exp(-wane_mat * delta.t) # probability of moving from maternal compartment to susceptible
    foi_prob <- 1 - exp( - beta_values[season_index[it]]*infect_matrix[, it-1]^alpha/ N_matrix[, it-1] *delta.t * dw) # probability of moving from S to I
    death_prob = 1 - exp(-deaths_rate/1000 / 26 * delta.t)
    death_prob_I = 1 - exp(-deaths_rate/100 / 26 * delta.t + measles_cfr * delta.t)
    MR1_prob = 1 - exp(-mcv1_cov_dat[, it-1] / MR1_time *  delta.t)
    MR2_prob = 1 - exp(-mcv2_cov_dat[, it-1] / MR2_time * delta.t)
    MR1C_prob = 1 - exp(-(sia_cov[, it-1]+sia_outbreak[, it-1])  * delta.t)
    
    # replace any negative values with zero just in case
    wane_prob = replace(wane_prob, wane_prob<0, 0)
    foi_prob = replace(foi_prob, foi_prob<0, 0)
    death_prob = replace(death_prob, death_prob<0, 0)
    death_prob_I = replace(death_prob_I, death_prob_I<0, 0)
    MR1_prob = replace(MR1_prob, MR1_prob<0, 0)
    MR2_prob = replace(MR2_prob, MR2_prob<0, 0)
    MR1C_prob = replace(MR1C_prob, MR1C_prob<0, 0)
    
    ## set up a transition matrix A
    nt = rbind(mat_matrix[,it-1],  suscep_matrix[,it-1], infect_matrix[,it-1],
               V1R_matrix[,it-1], V1C_matrix[,it-1], V2_matrix[,it-1], 
               R_matrix[,it-1], incid_deaths_matrix[,it-1])
    
    # create a list of empty transition matrices for each district
    TM_list = lapply(1:nrow(mat_matrix), matrix, data=0, nrow=dim(nt)[1], ncol=dim(nt)[1])
    for (k in 1:(nrow(mat_matrix))){
      Tm_a = TM_list[[k]]
      
      # transitions out of M compartment
      Tm_a[1, 2] = wane_prob
      Tm_a[1, 8] = death_prob[k,]
      Tm_a[1, 1] = max(0, 1 - wane_prob - death_prob[k,])      
      
      # transitions out of S compartment
      Tm_a[2, 3] = foi_prob[k]
      Tm_a[2, 4] = MR1_prob[k] * (1 - p_eff_V1r)
      Tm_a[2, 5] = MR1C_prob[k] * (1 - p_eff_V1c)
      Tm_a[2, 7] = MR1_prob[k] * p_eff_V1r + MR1C_prob[k] * p_eff_V1c
      Tm_a[2, 8] = death_prob[k,]
      Tm_a[2, 2] = max(0, 1 - Tm_a[2, 3] - Tm_a[2, 4] - Tm_a[2, 5] - Tm_a[2, 7] - Tm_a[2, 8])
      
      
      # transitions out of I compartment
      Tm_a[3, 7] = 1
      Tm_a[3, 8] = death_prob_I[k,]
      Tm_a[3, 3] = max(0, 1 - Tm_a[3,7] - Tm_a[3,8])
      
      # transitions out of V1R compartment
      Tm_a[4, 3] = foi_prob[[k]]
      Tm_a[4, 6] = MR2_prob[k] + MR1C_prob[k] * (p_eff_V1c)
      Tm_a[4, 8] = death_prob[k,]
      Tm_a[4, 4] = max(0, 1 - Tm_a[4, 3] - Tm_a[4, 6] - Tm_a[4, 8])
      
      # transitions out of V1C compartment
      Tm_a[5, 3] = foi_prob[[k]]
      Tm_a[5, 6] = MR2_prob[k]
      Tm_a[5, 8] = death_prob[k,]
      Tm_a[5, 5] = max(0, 1 - Tm_a[5, 3] - Tm_a[5, 6] - Tm_a[5, 8])
      
      # transitions out of V2 compartment
      Tm_a[6, 7] = 1
      Tm_a[6, 8] = death_prob[k,]
      Tm_a[6, 6] = max(0, 1 - Tm_a[6, 7] - Tm_a[6, 8])
      
      # transitions out of R compartment
      Tm_a[7, 8] = death_prob[k,]
      Tm_a[7, 7] = max(0, 1 - Tm_a[7,8])
      
      # transitions out of incid_deaths
      Tm_a[8, 8] = 1
      
      TM_list[[k]] = (Tm_a)  
    }
    
    ### do multinomial draws
    ## set up list to store matrices, one for each district
    TmN_i = list()
    for (kk in 1:nrow(pop_mat)){
      TmN_i[[kk]] = matrix(NA, nrow = nrow(nt), ncol = nrow(nt))
      for (jj in 1:nrow(nt)){
        TmN_i[[kk]][jj,] = rmultinom(1, size = nt[jj,kk], prob = TM_list[[kk]][jj,])
      }
      TmN_i[[kk]] = colSums(TmN_i[[kk]])
    }
    
    out = do.call("rbind", TmN_i)
    out = as.data.frame(out)
    out$district = rownames(pop_mat)
    out$time = it
    colnames(out) = c("M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
    
    # update M compartment with new births
    new_births = rep(NA, nrow(out))
    for (k in 1:nrow(out)){
      new_births[k] = rbinom(1, N_matrix[k, it-1], births_rate$cbr_province[k]/1000/26) 
    }
    new_births = as.data.frame(new_births)
    # update M compartment with new births
    out = out %>% mutate(M = M + new_births$new_births) 
    
    # update N
    N_matrix[, it] = out$M + out$S + out$I + out$V1R + out$V1C + out$V2 + out$R
    
    mat_matrix[,it] = out$M
    suscep_matrix[,it] = out$S 
    infect_matrix[,it] = out$I
    V1R_matrix[,it] = out$V1R
    V1C_matrix[,it] = out$V1C
    V2_matrix[,it] = out$V2
    incid_deaths_matrix[,it] = out$incid_deaths
    R_matrix[,it] = out$R
    
    all_out = rbind(all_out, out)
    
    # calculate cumulative infections by district
    cum_inf = as.data.frame(all_out) %>% ungroup() %>% group_by(district) %>% dplyr::mutate(cuminf = cumsum(I)) 
    cum_inf = cum_inf %>% filter(cuminf>=3) %>% arrange(time)
    cum_inf = cum_inf %>% slice(1)
    
    # update sia_outbreak ---> assign sia to district of outbreak 2 time steps after first time they make it in
    #dist_w_outbreaks = cum_inf$district
    for (dout in 1:nrow(cum_inf)){
      sia_outbreak[cum_inf$district[dout], cum_inf$time[dout]+2] = 0.8
    }
  }    
  
  return(all_out)
}


# including district-level outbreak response

multi_run_outbreak <- function(num_runs = n_runs, p_mat = pop_mat, ICs = initial_conditions, pars = parameters){
  
  dat <-matrix(nrow = 0, ncol = 11)
  colnames(dat) = c("run_index", "M", "S", "I", "V1R", "V1C", "V2", "R", "incid_deaths", "district", "time")
  
  for(n in 1:num_runs){ 
    ## run the simulation one time
    single.sim <- single_run_outbreak(pop_mat = p_mat, InC = ICs, params = pars)
    
    #add on a value for the run_index (simulation number)
    run_index = rep(n, nrow(single.sim))
    single.sim <- cbind(run_index, single.sim)
    dat <- rbind(single.sim, dat)
    print(n)
  }
  
  return(dat)  
}
