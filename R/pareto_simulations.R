# Pareto simulations to determine the lower bound (location parameter)
# required to achieve different levels of viral load suppression (VLS),
# defined as P(VL ≤ 1000 copies/mL).

# For a given shape and target VLS proportion, the function:
# 1) Solves for the Pareto location parameter such that
#    ppareto(log10(1000)) = target VLS.
# 2) Converts the location (log10 scale) back to the original scale
#    to obtain the lower bound of the VL distribution.
# The function is vectorized to evaluate multiple VLS targets efficiently.

# The script evaluates VLS levels from 10% to 99% across several
# shape parameters (including 1.73 from Johnson et al. and alternatives),
# reshapes results into tidy format, and combines them.

# Finally, it plots the relationship between:
# - Target suppression proportion (x-axis)
# - Required lower bound of the VL distribution (y-axis)
# across different Pareto shape assumptions.
pareto_vec_func <- Vectorize(pareto_func <- function(shape, vls_est){
  
  myfun <- function(location) {
    shape <- shape
    EnvStats::ppareto(log10(1000), shape=shape, location) - vls_est
  }
  
  tmp <- uniroot(myfun, lower=0.01, upper=100)
  
  mylocation <- tmp$root
  
  lower_limit <- 10^mylocation
  return(data.frame(scale = mylocation, lower_limit = lower_limit))
  
})


pareto_vec_func(shape = 1.73, vls_est = .5)

# create data frame of values to evaluate
# Johnson et al
vls_df_shape_johnson <- data.frame(vls_est = seq(.1,.99,.01), shape = rep(1.73))

lower_lim_est_johnson <- pareto_vec_func(shape = vls_df_shape_johnson$shape, vls_est = vls_df_shape_johnson$vls_est)

lower_lim_est_johnson <- as.data.frame(lower_lim_est_johnson)

lower_lim_est_johnson <- setNames(cbind(rownames(lower_lim_est_johnson), lower_lim_est_johnson, row.names = NULL), 
                                  c("indicator", seq(.1,.99,.01)))

lower_lim_johnson <- lower_lim_est_johnson %>%
  tidyr::pivot_longer(cols = 2:91)

scale_johnson <- lower_lim_johnson %>% filter(indicator == "scale") %>%
  select(name, value) %>%
  rename(perc_supp = name, scale_val = value) 

lower_lim_johnson <- lower_lim_johnson %>% filter(indicator == "lower_limit") %>%
  select(name, value) %>%
  rename(perc_supp = name, lower_limit_val = value) 

scale_lower_johnson <- left_join(scale_johnson, lower_lim_johnson, by = "perc_supp")
scale_lower_johnson <- scale_lower_johnson %>%
  mutate(perc_supp = as.numeric(perc_supp), scale_val = as.numeric(scale_val), lower_limit_val = as.numeric(lower_limit_val),
         shape_val = rep("1.73 (Johnson et al.)"))

# shape = 1.20
vls_df_shape_one_two <- data.frame(vls_est = seq(.1,.99,.01), shape = rep(1.2))

lower_lim_est_one_two <- pareto_vec_func(shape = vls_df_shape_one_two$shape, vls_est = vls_df_shape_one_two$vls_est)

lower_lim_est_one_two <- as.data.frame(lower_lim_est_one_two)

lower_lim_est_one_two <- setNames(cbind(rownames(lower_lim_est_one_two), lower_lim_est_one_two, row.names = NULL), 
                              c("indicator", seq(.1,.99,.01)))

lower_lim_one_two <- lower_lim_est_one_two %>%
  tidyr::pivot_longer(cols = 2:91)

scale_one_two <- lower_lim_one_two %>% filter(indicator == "scale") %>%
  select(name, value) %>%
  rename(perc_supp = name, scale_val = value) 

lower_lim_one_two <- lower_lim_one_two %>% filter(indicator == "lower_limit") %>%
  select(name, value) %>%
  rename(perc_supp = name, lower_limit_val = value) 

scale_lower_one_two <- left_join(scale_one_two, lower_lim_one_two, by = "perc_supp")
scale_lower_one_two <- scale_lower_one_two %>%
  mutate(perc_supp = as.numeric(perc_supp), scale_val = as.numeric(scale_val), lower_limit_val = as.numeric(lower_limit_val),
         shape_val = rep("1.20"))

# shape = 2.00
vls_df_shape_two <- data.frame(vls_est = seq(.1,.99,.01), shape = rep(2))

lower_lim_est_two <- pareto_vec_func(shape = vls_df_shape_two$shape, vls_est = vls_df_shape_two$vls_est)

lower_lim_est_two <- as.data.frame(lower_lim_est_two)

lower_lim_est_two <- setNames(cbind(rownames(lower_lim_est_two), lower_lim_est_two, row.names = NULL), 
                              c("indicator", seq(.1,.99,.01)))

lower_lim_two <- lower_lim_est_two %>%
  tidyr::pivot_longer(cols = 2:91)

scale_two <- lower_lim_two %>% filter(indicator == "scale") %>%
  select(name, value) %>%
  rename(perc_supp = name, scale_val = value) 

lower_lim_two <- lower_lim_two %>% filter(indicator == "lower_limit") %>%
  select(name, value) %>%
  rename(perc_supp = name, lower_limit_val = value) 

scale_lower_two <- left_join(scale_two, lower_lim_two, by = "perc_supp")
scale_lower_two <- scale_lower_two %>%
  mutate(perc_supp = as.numeric(perc_supp), scale_val = as.numeric(scale_val), lower_limit_val = as.numeric(lower_limit_val),
         shape_val = rep("2.00"))


# shape = 2.5
vls_df_shape_two_five <- data.frame(vls_est = seq(.1,.99,.01), shape = rep(2.5))

lower_lim_est_two_five <- pareto_vec_func(shape = vls_df_shape_two_five$shape, vls_est = vls_df_shape_two_five$vls_est)

lower_lim_est_two_five <- as.data.frame(lower_lim_est_two_five)

lower_lim_est_two_five <- setNames(cbind(rownames(lower_lim_est_two_five), lower_lim_est_two_five, row.names = NULL), 
                                  c("indicator", seq(.1,.99,.01)))

lower_lim_two_five <- lower_lim_est_two_five %>%
  tidyr::pivot_longer(cols = 2:91)

scale_two_five <- lower_lim_two_five %>% filter(indicator == "scale") %>%
  select(name, value) %>%
  rename(perc_supp = name, scale_val = value) 

lower_lim_two_five <- lower_lim_two_five %>% filter(indicator == "lower_limit") %>%
  select(name, value) %>%
  rename(perc_supp = name, lower_limit_val = value) 

scale_lower_two_five <- left_join(scale_two_five, lower_lim_two_five, by = "perc_supp")
scale_lower_two_five <- scale_lower_two_five %>%
  mutate(perc_supp = as.numeric(perc_supp), scale_val = as.numeric(scale_val), lower_limit_val = as.numeric(lower_limit_val),
         shape_val = rep("2.50"))

# shape = 3
vls_df_shape_three <- data.frame(vls_est = seq(.1,.99,.01), shape = rep(3))

lower_lim_est_three <- pareto_vec_func(shape = vls_df_shape_three$shape, vls_est = vls_df_shape_three$vls_est)

lower_lim_est_three <- as.data.frame(lower_lim_est_three)

lower_lim_est_three <- setNames(cbind(rownames(lower_lim_est_three), lower_lim_est_three, row.names = NULL), 
                                c("indicator", seq(.1,.99,.01)))

lower_lim_three <- lower_lim_est_three %>%
  tidyr::pivot_longer(cols = 2:91)

scale_three <- lower_lim_three %>% filter(indicator == "scale") %>%
  select(name, value) %>%
  rename(perc_supp = name, scale_val = value) 

lower_lim_three <- lower_lim_three %>% filter(indicator == "lower_limit") %>%
  select(name, value) %>%
  rename(perc_supp = name, lower_limit_val = value) 

scale_lower_three <- left_join(scale_three, lower_lim_three, by = "perc_supp")
scale_lower_three <- scale_lower_three %>%
  mutate(perc_supp = as.numeric(perc_supp), scale_val = as.numeric(scale_val), lower_limit_val = as.numeric(lower_limit_val),
         shape_val = rep("3.00"))


scale_lower_all <- bind_rows(scale_lower_johnson, scale_lower_one_two,
                             scale_lower_two, scale_lower_two_five)#, scale_lower_three)
# Figure S4 code
scale_lower_all %>%
  ggplot(aes(x = perc_supp, y = lower_limit_val, color = factor(shape_val))) + 
  geom_line(linewidth = 0.8) +
  scale_x_continuous(labels = scales::percent, breaks = seq(.1,1,.1)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,200,20)) +
  #scale_color_manual(values = c("green", "red", "blue", "black")) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept =  20, linetype = "dashed", color = "grey") +
  coord_cartesian(xlim = c(.4,1), ylim = c(0,100)) +
  theme_bw(base_size = 11) +
  labs(y = "Lower bound of VL values (copies/mL)", 
       x = "Percentage with VL \u22641000 copies/mL",
       color = "shape parameter") +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1), family = "sans"),
        #panel.grid = element_blank(),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.0), family = "sans"))

# Figure S2 code
# plot distributions using varying scales for illustrative countries specific countries alongside histograms of observed VL data
civ <- ssa_surveys_art %>% filter(country == "Côte d'Ivoire (2017-18)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.9871715), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.387628), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 1.539887), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 1.759592), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 1.923219), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 4)) +
  labs(x = "", color = "shape parameter", y = "Density",title = "A. Côte d'Ivoire (2017-18): 73.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.0), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
civ

nig <- ssa_surveys_art %>% filter(country == "Nigeria (2018)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.7558794), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.153083), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 1.311959), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 1.54796), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 1.728438), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 4.9)) +
  labs(x = "", color = "", y = "Density", title = "B. Nigeria (2018): 80.9% VLS") +
  #title = "A. Nigeria (2018): 80.9% \U2264 1000 copies/mL") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nig

les<-ssa_surveys_art %>% filter(country == "Lesotho (2016-17)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.5238001), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.8940994), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 1.052837), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 1.298115), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 1.49264), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 6.4)) +
  #geom_vline(xintercept = log10(200), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  #geom_vline(xintercept = log10(1000), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  labs(x = "", color = "shape parameter", y = "Density", title = "C. Lesotho (2016-2017): 87.7% VLS") +
  #title = "B. Lesotho (2016-2017)") +
  #title = "B. Lesotho (2016-2017): 87.7% \U2264 1000 copies/mL") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        legend.title = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

les

nam<-ssa_surveys_art %>% filter(country == "Namibia (2017)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.391511), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.7306145), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 0.8841047), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 1.128806), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 1.32854), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 6.2)) +
  #geom_vline(xintercept = log10(200), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  #geom_vline(xintercept = log10(1000), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title ="D. Namibia (2017): 91.3% VLS") +
  #title = "C. Namibia (2017): 91.3% \U2264 1000 copies/mL") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nam

esw<-ssa_surveys_art %>% filter(country == "Eswatini (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.1966564), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.4531612), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 0.5848739), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 0.8110876), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 1.008686), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 7.2)) +
  #geom_vline(xintercept = log10(200), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  #geom_vline(xintercept = log10(1000), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title = "E. Eswatini (2021): 96.2% VLS") +
  #title = "B. Lesotho (2016-2017)") +
  #title = "B. Lesotho (2016-2017): 87.7% \U2264 1000 copies/mL") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
esw

bots<-ssa_surveys_art %>% filter(country == "Botswana (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.2, location = 0.1196038), 
                aes(colour = "1.20"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.3209482), 
                aes(colour = "1.73 (Johnson et al.)"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.0, location = 0.4340284), 
                aes(colour = "2.00"), linewidth = 1) +
  stat_function(fun = EnvStats::dpareto, args = list(shape = 2.5, location = 0.6388788), 
                aes(colour = "2.50"), linewidth = 1) +
  #stat_function(fun = EnvStats::dpareto, args = list(shape = 3, location = 0.8267403), 
  #              aes(colour = "3.00"), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(1.5, 7)) +
  #geom_vline(xintercept = log10(200), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  #geom_vline(xintercept = log10(1000), col = "darkgrey", linewidth = 0.5, linetype = "dashed") +
  labs(x = "Viral load count (copies/mL)", color = "shape parameter", y = "Density", title = "F. Botswana (2021): 97.9% VLS") +
  #title = "D. Botswana (2021): 97.9% \U2264 1000 copies/mL") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        legend.title = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

bots

cowplot::plot_grid(civ,nig,les,nam,esw,bots, rel_widths = c(.8,.8,1.1))
