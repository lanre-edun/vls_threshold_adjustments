# source script that reads data and functions required for analyses
source(file = "C:/Users/oedun/OneDrive - Imperial College London/R projects/viremia_thresholds/R/data_preparation.R")
#' select analytic variables and merge all surveys for plotting

table(ssa_surveys$country, useNA = "always")
table(ssa_surveys$country, ssa_surveys$viral_sup1000, useNA = "always")
table(ssa_surveys$country, ssa_surveys$viral_sup400, useNA = "always")
table(ssa_surveys$country, ssa_surveys$viral_sup200, useNA = "always")
table(ssa_surveys$country, ssa_surveys$viral_sup50, useNA = "always")
table(ssa_surveys_art$resultvlc_numeric, useNA = "always")

# Code for Figure 2
# using new and old parameters
# plot distributions for selected illustrative countries
civ <- ssa_surveys_art %>% filter(country == "Côte d'Ivoire (2017-18)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) +  
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 2.137698),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 4.573402)},
                aes(colour = "Reverse Weibull"), linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.387628), 
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 2.185994),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 4.464708)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.9871715), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  labs(x = "", color = "", y = "Density",title = "A. Côte d'Ivoire (2017-18): 73.7% VLS") +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 4.1)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
civ

nig <- ssa_surveys_art %>% filter(country == "Nigeria (2018)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.65945),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 5.207989)},
                aes(colour = "Reverse Weibull"),linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.153083),
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.725519),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 5.046664)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.7558794), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 5.1)) +
  labs(x = "", color = "", y = "Density", title = "B. Nigeria (2018): 80.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nig

les <- ssa_surveys_art %>% filter(country == "Lesotho (2016-17)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.257307),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 6.176566)},
                aes(colour = "Reverse Weibull"),
                linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.8940994), 
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.33149),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 5.927286)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.5238001), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 6.6)) +
  labs(x = "", color = "", y = "Density", title = "C. Lesotho (2016-2017): 87.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
les

nam <- ssa_surveys_art %>% filter(country == "Namibia (2017)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.048608),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 7.043988)},
                aes(colour = "Reverse Weibull"),
                linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.7306145), 
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.123896),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 6.709182)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.391511), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 6.3)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title ="D. Namibia (2017): 91.3% VLS") +
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
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.74435),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 9.53912)},
                aes(colour = "Reverse Weibull"),
                linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.4531612), 
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 0.8160148),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 8.929944)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.1966564), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 7.4)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title = "E. Eswatini (2021): 96.2% VLS") +
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
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.6111164),
                aes(colour = "Weibull"), linewidth = .8) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 11.83286)},
                aes(colour = "Reverse Weibull"),
                linewidth = .8) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.3209482), 
                aes(colour = "Pareto"), linewidth = .8) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 0.6787366),
                aes(colour = "Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 10.94185)},
                aes(colour = "Reverse Weibull"), linewidth = .8, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.1196038), 
                aes(colour = "Pareto"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#4292C6", "orange","red"))+ 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6)) +
  ggbreak::scale_y_break(c(.8, 7.1)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title = "F. Botswana (2021): 97.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
bots

cowplot::plot_grid(civ,nig,les,nam,esw,bots, rel_widths = c(.7,.7,1))
cowplot::plot_grid(civ,nig,les,nam,esw,bots)

# code for Figure S6
# gamma, frechet and lognormal
# plot distributions for all included countries
civ <- ssa_surveys_art %>% filter(country == "Côte d'Ivoire (2017-18)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 2.137698),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 4.573402)},
                aes(colour = "Reverse Weibull"), linewidth = 1,linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.387628), 
                aes(colour = "Pareto"), linewidth = 1,linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = 0.5355143), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 1.586612), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 2.79633), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red"))+ 
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "", color = "", y = "Density",title = "A. Côte d'Ivoire (2017-18): 73.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
civ

nig <- ssa_surveys_art %>% filter(country == "Nigeria (2018)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.65945),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 5.207989)},
                aes(colour = "Reverse Weibull"),linewidth = 1, linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.153083),
                aes(colour = "Pareto"), linewidth = 1, linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = 0.3213722), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 1.303844), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 2.19989), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red"))+ 
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "", color = "", y = "Density", title = "B. Nigeria (2018): 80.9% VLS") +
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
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.257307),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 6.176566)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1, linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.8940994), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = 0.06680862), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 1.00765), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 1.695312), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red")) + 
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "", color = "", y = "Density", title = "C. Lesotho (2016-2017): 87.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
les

nam<-ssa_surveys_art %>% filter(country == "Namibia (2017)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  geom_histogram(alpha = 0.3, position = "identity", 
                 aes(y = after_stat(density)), 
                 color = "black", binwidth = 0.1) + 
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.048608),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 7.043988)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1, linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.7306145), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = -0.1121624), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 0.8262175), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 1.431562), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red")) + 
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title ="D. Namibia (2017): 91.3% VLS") +
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
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.74435),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 9.53912)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1, linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.4531612), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = -0.4804875), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 0.5225738), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 1.043111), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red")) + 
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title = "E. Eswatini (2021): 96.2% VLS") +
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
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.6111164),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dotted") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 11.83286)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1, linetype = "dotted") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.3209482), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dotted") +
  stat_function(fun = dlnorm, args = list(sdlog = 0.89, meanlog = -0.712478), 
                aes(colour = "Lognormal"), linewidth = 1) + 
  stat_function(fun = dfrechet, args = list(shape = 1.86, scale = 0.3773674), 
                aes(colour = "Frechet"), linewidth = 1) + 
  stat_function(fun = dgamma, args = list(shape = 0.81, scale = 0.8707441), 
                aes(colour = "Gamma"), linewidth = 1) + 
  scale_color_manual(values = c( "black", "yellow", "green","#4292C6", "orange","red")) +
  scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
                                expression(10^{5}),expression(10^{6}))) +
  coord_cartesian(ylim = c(0,.5), xlim = c(1.7, 6)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "Density", title = "F. Botswana (2021): 97.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
bots

cowplot::plot_grid(civ,nig,les,nam,esw,bots, rel_widths = c(.7,.7,1))


# number of PLHIV on ART in included surveys
allsurveys <- ssa_surveys_art %>% 
  group_by(country) %>%
  summarise(n = n())

# review max and min VLS ests
PHIA_vls_estimates[which.max(PHIA_vls_estimates$estimate_1000),"country"]
PHIA_vls_estimates[which.max(PHIA_vls_estimates$estimate_1000),"estimate_1000"]

PHIA_vls_estimates[which.min(PHIA_vls_estimates$estimate_1000),"country"]
PHIA_vls_estimates[which.min(PHIA_vls_estimates$estimate_1000),"estimate_1000"]


vls_est <- data.frame(survey = PHIA_vls_estimates$survey,
           VLS = paste0(round(PHIA_vls_estimates$estimate_1000,3)*100, " (", 
                        round(PHIA_vls_estimates$lower_ci_1000,3)*100, ", ",
                        round(PHIA_vls_estimates$upper_ci_1000,3)*100, ")"))


# code for appendix table S1
table(ssa_surveys_art$country, useNA = "always")

# lower limit of detect of surveys
head(table(ssa_surveys_art$resultvlc, useNA = "always"), n = 20)
tail(table(ssa_surveys_art$resultvlc, useNA = "always"), n = 20)
ssa_surveys_art %>% 
  filter(resultvlc %in% c("Target Not Detected")) %>%
  select(country) %>%
  distinct()

table(ssa_surveys_art$country)
ssa_surveys_art %>% 
  group_by(country) %>%
  filter(resultvlc %in% c("Target Not Detected")) %>%
  summarise(n())

# code for table S2
# countries reporting different thresholds in spectrum
spectrum_vl <- read.csv(file = "C:/Users/oedun/OneDrive - Imperial College London/R projects/viremia_thresholds/data/spectrum-vl-entered-data-2024_2024-11-19.csv")

names(spectrum_vl)
table(spectrum_vl$is_restricted, useNA = "always")
table(spectrum_vl$year, useNA = "always")

# four countries with multiple spectrum files
spectrum_vl %>% 
  filter(is_restricted == FALSE & year == "2023") %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  filter(n > 2)

# excluding the four countries we have 110 countries reporting in 2023, 21 reported at other thresholds
spectrum_vl %>% 
  filter(is_restricted == FALSE & year == "2023") %>%
  filter(!country %in% c("Ethiopia","Kenya","Zimbabwe","Republic of Moldova")) %>%
  group_by(vls_threshold) %>%
  summarise(n_distinct(country))

# all four countries report at <1000
spectrum_vl %>% 
  filter(is_restricted == FALSE & year == "2023") %>%
  filter(country %in% c("Ethiopia","Kenya","Zimbabwe","Republic of Moldova")) %>%
  select(vls_threshold, country, vl_suppressed_female15pl)


# accounting for countries with actual data
# excluding the four countries we have 91 countries reporting in 2023, 13 reported at other thresholds
spectrum_vl %>% 
  filter(is_restricted == FALSE & year == "2023") %>%
  filter(!country %in% c("Ethiopia","Kenya","Zimbabwe","Republic of Moldova")) %>%
  filter(!is.na(vl_suppressed_female15pl)|!is.na(vl_suppressed_male15pl)|!is.na(vl_suppressed_child)) %>%
  group_by(vls_threshold) %>%
  summarise(n_distinct(country))

# grouped by year
spectrum_vl %>% 
  filter(is_restricted == FALSE) %>%
  filter(!country %in% c("Ethiopia","Kenya","Zimbabwe","Republic of Moldova")) %>%
  filter(!is.na(vl_suppressed_female15pl)|!is.na(vl_suppressed_male15pl)|!is.na(vl_suppressed_child)) %>%
  group_by(year, vls_threshold) %>%
  summarise(n_distinct(country)) %>% print(n = Inf)

spectrum_vl %>% 
  filter(is_restricted == FALSE & year == "2023") %>%
  filter(!country %in% c("Ethiopia","Kenya","Zimbabwe","Republic of Moldova")) %>%
  filter(is.na(vl_suppressed_female15pl)&is.na(vl_suppressed_male15pl)&is.na(vl_suppressed_child)) %>%
  select(vls_threshold, country, vl_suppressed_female15pl, vl_suppressed_male15pl, vl_suppressed_child,
         vl_tested_child,vl_tested_female15pl,vl_tested_male15pl)

# Figure S9
library(dplyr)
library(ggplot2)
library(EnvStats)

civ <- ssa_surveys_art %>% filter(country == "Côte d'Ivoire (2017-18)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 2.137698),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x) {1 - pweibull(6 - x, shape = 2.81, scale = 4.573402) },
                aes(colour = "Reverse Weibull"), linewidth = 1) +
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 1.387628), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 2.185994),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 4.464708)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.9871715), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "", color = "", y = "CDF",title = "A. Côte d'Ivoire (2017-18): 73.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
civ

nig <- ssa_surveys_art %>% filter(country == "Nigeria (2018)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 1.65945),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.81, scale = 5.207989)},
                aes(colour = "Reverse Weibull"),linewidth = 1) + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 1.153083),
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 1.725519),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 5.046664)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.7558794), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "", color = "", y = "CDF", title = "B. Nigeria (2018): 80.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nig

les <- ssa_surveys_art %>% filter(country == "Lesotho (2016-17)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 1.257307),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.81, scale = 6.176566)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 0.8940994), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 1.33149),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 5.927286)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.5238001), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "", color = "", y = "CDF", title = "C. Lesotho (2016-2017): 87.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
les

nam <- ssa_surveys_art %>% filter(country == "Namibia (2017)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 1.048608),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.81, scale = 7.043988)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 0.7306145), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 1.123896),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 6.709182)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.391511), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "Viral load count (copies/mL)", color = "", y = "CDF", title ="D. Namibia (2017): 91.3% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nam

esw <- ssa_surveys_art %>% filter(country == "Eswatini (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 0.74435),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.81, scale = 9.53912)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 0.4531612), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 0.8160148),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 8.929944)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.1966564), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "Viral load count (copies/mL)", color = "", y = "CDF", title = "E. Eswatini (2021): 96.2% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
esw

bots <- ssa_surveys_art %>% filter(country == "Botswana (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_ecdf(geom = "step", aes(colour = "Survey"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.85, scale = 0.6111164),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.81, scale = 11.83286)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.73, location = 0.3209482), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = pweibull, args = list(shape = 0.91, scale = 0.6787366),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){1 - pweibull(6 - x, shape = 2.98, scale = 10.94185)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::ppareto, args = list(shape = 1.20, location = 0.1196038), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(breaks = 1:6,
                     labels = c(expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4), 
                                expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1)) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey") +
  labs(x = "Viral load count (copies/mL)", color = "", y = "CDF", title = "F. Botswana (2021): 97.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
bots

cowplot::plot_grid(civ,nig,les,nam,esw,bots, rel_widths = c(.7,.7,1))

# Figure S10 showing PDFs
civ <- ssa_surveys_art %>% filter(country == "Côte d'Ivoire (2017-18)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 2.137698),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 4.573402)},
                aes(colour = "Reverse Weibull"), linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.387628), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 2.185994),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 4.464708)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.9871715), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  #scale_x_continuous(labels = c(expression(10^{2}), expression(10^{3}),expression(10^{4}),
  #                              expression(10^{5}),expression(10^{6}))) +
  #coord_cartesian(ylim = c(0,1), xlim = c(1.7, 6)) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 1)) +
  labs(x = "", color = "", y = "PDF",title = "A. Côte d'Ivoire (2017-18): 73.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
civ

nig <- ssa_surveys_art %>% filter(country == "Nigeria (2018)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.65945),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 5.207989)},
                aes(colour = "Reverse Weibull"),linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 1.153083),
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.725519),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 5.046664)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.7558794), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 1.1)) +
  labs(x = "", color = "", y = "PDF", title = "B. Nigeria (2018): 80.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nig

les <- ssa_surveys_art %>% filter(country == "Lesotho (2016-17)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.257307),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 6.176566)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.8940994), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.33149),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 5.927286)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.5238001), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 5.2)) +
  labs(x = "", color = "", y = "PDF", title = "C. Lesotho (2016-2017): 87.7% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
les

nam <- ssa_surveys_art %>% filter(country == "Namibia (2017)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 1.048608),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 7.043988)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.7306145), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 1.123896),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 6.709182)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.391511), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 4.2)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "PDF", title ="D. Namibia (2017): 91.3% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))

nam

esw <- ssa_surveys_art %>% filter(country == "Eswatini (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.74435),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 9.53912)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.4531612), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 0.8160148),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 8.929944)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.1966564), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 2.55)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "PDF", title = "E. Eswatini (2021): 96.2% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
esw

bots <- ssa_surveys_art %>% filter(country == "Botswana (2021)") %>%
  ggplot(aes(log10(resultvlc_numeric))) + 
  stat_density(aes(colour = "Survey"), geom = "line", linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.85, scale = 0.6111164),
                aes(colour = "Weibull"), linewidth = 1) + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.81, scale = 11.83286)},
                aes(colour = "Reverse Weibull"),
                linewidth = 1) + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.73, location = 0.3209482), 
                aes(colour = "Pareto"), linewidth = 1) +
  stat_function(fun = dweibull, args = list(shape = 0.91, scale = 0.6787366),
                aes(colour = "Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = function(x){dweibull(6-x, shape = 2.98, scale = 10.94185)},
                aes(colour = "Reverse Weibull"), linewidth = 1, linetype = "dashed") + 
  stat_function(fun = EnvStats::dpareto, args = list(shape = 1.20, location = 0.1196038), 
                aes(colour = "Pareto"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Survey" = "black",
                                "Pareto" = "#4292C6",
                                "Reverse Weibull" = "orange",
                                "Weibull" = "red"),
                     breaks = c("Survey", "Weibull", "Reverse Weibull", "Pareto")) + 
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2),
               expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  coord_cartesian(xlim = c(1, 6), ylim = c(0,1.5)) +
  ggbreak::scale_y_break(c(.75, 17)) +
  labs(x = "Viral load count (copies/mL)", color = "", y = "PDF", title = "F. Botswana (2021): 97.9% VLS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = rel(1.0), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans",face = "bold"),
        legend.text = element_text(size = rel(1.2), family = "sans"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1), family = "sans", face = "bold"))
bots

cowplot::plot_grid(civ,nig,les,nam,esw,bots, rel_widths = c(.7,.7,1))
