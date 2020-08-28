# Aug 21 2020
# Aug 27 2020
# Hikari Nagoe
# Fv / Fm - T

library(tidyverse)
library(readxl)
library(rstan)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggpubr)
library(knitr)

rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

# ラベル ----
xlabF = expression("Temperature" ~ ({}^degree*C ))
ylabF = "Maximum Quantum Yield (Fv/Fm)" 

# read "publushed parameters" ----

priors = read_xlsx("published_parameters.xlsx")

priors = priors %>% group_nest(reference, type, species) %>%
  unnest(cols = c(data)) %>%
  filter(str_detect(parameter, "hd", negate = T)) %>%
  group_by(parameter, type) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  arrange(type, parameter)

oxygen_priors = priors %>% filter(type == "oxygen")
pam_priors = priors %>% filter(type == "pam")


# 基本データ読み込み----

RNG = "A1:G831"
sheet = "200821_wakame_168h"
dset = read_xlsx(dir(pattern = "200821_wakame_168h*.*xlsx",
                     full.names = TRUE),
                 range = RNG,
                 sheet = sheet)

#1個だけ項目名を変換したい
names(dset)[ which( names(dset)=="Fv/Fm" ) ] = "FvFm"

# data24 = dset %>% 
#   group_by(Growth) %>% 
#   filter(Hour == 24) %>% 
#   select(Hour,Growth,temperature = Temp,fvfm =FvFm)

# H,Gごとにデータ読み込む ----

data24M = dset %>% 
  filter(Hour == 24,Growth == "Mature") %>% 
  select(Hour,Growth,temperature = Temp,fvfm =FvFm)

data24J = dset %>% 
  filter(Hour == 24,Growth == "Juvenile") %>% 
  select(Hour,Growth,temperature = Temp,fvfm =FvFm)

data48M = dset %>% 
  filter(Hour == 48,Growth == "Mature") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)

data48J = dset %>% 
  filter(Hour == 48,Growth == "Juvenile") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)

data72M = dset %>% 
  filter(Hour == 72,Growth == "Mature") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)

data72J = dset %>% 
  filter(Hour == 72,Growth == "Juvenile") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)

data168M = dset %>% 
  filter(Hour == 168,Growth == "Mature") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)

data168J = dset %>% 
  filter(Hour == 168,Growth == "Juvenile") %>% 
  select(Hour,temperature = Temp,fvfm =FvFm)


# ここでstan ----

stan_function =
  "
real fvfmmodel (real ps, real ha, real eta, real ktopt, real temperature) {
real invkelvin = 1.0 / (temperature + 273.15);
real gas_constant = 8.314/1000.0;
real tmp0 = (1.0 / ktopt - invkelvin);
real tmp1 = ha / gas_constant * tmp0 ;
real tmp2 = (ha * eta) / gas_constant * tmp0;
tmp1 = (ha * eta) * exp(tmp1);
tmp2 = 1 - exp(tmp2);
tmp2 = (ha * eta) - ha * tmp2;
return ps * tmp1 / tmp2;
}
"
stanvars = stanvar(scode = stan_function, block = "functions")

# brms model 作成 ----
brmsmodel = brmsformula(fvfm ~ fvfmmodel(PS, HA, ET, KT, temperature),
                        PS ~ 1,
                        HA ~ 1,
                        ET ~ 1,
                        KT ~ 1,
                        zi ~ temperature,
                        nl = TRUE,
                        family = brms::zero_inflated_beta(link = "identity",
                                                          link_phi = "log",
                                                          link_zi= "logit"))

# パラメーターの設定 ----
prs = pam_priors %>%
  mutate(mean = ifelse(parameter == "topt", mean + 273.15, mean)) %>%
  mutate(prior = str_glue("normal({mean}, {sd})")) %>%
  select(parameter, prior)

PRIORS = set_prior(prs %>% filter(parameter == "eta") %>% pull(prior),
                   class = "b", lb = 0, nlpar = "ET") +
  set_prior(prs %>% filter(parameter == "topt") %>% pull(prior),
            class = "b", lb = 0, nlpar = "KT") +
  set_prior(prs %>% filter(parameter == "ymax") %>% pull(prior),
            class = "b", lb = 0, nlpar = "PS") +
  set_prior(prs %>% filter(parameter == "ha") %>% pull(prior),
            class = "b", lb = 0, nlpar = "HA")

# # 24
# まとめてやりたかったけどできなかった。
# brmsmodel2 = brmsformula(fvfm ~ fvfmmodel(PS, HA, ET, KT, temperature),
#                         PS ~ (1|Growth),
#                         HA ~ (1|Growth),
#                         ET ~ (1|Growth),
#                         KT ~ (1|Growth),
#                         zi ~ temperature,
#                         nl = TRUE,
#                         family = brms::zero_inflated_beta(link = "identity",
#                                                           link_phi = "log",
#                                                           link_zi= "logit"))
# 
# brmout24 = brm(brmsmodel2, data = data24,
#                  stanvars = stanvars, prior = PRIORS,
#                  iter = 3000, chains = 4, cores = 4, seed = 1500,
#                  control = list(adapt_delta = 0.999, max_treedepth = 13))
# 
# summary(brmout24)
# 
# posterior_out24 = as.array(brmout24)
# mcmc_rank_overlay(posterior_out24)
# mcmc_trace(brmout24)
# 
# expose_functions(brmout24, vectorize = TRUE)
# y24 = brmout24$data %>% pull(fvfm)
# yrep24 = posterior_predict(brmout24, nsamples = 500)
# ppc_dens_overlay(y24, yrep24)
# ppc_hist(y24, yrep24[1:8,])
# xval24 = data24$temperature
# ppc_error_scatter_avg_vs_x(y24, yrep24, xval24) +
#   geom_hline(yintercept = 0, linetype = "dashed")
# 
# pdata2400 = data24 %>%
#   expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
#   add_predicted_draws(brmout24) %>%
#   group_by(temperature) %>%
#   mean_hdci()
# 
# edata24 = data24 %>%
#   expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
#   add_linpred_draws(brmout24) %>%
#   group_by(temperature) %>%
#   mean_hdci()
# 


# 24M ----
brmout24M = brm(brmsmodel, data = data24M,
                stanvars = stanvars, prior = PRIORS,
                iter = 10000, chains = 4, cores = 4, seed = 5050,
                control = list(adapt_delta = 0.9999, max_treedepth = 14))

summary(brmout24M)

posterior_out24M = as.array(brmout24M)
mcmc_rank_overlay(posterior_out24M)
mcmc_trace(brmout24M)

expose_functions(brmout24M, vectorize = TRUE)
y24M = brmout24M$data %>% pull(fvfm)
yrep24M = posterior_predict(brmout24M, nsamples = 500)
ppc_dens_overlay(y24M, yrep24M)
ppc_hist(y24M, yrep24M[1:8,])
xval24M = data24M$temperature
ppc_error_scatter_avg_vs_x(y24M, yrep24M, xval24M) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata24M = data24M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout24M) %>%
  group_by(temperature) %>%
  mean_hdci()

edata24M = data24M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout24M) %>%
  group_by(temperature) %>%
  mean_hdci()


# 24J ----
brmout24J = brm(brmsmodel, data = data24J,
                stanvars = stanvars, prior = PRIORS,
                iter = 10000, chains = 4, cores = 4, seed = 5050,
                control = list(adapt_delta = 0.999, max_treedepth = 13))


posterior_out24J = as.array(brmout24J)
mcmc_rank_overlay(posterior_out24J)

expose_functions(brmout24J, vectorize = TRUE)
y24J = brmout24J$data %>% pull(fvfm)
yrep24J = posterior_predict(brmout24J, nsamples = 500)
ppc_dens_overlay(y24J, yrep24J)
ppc_hist(y24J, yrep24J[1:8,])
xval24J = data24J$temperature
ppc_error_scatter_avg_vs_x(y24J, yrep24J, xval24J) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata24J = data24J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout24J) %>%
  group_by(temperature) %>%
  mean_hdci()

edata24J = data24J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout24J) %>%
  group_by(temperature) %>%
  mean_hdci()

# 作図 24 ----
cols = c("L1" = "seagreen","L2" = "salmon")

plot24 = ggplot() +
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata24M,alpha = 0.2,fill = "seagreen") +
  geom_line(aes(x = temperature,y = .value,colour = "L1"),data = edata24M) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L1"),data = data24M) +  
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata24J,alpha = 0.2,fill = "salmon") +
  geom_line(aes(x = temperature,y = .value,color = "L2"),data = edata24J) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L2"),data = data24J) +
  ggtitle("a")+
  annotate("text",x = 35, y = 0.95,label = "24-h")+
  scale_x_continuous(xlabF,limits = c(4, 36),breaks = seq(4, 36, by = 4)) +
  scale_y_continuous(ylabF,limits = c(0, 1.0),breaks = seq(0,1.0, by = 0.2)) +
  scale_colour_manual(name="",values=cols,labels = c("Mature","Juvenile"))+
  theme_pubr()+
  theme(legend.position = c(0.1,0.2),
        legend.background = element_blank(),
        plot.title = element_text(vjust = -6,hjust = 0.02))

plot24  

# 48M ----
brmout48M = brm(brmsmodel, data = data48M,
                stanvars = stanvars, prior = PRIORS,
                iter = 10000, chains = 4, cores = 4, seed = 5050,
                control = list(adapt_delta = 0.99999, max_treedepth = 14))

summary(brmout48M)

posterior_out48M = as.array(brmout48M)
mcmc_rank_overlay(posterior_out48M)
mcmc_trace(brmout48M)

expose_functions(brmout48M, vectorize = TRUE)
y48M = brmout48M$data %>% pull(fvfm)
yrep48M = posterior_predict(brmout48M, nsamples = 500)
ppc_dens_overlay(y48M, yrep48M)
ppc_hist(y48M, yrep48M[1:8,])
xval48M = data48M$temperature
ppc_error_scatter_avg_vs_x(y48M, yrep48M, xval48M) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata48M = data48M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout48M) %>%
  group_by(temperature) %>%
  mean_hdci()

edata48M = data48M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout48M) %>%
  group_by(temperature) %>%
  mean_hdci()


# 48J ----
brmout48J = brm(brmsmodel, data = data48J,
                stanvars = stanvars, prior = PRIORS,
                iter = 7000, chains = 4, cores = 4, seed = 3456,
                control = list(adapt_delta = 0.999, max_treedepth = 12))


posterior_out48J = as.array(brmout48J)
mcmc_rank_overlay(posterior_out48J)

expose_functions(brmout48J, vectorize = TRUE)
y48J = brmout48J$data %>% pull(fvfm)
yrep48J = posterior_predict(brmout48J, nsamples = 500)
ppc_dens_overlay(y48J, yrep48J)
ppc_hist(y48J, yrep48J[1:8,])
xval48J = data48J$temperature
ppc_error_scatter_avg_vs_x(y48J, yrep48J, xval48J) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata48J = data48J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout48J) %>%
  group_by(temperature) %>%
  mean_hdci()

edata48J = data48J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout48J) %>%
  group_by(temperature) %>%
  mean_hdci()

# 作図 48 ----

plot48 = ggplot() +
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata48M,alpha = 0.2,fill = "seagreen") +
  geom_line(aes(x = temperature,y = .value,colour = "L1"),data = edata48M) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L1"),data = data48M) +  
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata48J,alpha = 0.2,fill = "salmon") +
  geom_line(aes(x = temperature,y = .value,color = "L2"),data = edata48J) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L2"),data = data48J) +
  ggtitle("a")+
  annotate("text",x = 35, y = 1.0,label = "24-h")+
  scale_x_continuous(xlabF,limits = c(4, 36),breaks = seq(4, 36, by = 4)) +
  scale_y_continuous(ylabF,limits = c(0, 0.8),breaks = seq(0,0.8, by = 0.2))+
  theme_pubr()+
  theme(plot.title = element_text(vjust = -6,hjust = 0.02))

plot48

# 72M ----
brmout72M = brm(brmsmodel, data = data72M,
                stanvars = stanvars, prior = PRIORS,
                iter = 7000, chains = 4, cores = 4, seed = 3456,
                control = list(adapt_delta = 0.999, max_treedepth = 13))

summary(brmout72M)

posterior_out72M = as.array(brmout72M)
mcmc_rank_overlay(posterior_out72M)
mcmc_trace(brmout72M)

expose_functions(brmout72M, vectorize = TRUE)
y72M = brmout72M$data %>% pull(fvfm)
yrep72M = posterior_predict(brmout72M, nsamples = 500)
ppc_dens_overlay(y72M, yrep72M)
ppc_hist(y72M, yrep72M[1:8,])
xval72M = data72M$temperature
ppc_error_scatter_avg_vs_x(y72M, yrep72M, xval72M) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata72M = data72M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout72M) %>%
  group_by(temperature) %>%
  mean_hdci()

edata72M = data72M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout72M) %>%
  group_by(temperature) %>%
  mean_hdci()


# 72J ----
brmout72J = brm(brmsmodel, data = data72J,
                stanvars = stanvars, prior = PRIORS,
                iter = 7000, chains = 4, cores = 4, seed = 3456,
                control = list(adapt_delta = 0.999, max_treedepth = 12))


posterior_out72J = as.array(brmout72J)
mcmc_rank_overlay(posterior_out72J)

expose_functions(brmout72J, vectorize = TRUE)
y72J = brmout72J$data %>% pull(fvfm)
yrep72J = posterior_predict(brmout72J, nsamples = 500)
ppc_dens_overlay(y72J, yrep72J)
ppc_hist(y72J, yrep72J[1:8,])
xval72J = data72J$temperature
ppc_error_scatter_avg_vs_x(y72J, yrep72J, xval72J) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata72J = data72J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout72J) %>%
  group_by(temperature) %>%
  mean_hdci()

edata72J = data72J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout72J) %>%
  group_by(temperature) %>%
  mean_hdci()

# 作図 72 ----

plot72 = ggplot() +
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata72M,alpha = 0.2,fill = "seagreen") +
  geom_line(aes(x = temperature,y = .value,colour = "L1"),data = edata72M) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L1"),data = data72M) +  
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata72J,alpha = 0.2,fill = "salmon") +
  geom_line(aes(x = temperature,y = .value,color = "L2"),data = edata72J) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L2"),data = data72J) +
  ggtitle("a")+
  annotate("text",x = 35, y = 1.0,label = "24-h")+
  scale_x_continuous(xlabF,limits = c(4, 36),breaks = seq(4, 36, by = 4)) +
  scale_y_continuous(ylabF,limits = c(0, 0.8),breaks = seq(0,0.8, by = 0.2)) +
  theme_pubr()+
  theme(plot.title = element_text(vjust = -6,hjust = 0.02))

plot72

# 168M ----
brmout168M = brm(brmsmodel, data = data168M,
                 stanvars = stanvars, prior = PRIORS,
                 iter = 7000, chains = 4, cores = 4, seed = 3456,
                 control = list(adapt_delta = 0.999, max_treedepth = 13))

summary(brmout168M)

posterior_out168M = as.array(brmout168M)
mcmc_rank_overlay(posterior_out168M)
mcmc_trace(brmout168M)

expose_functions(brmout168M, vectorize = TRUE)
y168M = brmout168M$data %>% pull(fvfm)
yrep168M = posterior_predict(brmout168M, nsamples = 500)
ppc_dens_overlay(y168M, yrep168M)
ppc_hist(y168M, yrep168M[1:8,])
xval168M = data168M$temperature
ppc_error_scatter_avg_vs_x(y168M, yrep168M, xval168M) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata168M = data168M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout168M) %>%
  group_by(temperature) %>%
  mean_hdci()

edata168M = data168M %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout168M) %>%
  group_by(temperature) %>%
  mean_hdci()

# 168J ----
brmout168J = brm(brmsmodel, data = data168J,
                 stanvars = stanvars, prior = PRIORS,
                 iter = 7000, chains = 4, cores = 4, seed = 3456,
                 control = list(adapt_delta = 0.999, max_treedepth = 12))


posterior_out168J = as.array(brmout168J)
mcmc_rank_overlay(posterior_out168J)

expose_functions(brmout168J, vectorize = TRUE)
y168J = brmout168J$data %>% pull(fvfm)
yrep168J = posterior_predict(brmout168J, nsamples = 500)
ppc_dens_overlay(y168J, yrep168J)
ppc_hist(y168J, yrep168J[1:8,])
xval168J = data168J$temperature
ppc_error_scatter_avg_vs_x(y168J, yrep168J, xval168J) +
  geom_hline(yintercept = 0, linetype = "dashed")

pdata168J = data168J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_predicted_draws(brmout168J) %>%
  group_by(temperature) %>%
  mean_hdci()

edata168J = data168J %>%
  expand(temperature = seq(min(temperature), max(temperature), by = 1)) %>%
  add_linpred_draws(brmout168J) %>%
  group_by(temperature) %>%
  mean_hdci()

# 作図 168 ----

plot168 = ggplot() +
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata168M,alpha = 0.2,fill = "seagreen") +
  geom_line(aes(x = temperature,y = .value,colour = "L1"),data = edata168M) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L1"),data = data168M) +  
  geom_ribbon(aes(ymin = .lower,ymax = .upper,x = temperature),data = edata168J,alpha = 0.2,fill = "salmon") +
  geom_line(aes(x = temperature,y = .value,color = "L2"),data = edata168J) +
  geom_point(aes(x = temperature, y = fvfm,colour = "L2"),data = data168J) +
  ggtitle("a")+
  annotate("text",x = 35, y = 1.0,label = "24-h")+
  scale_x_continuous(xlabF,limits = c(4, 36),breaks = seq(4, 36, by = 4)) +
  scale_y_continuous(ylabF,limits = c(0, 0.8),breaks = seq(0,0.8, by = 0.2)) +
  theme_pubr()+
  theme(plot.title = element_text(vjust = -6,hjust = 0.02))

plot168