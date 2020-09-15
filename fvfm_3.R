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

RNG = "A1:G941"
sheet = "200821_wakame_168h"
dset = read_xlsx(dir(pattern = "200821_wakame_168h_2*.*xlsx",
                     full.names = TRUE),
                 range = RNG,
                 sheet = sheet)

#1個だけ項目名を変換したい
names(dset)[ which( names(dset)=="Fv/Fm" ) ] = "FvFm"

dset = dset %>% 
  select(Hour,Growth,temperature = Temp,fvfm =FvFm)

# data24 = dset %>% 
#   group_by(Growth) %>% 
#   filter(Hour == 24) %>% 
#   select(Hour,Growth,temperature = Temp,fvfm =FvFm)

# H,Gごとにデータ読み込む ----

data24M = dset %>% 
  filter(Hour == 24,Growth == "Mature")

data24J = dset %>% 
  filter(Hour == 24,Growth == "Juvenile")

data48M = dset %>% 
  filter(Hour == 48,Growth == "Mature") 

data48J = dset %>% 
  filter(Hour == 48,Growth == "Juvenile")

data72M = dset %>% 
  filter(Hour == 72,Growth == "Mature")

data72J = dset %>% 
  filter(Hour == 72,Growth == "Juvenile") 

data168M = dset %>% 
  filter(Hour == 168,Growth == "Mature")
data168J = dset %>% 
  filter(Hour == 168,Growth == "Juvenile")

dset %>% 
  ggplot() + 
  geom_point(aes(x=Temp, y = FvFm)) +
  facet_grid(rows = vars(Hour),
             cols = vars(Growth))

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
                        nl = TRUE,
                        family = brms::Beta(link = "identity",
                                            link_phi = "log"))

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

# 仮当てはめ
brmout = brm(brmsmodel, 
                data = dset %>% filter(temperature < 30),
                stanvars = stanvars, prior = PRIORS,
                iter = 10, chains = 4, cores = 4, seed = 2020,
                control = list(adapt_delta = 0.9999, max_treedepth = 12))

# Mature
dset2 = dset %>% 
  filter(Growth != "Juvenile") %>% 
  filter(temperature < 30) %>% 
  select(Hour,temperature,fvfm) %>% 
  group_nest(Hour) %>% 
  mutate(bfit = map(data, function(df) {
    update(brmout, 
           newdata = df, 
           stanvars = stanvars, prior = PRIORS,
           iter = 5000, chains = 4, cores = 4, seed = 2020,
           control = list(adapt_delta = 0.9999, max_treedepth = 12))
  }))

# Juvenile
dset2J = dset %>% 
  filter(Growth == "Juvenile") %>% 
  filter(temperature < 30) %>% 
  select(Hour,temperature,fvfm) %>% 
  group_nest(Hour) %>% 
  mutate(bfit = map(data, function(df) {
    update(brmout, 
           newdata = df, 
           stanvars = stanvars, prior = PRIORS,
           iter = 5000, chains = 4, cores = 4, seed = 2020,
           control = list(adapt_delta = 0.9999, max_treedepth = 12))
  }))

expose_functions(brmout24M, vectorize = TRUE)

dset2 = dset2 %>% 
  mutate(summary = map(bfit, summary)) %>% 
  mutate(ppc_dens_overlay = map(bfit, function(bfit) {
    y = bfit$data %>% pull(fvfm)
    yrep = posterior_predict(bfit, nsamples = 500)
    ppc_dens_overlay(y, yrep)
  })) %>% 
    mutate(mcmc_rank_overlay = map(bfit, function(bfit){
      posterior_out = as.array(bfit)
      mcmc_rank_overlay(posterior_out)      
    }))

dset2 = dset2 %>% 
  mutate(predict = map(bfit, function(bfit) {
    bfit$data %>% 
      expand(temperature = seq(min(temperature), 
                               max(temperature),
                               by = 1)) %>% 
      add_predicted_draws(bfit) %>% 
      group_by(temperature) %>% 
      mean_hdci()
  })) %>% 
  mutate(expect = map(bfit, function(bfit) {
    bfit$data %>% 
      expand(temperature = seq(min(temperature), 
                               max(temperature),
                               by = 1)) %>% 
      add_linpred_draws(bfit) %>% 
      group_by(temperature) %>% 
      mean_hdci()
  }))


dset2J = dset2J %>% 
  mutate(summary = map(bfit, summary)) %>% 
  mutate(ppc_dens_overlay = map(bfit, function(bfit) {
    y = bfit$data %>% pull(fvfm)
    yrep = posterior_predict(bfit, nsamples = 500)
    ppc_dens_overlay(y, yrep)
  })) %>% 
  mutate(mcmc_rank_overlay = map(bfit, function(bfit){
    posterior_out = as.array(bfit)
    mcmc_rank_overlay(posterior_out)      
  }))

dset2J = dset2J %>% 
  mutate(predict = map(bfit, function(bfit) {
    bfit$data %>% 
      expand(temperature = seq(min(temperature), 
                               max(temperature),
                               by = 1)) %>% 
      add_predicted_draws(bfit) %>% 
      group_by(temperature) %>% 
      mean_hdci()
  })) %>% 
  mutate(expect = map(bfit, function(bfit) {
    bfit$data %>% 
      expand(temperature = seq(min(temperature), 
                               max(temperature),
                               by = 1)) %>% 
      add_linpred_draws(bfit) %>% 
      group_by(temperature) %>% 
      mean_hdci()
  }))

xlabF = expression("Temperature" ~ ({}^degree*C ))
ylabF = "Maximum Quantum Yield (Fv/Fm)" 


plot = ggplot() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper,x = temperature),
              alpha = 0.2,
              data = dset2 %>% unnest(predict),
              fill = "seagreen") +
  geom_line(aes(x = temperature, y = .value,colour = "seagreen"),
            data = dset2 %>% unnest(expect)) +
  geom_point(aes(x = temperature, y = fvfm, colour = "seagreen"),
             data = dset2 %>% unnest(data)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper,x = temperature),
              alpha = 0.2,
              data = dset2J %>% unnest(predict),
              fill = "salmon") +
  geom_line(aes(x = temperature, y = .value,colour = "salmon"),
            data = dset2J %>% unnest(expect)) +
  geom_point(aes(x = temperature, y = fvfm,colour = "salmon"),
             data = dset2J %>% unnest(data)) +
  facet_wrap("Hour")+
  theme_pubr()

plot

  ggtitle("a")+
  annotate("text",x = 35, y = 0.95,label = "24-h")+
  scale_x_continuous(xlabF,limits = c(4, 36),breaks = seq(4, 36, by = 4)) +
  scale_y_continuous(ylabF,limits = c(0, 1.0),breaks = seq(0,1.0, by = 0.2)) +
  theme(legend.position = c(0.1,0.2),
        legend.background = element_blank(),
        plot.title = element_text(vjust = -6,hjust = 0.02))

ggsave("plotM.jpg",plot = plotM)
