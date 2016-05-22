library('ggplot2'); library('dplyr'); library('tidyr');
theme_set(theme_bw(12))

### load data
load('figure_nonlinear_fitness.rdata')
result_stationary_p = result_stationary %>% 
	filter(name == "antibody") %>%
	mutate(logM_mean = signif(logM_mean,3), beta_mean = signif(beta_mean,2))

## mean binding
qplot(data = result_stationary_p %>% filter(logM_mean >0), 
	  x=sa, y = E_mean, color = factor(logM_mean), 
	  facets=beta_mean~Estar_mean, geom="line") + 
	geom_line(data = result_stationary_p %>% filter(logM_mean ==0), color="black")
ggsave('figures/nonlinear_fitness_E.pdf', width=6, height=5)


## fitness variances

# linear approximation
result_stationary_lin = result_stationary_p %>% 
	filter(logM_mean == 0, beta_mean == 0.1, Estar_mean == 0) %>%
	select(sa, sv, M2_mean) %>% rename(M2lin = M2_mean) %>%
	inner_join(result_stationary_p %>% filter(logM_mean==0), by=c("sa","sv"))

qplot(data = result_stationary_p %>% filter(logM_mean >0),
	  x=sa, y = F2_mean, color = factor(logM_mean),
	  facets=beta_mean~Estar_mean, geom="line") + 
	geom_line(data = result_stationary_p %>% filter(logM_mean ==0), color="black") +
	geom_line(data = result_stationary_lin, aes(y=sa^2*M2lin),
			  linetype=2, color="black")
ggsave('figures/nonlinear_fitness_F2.pdf', width=6, height=5)
