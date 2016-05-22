library('dplyr'); library('readr');library('tidyr');library('ggplot2')
theme_set(theme_bw(12))

load('figure_timeshift.rdata')

### data wrangling, put M2A M2V as new column
result_stationary_m2 = result_stationary %>%
	filter(t == 0) %>%
	select(name, sa, sv, theta_a, theta_v, lhat, M2_mean) %>%
	spread(name, M2_mean) %>%
	rename(M2A = antibody, M2V = virus) %>%
	inner_join(result_stationary %>% filter(name=='virus'), by = c("sa", "sv", "theta_a", "theta_v", "lhat"))

# write theory
result_stationary1 = result_stationary_m2 %>%
  filter(t <= 0, theta_a != theta_v) %>%
  mutate(E_theory = 1/2*M2V*sv/(theta_a-theta_v)*exp(2*theta_a*t) +
    1/2*(M2A*sa*(theta_a-theta_v)-2*(theta_a*M2V*sv))/(theta_a^2-theta_v^2)*exp(2*theta_v*t))
result_stationary2 = result_stationary_m2 %>%
  filter(t > 0, theta_a != theta_v) %>%
  mutate(E_theory = 1/2*M2A*sa/(theta_a-theta_v)*exp(-2*theta_v*t) +
    1/2*(M2V*sv*(theta_v-theta_a)-2*(theta_v*M2A*sa))/(theta_a^2-theta_v^2)*exp(-2*theta_a*t))
result_stationary3 = result_stationary_m2 %>%
  filter(t <= 0, theta_a == theta_v) %>%
  mutate(E_theory = ((M2A*sa-M2V*sv)/4/theta_a+M2V*sv*t)*exp(2*theta_a*t))
result_stationary4 = result_stationary_m2 %>%
  filter(t > 0, theta_a == theta_v) %>%
  mutate(E_theory = ((M2A*sa-M2V*sv)/4/theta_a+M2A*sa*t)*exp(-2*theta_a*t))
result_stationary_theory = bind_rows(result_stationary1, result_stationary2, result_stationary3, result_stationary4)

result_stationary_main = result_stationary_theory %>% filter(theta_a == 1/50, theta_v == 1/50, lhat == 50, (sa == 2 & sv ==2) | (sa==0 & sv==2) | (sa==2 & sv==0))
									
qplot(data = result_stationary_main, x = t, y = E_mean, color = paste0(sa,sv), geom="line") +
	geom_line(aes(y = E_theory), linetype=2) + guides(color = "none")
ggsave('figures/timeshift_model_E.pdf', height = 3, width = 3)

qplot(data = result_stationary_main, x = t, y = -sv*E_mean, color = paste0(sa,sv), geom="line") +
	geom_line(aes(y = -sv*E_theory), linetype=2) + 
	geom_segment(aes(x = 0, xend = -100, y = -sv*E_mean, yend = -sv*E_mean + (sv^2*M2V+2*theta_v*sv*E_mean)*100), 
				 data = result_stationary_main %>% filter(t==0), linetype=3) +
	geom_segment(aes(x = 0, xend = 100, y = -sv*E_mean, yend = -sv*E_mean - (sa*sv*M2A-2*theta_a*sv*E_mean)*100),
				 data = result_stationary_main %>% filter(t==0), linetype=3) +
	guides(color = "none") +
	coord_cartesian(ylim = c(-2.5, 4))
ggsave('figures/timeshift_model_F.pdf', height = 3, width = 3)

qplot(data = result_stationary_theory %>% filter(lhat == 50, sa >0 | sv > 0),
	  x = t, y = E_mean, color = paste0(sa,sv), geom="line", facets= theta_a~theta_v) + 
	geom_line(aes(y = E_theory), linetype=2)
ggsave('figures/timeshift_model_SI.pdf', height = 5, width = 8)


