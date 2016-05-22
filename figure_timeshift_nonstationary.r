library('ggplot2'); library('dplyr'); library('tidyr');
theme_set(theme_bw(12))

### load data
load('figure_timeshift_nonstationary.rdata')

### data wrangling
output_p = output %>%
	rename(F = F_mean, E = E_mean, Ehat = Ehat_mean_mean, M2 = M2_mean, M2hat = M2hat_mean) #%>%

output_t0 = output_p %>% filter(t==0) %>% 
	select(name, sa, sv, theta_a, theta_v, F) %>%
	spread(name, F) %>% 
	rename(Fv=virus) %>% 
	inner_join(output_p %>% filter(t==0), by=c("sa","sv","theta_a","theta_v"))

qplot(data = output_p %>% filter(name=="virus"), x = t, y = F, color = factor(paste(sa,sv)), geom="line") + 
	geom_segment(aes(x = 0, xend = -25, y = Fv, yend = Fv + (sv^2*M2-2*theta_v*Fv)*25, 
					 color = factor(paste(sa,sv))), 
				 data=output_t0 %>% filter(name=="virus"),
				 linetype=2) +
	geom_segment(aes(x = 0, xend = 25, y = Fv, yend = Fv - (sa*sv*M2+sa*sv*M2hat+2*theta_a*Fv)*25,
					 color=factor(paste(sa,sv))),
				 data=output_t0 %>% filter(name=="antibody"),
				 linetype=2) + 
	guides(color = "none")
#+
#	coord_cartesian(ylim=c(-2,2.5))
	
ggsave('figures/timeshift_nonstationary.pdf', height=3, width=4)

qplot(data = output_p %>% filter(name=="virus"), x = t, y = E, color = factor(paste(sa,sv)), geom="line") + 
	geom_segment(aes(x = 0, xend = -25, y = E, yend = E - (sv*M2-2*theta_v*E)*25, 
					 color = factor(paste(sa,sv))), 
				 data=output_t0 %>% filter(name=="virus"),
				 linetype=2) +
	geom_segment(aes(x = 0, xend = 25, y = E, yend = E + (sa*M2+2*theta_a*E)*25,
					 color=factor(paste(sa,sv))),
				 data=output_t0 %>% filter(name=="antibody"),
				 linetype=2) + 
	guides(color = "none")
ggsave('figures/timeshift_nonstationary_E.pdf', height = 3, width = 4)

qplot(data = output_p %>% filter(name=="virus"), x = t, y = Ehat, color = factor(paste(sa,sv)), geom="line") + 
	geom_segment(aes(x = 0, xend = 25, y = Ehat, yend = Ehat + (sa*M2hat+2*theta_a*Ehat)*25,
					 color=factor(paste(sa,sv))),
				 data=output_t0 %>% filter(name=="antibody"),
				 linetype=2) + 
	guides(color = "none")
ggsave('figures/timeshift_nonstationary_Ehat.pdf', height = 3, width = 4)
