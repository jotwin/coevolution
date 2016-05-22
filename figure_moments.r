library('ggplot2'); library('dplyr'); library('tidyr');
theme_set(theme_bw(12))

### load data
load('figure_moments.rdata')
n = 1000
l = 50
lhat = 50
theta_a = 1/50
theta_v = 1/50

### data wrangling
result_stationary_av = result_stationary %>% 
	filter(name == "antibody") %>%
	transmute(sa=sa, sv=sv, E=E_mean, M2A=M2_mean, M3A=M3_mean, M4A=M4_mean)
result_stationary_vv = result_stationary %>% 
	filter(name == "virus") %>%
	transmute(sa=sa, sv=sv, M2V=M2_mean, M3V=M3_mean, M4V=M4_mean)
result_stationary_var = inner_join(result_stationary_av, result_stationary_vv, by = c("sa", "sv")) %>%
	mutate(region = 1)
result_stationary_ac = result_stationary %>% 
	filter(name == "antibody") %>%
	transmute(sa=sa, sv=sv, E=Ehat_mean, M2A=M2hat_mean, M3A=M3hat_mean, M4A=M4hat_mean)
result_stationary_vc = result_stationary %>% 
	filter(name == "virus") %>%
	transmute(sa=sa, sv=sv, M2V=M2hat_mean, M3V=M3hat_mean, M4V=M4hat_mean)
result_stationary_con = inner_join(result_stationary_ac, result_stationary_vc, by = c("sa", "sv")) %>%
	mutate(region = 0)
result_stationary_2 = bind_rows(result_stationary_var, result_stationary_con)

# moment equations
result_stationary_p = result_stationary_2 %>%
	mutate(theta_t = theta_a + theta_v*region,
		   	ell = l*region + lhat*(1-region),
			E_M2 = (sa*M2A-region*sv*M2V)/theta_t/2,
			M2A_M3 = (4*theta_a+sa*M3A)/(1+4*theta_t),
			M2V_M3 = (4*theta_v-sv*M3V)/(1+4*theta_t),
			M2A_neutral = 4*theta_a/(1+4*theta_t),
			M2V_neutral = region*4*theta_v/(1+4*theta_t),
			E_M2_neutral = (sa*M2A_neutral-sv*M2V_neutral)/theta_t/2,
			M4A_theory = 24*2*5*theta_a^2/(3+28*theta_t),
			M4V_theory = 24*2*5*theta_v^2/(3+28*theta_t),
			M22A_theory = 8*2*6*theta_a^2/(3+28*theta_t),
			M22V_theory = 8*2*6*theta_v^2/(3+28*theta_t),
			M3A_theory = -(8/3*theta_a*E_M2_neutral/ell + (M22A_theory-M4A_theory/3)*sa)/(1+2*theta_t),
			M3V_theory = -(8/3*theta_v*E_M2_neutral/ell - (M22V_theory-M4V_theory/3)*sv)/(1+2*theta_t),
			M2A_theory = (4*theta_a+sa*M3A_theory)/(1+4*theta_t),
			M2V_theory = (4*theta_v-sv*M3V_theory)/(1+4*theta_t),
			E_theory = (sa*M2A_theory-region*sv*M2V_theory)/theta_t/2)

## FIGURE S1
qplot(data = result_stationary_p %>% filter(sv==0 | sv==1 | sv==2 | sv==3) %>%
	  	mutate(region = factor(region, levels=c(1,0))), x=sa, y=E, color=factor(sv), facets=.~region) + geom_line(aes(y=E_theory)) +
	geom_line(aes(y=E_M2), linetype=2)
ggsave("figures/stationary_E.pdf", height=3, width=9)

## FIGURE S2
# M2A
qplot(data = result_stationary_p %>% filter(sv==0 | sv==1 | sv==2 | sv==3, region == 1), x=sa, y=M2A, color=factor(sv)) +
	geom_line(aes(y=M2A_theory)) +
	geom_line(aes(y=M2A_M3), linetype=2) + guides(color="none")
ggsave("figures/stationary_M2A.pdf", height=3, width=4)
# M2Ahat
qplot(data = result_stationary_p %>% filter(sv==0 | sv==1 | sv==2 | sv==3, region == 0), x=sa, y=M2A, color=factor(sv)) +
	geom_line(aes(y=M2A_theory)) +
	geom_line(aes(y=M2A_M3), linetype=2) + guides(color="none")
ggsave("figures/stationary_M2Ahat.pdf", height=3, width=4)
# M2V
qplot(data = result_stationary_p %>% filter(sa==0 | sa==1 | sa==2 | sa==3, region == 1), x=sv, y=M2V, color=factor(sa)) +
	geom_line(aes(y=M2V_theory)) +
	geom_line(aes(y=M2V_M3), linetype=2) + guides(color="none")
ggsave("figures/stationary_M2V.pdf", height=3, width=4)
# E covariance
qplot(data = result_stationary %>% filter(sv==0 | sv==1 | sv==2 | sv==3, name == "antibody"), x=sa, y=EEcov_mean, color=factor(sv))  + guides(color="none") 
ggsave("figures/stationary_EEcov.pdf", height=3, width=4.2)
