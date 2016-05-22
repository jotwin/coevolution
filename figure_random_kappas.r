library('ggplot2'); library('dplyr'); library('tidyr');
theme_set(theme_bw(10))

### load data
load('figure_stationary_random_kappas.rdata')
#load('figure_random_kappas_unscaled.rdata')
theta_a = 1/50
theta_v = 1/50

datav = result_stationary %>% 
	filter(name == "antibody") %>% 
	mutate(E_theory = (sa*theta_a/(1+4*(theta_a+theta_v))-sv*theta_v/(1+4*(theta_a+theta_v)))/(theta_a+theta_v)*2,
		   E_M2 = E_mean/sqrt(M2_mean/theta_a/4*(1+4*(theta_a+theta_v)))) %>%
	#select(gamma_a, gamma_b, E_theory, E_M2) %>%
	mutate(region = "variable")
datac = result_stationary %>% 
	filter(name == "antibody") %>% 
	mutate(E_theory = 2*sa/(1+4*theta_a),
		   E_M2 = Ehat_mean/sqrt(M2hat_mean/theta_a/4*(1+4*(theta_a)))) %>%
	#select(gamma_a, gamma_b, E_theory, E_M2) %>%
	mutate(region = "conserved")
data = bind_rows(datav, datac)

qplot(data = data, x = E_theory, y = E_M2, shape = region, color = paste0(gamma_a, gamma_b)) + geom_abline() + coord_equal()

ggsave('figures/stationary_random_kappa.pdf', width=7, height=5)

