library('dplyr'); library('readr');library('tidyr');library('ggplot2')
theme_set(theme_bw(12))

load('/data/jakub/figure_fixation.rdata')
#load('/data/jakub/figure_fixation_seffect.rdata')
l = 50
lhat = 50
# multiply values by the lineage frequency
output = output %>% mutate(
	E = E*frequency, Ehat = Ehat*frequency,
	M2 = M2*frequency, M2hat = M2hat*frequency,
	M3 = M3*frequency, M3hat = M3hat*frequency,
	LA1V1 = LA1V1*frequency)

output_pop = output %>%
	filter(t==0, name!='virus') %>% 
	group_by(theta_v, sa, sv, theta_a, bnab, id) %>% 
	summarize(Em = sum(E), Ehatm = sum(Ehat),
			  M2m = sum(M2), M2hatm = sum(M2hat),
			  M3m = sum(M3), M3hatm = sum(M3hat),
			  Q02 = sum(frequency*(name=='antibody resident' | !bnab)), Q02hat = sum(frequency*(name=='antibody invader' & bnab))) 
output_v = output %>% 
	filter(t==0, name == "virus") %>% 
	select(M2, sa, sv, theta_a, theta_v, bnab, id) %>%
	rename(M2V = M2)

output.mean = output %>% 
	filter(t==0, name == "antibody invader") %>% select(-name) %>%
	inner_join(output_pop, by = c("theta_v","sa","sv","theta_a","bnab", "id")) %>%
	inner_join(output_v, by = c("theta_v","sa","sv","theta_a","bnab", "id")) %>%
	mutate(LA12 = (E-Em*frequency)*M2m, 
		   LA12hat = (Ehat-Ehatm*frequency)*M2hatm,
		   Q02 = frequency*Q02, Q02hat = frequency*Q02hat,
		   Q2 = frequency*!bnab, Q2hat = frequency*bnab) %>%
	group_by(theta_v, sa, sv, theta_a, bnab) %>% 
	summarize_each(funs(mean))

output.pfix = output %>% 
	filter(t>0, name == "antibody invader") %>% 
	group_by(theta_v, sa, sv, theta_a, name, bnab) %>% 
	summarize(pfix = n()/100000) %>% 
	inner_join(output.mean, by=c("theta_v","sa","sv","theta_a","bnab")) %>%
	mutate(SA = sa/sqrt(l), 
		   L1.theory = SA*(E-Em*frequency)/(1+2*(theta_a+theta_v)) + SA*(Ehat-Ehatm*frequency)/(1+2*theta_a),
		   L2.neutral = SA^2*(4*theta_a*l*(Q2-Q02))/(1+2*(theta_a+theta_v))/(3+4*(theta_a+theta_v)) +
		   	SA^2*(4*theta_a*lhat*(Q2hat-Q02hat))/(1+2*(theta_a))/(3+4*(theta_a)),
		   L2.theory = SA^2*((M2-M2m*frequency )+ 4*theta_a*l*(Q2-Q02))/(1+2*(theta_a+theta_v))/(3+4*(theta_a+theta_v)) +
		   	SA^2*((M2hat-M2hatm*frequency) + 4*theta_a*lhat*(Q2hat-Q02hat))/(1+2*(theta_a))/(3+4*(theta_a)),
		   L3.theory = -SA^3*LA12/(1+2*(theta_a+theta_v))/(3+4*(theta_a+theta_v)) +
		   	-SA^3*LA12hat/(1+2*(theta_a))/(3+4*(theta_a)),
		   LAV = -SA*sv/sqrt(l)*(LA1V1-M2V*frequency)/(1+2*(theta_a+theta_v))/(1+4*(theta_a+theta_v)),
		   pfix.theory.1 = frequency + L1.theory + L2.theory + L3.theory + LAV,
		   pfix.theory.2 = frequency + L1.theory + L2.neutral + LAV)
save(output.pfix, file='figure_fixation.rdata')
qplot(data = output.pfix %>% filter(theta_a==1/50), x = M2V/l, y = pfix, color = bnab, facets=.~sa) +
	geom_line(data = output.pfix %>% filter(theta_a==1/50, sa<2), aes(y=pfix.theory.1)) +
	geom_line(aes(y=pfix.theory.2), linetype=2) + 
	guides(color="none")

ggsave('figures/pfix.pdf', height = 4, width = 8)
qplot(data = output.pfix, x = sv, y = pfix, color = bnab) + geom_line(aes(y=pfix.theory))
