library('dplyr'); library('readr'); library('tidyr'); library('ggplot2')
theme_set(theme_bw(10))
load('figure_timecorrelation.rdata')

acf.to.df = function (dts, lag.max)
{
	d.acf = acf(dts, lag.max=lag.max, plot=FALSE)
	d = data_frame()
	for (i in 1:dim(d.acf$acf)[2])
		for (j in 1:dim(d.acf$acf)[3])
		{
			if (i<j) 
				d.name = paste(d.acf$snames[i], "&", d.acf$snames[j])
			else
				d.name = paste(d.acf$snames[j], "&", d.acf$snames[i])
			d.1 = data_frame(t = d.acf$lag[,i,j], variable = d.name, correlation = d.acf$acf[,i,j])
			d = bind_rows(d, d.1)
		}
	return (d)
}

d.antibody = output %>% filter(name=='antibody') %>% transmute(t=t, E=E, Ehat=Ehat, M2A=M2, M2Ahat=M2hat)
d.virus = output %>% filter(name=='virus') %>% transmute(t=t, M2V=M2)
d = inner_join(d.antibody, d.virus, by=c("t"))
d.acf = acf.to.df(ts(d %>% select(-t), frequency=1), 50)

qplot(data=d.acf %>% filter(variable == "E & E" | variable == "Ehat & Ehat" | variable == "M2A & M2A" |variable == "M2Ahat & M2Ahat" | variable == "M2V & M2V") , x = t, y=correlation, geom="line", color = variable) + 
	geom_line(data=d.acf %>% filter(variable=="E & E"), aes(y=exp(-2*(1/25*t))), linetype=2) +
	geom_line(data=d.acf %>% filter(variable=="Ehat & Ehat"), aes(y=exp(-2*(1/50*t))), linetype=2) +
	geom_line(data=d.acf %>% filter(variable=="E & E"), aes(y=exp(-(t))), linetype=2, color="black")
ggsave('figures/stationary_timecorrelation.pdf', height=4, width=7)
