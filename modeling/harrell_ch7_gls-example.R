!#/usr/bin/R
  
# play around with GSL as alternative to M-M 
  
library(rms)
getHdata(cdystonia)
attach(cdystonia)

# make unique subject id 
uid <- with(cdystonia, factor(paste(site,id)))

# randomized drug/placebo; response variable is qttative score twstrs
# treatment at t=0, followup at t=2,4,12,16 weeks
# tabulate patterns of subjects' time points
# note: incomplete data points and censored obs 
table(tapply(week, uid, 
             function(w) paste(sort(unique(w)),collapse=' ')))
# plot raw data, superposing subjects
x1 <- xlab('Week'); y1 <- ylab('Outcome-total score')
ggplot(cdystonia, aes(x=week, y=twstrs, color = factor(id))) + 
  geom_line() + 
  x1 + y1 + facet_grid(treat ~ site) +
  guides(color = FALSE) 

# show quartiles, stratified by dose 
ggplot(cdystonia, 
       aes(x = week, y = twstrs)) + 
  x1 + y1 + 
  ylim(0, 70) + 
  stat_summary(fun.data = "median_hilow",
               conf.int = 0.5, geom = 'smooth') +
  facet_wrap(~ treat, nrow = 2)

# rearrange data so Yio is a baseline covariate (0th time point)
baseline <- subset(data.frame(cdystonia, uid), 
                   week == 0,
                   -week)
baseline <- upData(baseline, rename = c(twstrs= 'twstrs0'),
                   print = FALSE)
followup <- subset(data.frame(cdystonia, uid), week > 0,
                   c(uid, week, twstrs))
rm(uid)
both <- merge(baseline, followup, by = 'uid')
dd <- datadist(both)
options(datadist = 'dd')
  
# generalized least squares modeling
require(nlme)
# Six correlation patterns to try 
cp <- list(corCAR1, corExp, corCompSymm, corLin, corGaus, corSpher)

z <- vector('list', length(cp))
for (k in 1:length(cp)) {
  z[[k]] <- gls(twstrs ~ treat * rcs (week,3) + rcs(twstrs0, 3) + rcs(age, 4) * sex, data=both,
                correlation = cp[[k]](form = ~week | uid))
}
anova (z[[1]], z[[2]], z[[3]], z[[4]], z[[5]], z[[6]])
#corCAR1 looked best based on AIC 
a <- Gls(twstrs ~ treat * rcs (week,3) + rcs(twstrs0, 3) + rcs(age, 4) * sex, data=both,
    correlation = corCAR1(form = ~week | uid))
print(a)
# Todo: Example variogram

# Todo: Model assumptions
