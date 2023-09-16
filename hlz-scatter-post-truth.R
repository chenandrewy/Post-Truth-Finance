# Setup -------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(data.table)
library(ggplot2)
library(zoo)
library(extrafont)
library(latex2exp)

dir.create('../results/')


MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)


chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    
    # Font sizes
    , axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 18),
    
    # Tweaking legend
    legend.position = c(0.7, 0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 


# Simulate hlz ----------------------------------------------------------------

nport = 1e4
nsim = 500
n = nport*nsim

set.seed(121)

# hlz baseline
v  = runif(n) > 0.444 # Table 5, rho = 0.2, p_0
se = 1500/sqrt(12)/sqrt(240) # page 29, N = 240 and page 30, annual volatilty of 15% 
mu = rexp(n, 1/55.5); mu[!v] = 0 # Table 5, rho = 0.2, \lambda
theta = mu / se
theta_scatter = theta; theta_scatter[!v] = rnorm(sum(!v), 0, 0.1)
mu_scatter = mu; mu_scatter[!v] = rnorm(sum(!v), 0, 0.1)*se
pubnoise = runif(n)

# Fabian Winkler's fast method for common correlations
rho = 0.2
sigc = sqrt(rho)
sige = sqrt(1-sigc^2)

# simulate common component (in blocks)
c = matrix(rnorm(nsim, 0, sigc), nrow = nsim, ncol = nport) %>% t() %>% 
  as.vector
# add idiosyncratic noise
e = rnorm(nport*nsim, 0, sige)
Z = c + e

# sanity check
# c2 = rnorm(n, 0, sigc)
# z1 = c2 + rnorm(n, 0, sige)
# z2 = c2 + rnorm(n, 0, sige)
# cor(z2,z1)
# var(z1)
# var(z2)

# assemble into data table
dat = data.table(
  Z, theta, v, pubnoise, theta_scatter
) %>% 
  mutate(
    t = theta + Z
    , tabs = abs(t)
    , pub = case_when(
      tabs <= 1.96 ~ FALSE
      , (1.96 < tabs) & (tabs <= 2.57) ~ pubnoise < 0.5
      , 2.57 < tabs  ~ TRUE
    )
    , mu = mu
    , mu_scatter = mu_scatter
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE), labels = c('True Predictor', 'False Predictor')
    )
  ) %>%   
  mutate(
    tselect = tabs
  )

# fit shrinkage
datsum = dat %>% 
  mutate(
    tgroup = ntile(tselect,1000)
  ) %>% 
  group_by(tgroup) %>% 
  summarize(
    tselect_left = min(tselect), tselect_right = max(tselect)
    , Etselect = mean(tselect), Etheta = mean(theta), n = dplyr::n()
    , nfalse = sum(v == 'False Predictor') 
  ) 


# fit "cdf"
datsum = datsum %>% 
  arrange(-Etselect) %>% 
  mutate(
    nfalse_cum = cumsum(nfalse)
    , n_cum = cumsum(n)
    , fdr_tselect_left = nfalse_cum / n_cum * 100
  ) %>% 
  arrange(Etselect)

# find hurdles
hurdle_05 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 5)])
hurdle_01 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 1)])
hurdle_bonf05 = qnorm(1-0.05/300/2)


# find pub stats
pubstat = dat %>% 
  filter(tselect>2) %>% 
  summarize(
    Et = mean(tselect)
    , Etheta = mean(theta)
    , FDR = mean(theta <= 0)    
    , shrink = 1-mean(theta) / mean(tselect)
  )

pubplot = tibble(
  tselect = seq(2,10,0.1), theta = pubstat$Et, group = 'naive' 
)  %>% rbind(
  tibble(
    tselect = seq(2,10,0.1), theta = pubstat$Etheta, group = 'shrinkage' 
  )   
)




# Plot scatter bonf -------------------------------------------------------
# for ease of presentation
ntotal = 300/mean(dat$pub)
set.seed(430)
bonf_05 = qnorm(1-0.05/2/ntotal)

# settings for both panels here
nplot = 1500
texty = 250
textsize = 7

linesize = 1.1

set.seed(11)
small = dat[sample(1:n,nplot),]

small %>% filter(tabs > hurdle_01) %>% summarize(sum(v=='False Predictor'))

# holm algo
holm_05 = small %>% select(tabs) %>% 
  arrange(desc(tabs)) %>% 
  mutate(
    pval = 2*pnorm(-tabs)
    , k = row_number()
    , signif = pval < 0.05/(ntotal + 1 - k)
  ) %>% 
  filter(signif == F) %>% 
  filter(row_number() == 1) %>% 
  pull(tabs)



p1 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  theme(
    legend.position = c(.80, .15)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)"))

# HURDLES
p2 = p1 +
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) +  
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) +
  geom_vline(xintercept = bonf_05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=holm_05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) 


ggsave('../results/hlz-scatter-bonf.pdf', p2, width = 12, height = 8, device = cairo_pdf)




# Plot post-truth ---------------------------------------------------------

# for ease of presentation
ntotal = 300/mean(dat$pub)
set.seed(430)
bonf_05 = qnorm(1-0.05/2/ntotal)

# settings for both panels here
nplot = 1500
texty = 250
textsize = 7

linesize = 1.1

set.seed(11)
small = dat[sample(1:n,nplot),]

# post-truth labels
small = small %>% 
  mutate(
    v = case_when(
      tabs > bonf_05 ~ TRUE
      , tabs <= bonf_05 ~ FALSE
    )
  ) %>% mutate(
    v = factor(
      v, levels = c(TRUE,FALSE), labels = c('True Predictor', 'False Predictor')
    )
  ) 


# holm algo
holm_05 = small %>% select(tabs) %>% 
  arrange(desc(tabs)) %>% 
  mutate(
    pval = 2*pnorm(-tabs)
    , k = row_number()
    , signif = pval < 0.05/(ntotal + 1 - k)
  ) %>% 
  filter(signif == F) %>% 
  filter(row_number() == 1) %>% 
  pull(tabs)



p1 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  theme(
    legend.position = c(.80, .15)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)"))

# HURDLES
p2 = p1 +
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) +  
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) +
  geom_vline(xintercept = bonf_05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=holm_05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) 


ggsave('../results/hlz-post-truth.pdf', p2, width = 12, height = 8, device = cairo_pdf)



