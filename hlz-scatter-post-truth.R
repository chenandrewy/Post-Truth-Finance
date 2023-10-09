# Setup -------------------------------------------------------------------
rm(list = ls())

# uncomment for pretty fonts
# install.packages('extrafont')
# extrafont::font_import()

library(tidyverse)
library(data.table)
library(ggplot2)
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
    # , axis.title.x = element_text(size = 26),
    # axis.title.y = element_text(size = 26),
    # axis.text.x = element_text(size = 22),
    # axis.text.y = element_text(size = 22),
    # legend.text = element_text(size = 18),

    , axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    legend.text = element_text(size = 24),
    
    # Tweaking legend
    legend.position = c(.80, .15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 

# Simulate hlz ----------------------------------------------------------------

# true / false labels (for v)
label_true = 'True: [Exp Ret]>0'
label_false =  'False: [Exp Ret]=0'

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
mu_scatter = mu; mu_scatter[!v] = rnorm(sum(!v), 0, 0.05)*se

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
  Z, theta, v, theta_scatter
) %>% 
  mutate(
    t = theta + Z
    , tabs = abs(t)
    , mu = mu
    , mu_scatter = mu_scatter
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE)
      , labels = c(label_true, label_false)
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
    , nfalse = sum(v == label_false) 
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
hurdle_bonf05 = qnorm(1-0.05/300/2) # assumes everything is published, as in HLZ's conclusion text


# Plot setup ----------------------------------------------

# find a subset for plotting
nplot = 1500 # close to HLZ's estimate of "total factors" (ignoring unidentified scale)
set.seed(11)
small = dat[sample(1:n,nplot),]
small %>% filter(tabs > hurdle_01) %>% summarize(sum(v==label_false), n())


# plot sizing (shared by below)
texty = 250
textsize = 7
linesize = 1.1


# Plot slowly ----------------------------------------------------

p0 = ggplot(small %>% 
              filter(v==label_false), aes(x=tselect,y=mu_scatter)) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values=c(MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-0.jpg', width = 12, height = 8)

p1 = ggplot(small %>% 
              filter(v==label_false), aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values=c(MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-1.jpg', width = 12, height = 8)

p2 =  ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) + 
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)")) +
  theme(legend.position = 'none')

ggsave('../results/slow-2.jpg', width = 12, height = 8)

p3 = p2 + chen_theme

ggsave('../results/slow-3.jpg', width = 12, height = 8)

p4 = p3 +
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +  
  annotate(geom="text", 
           label="Discoveries ->", 
           x=3+1.1, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 0, size = textsize, color = 'black'
  ) 

ggsave('../results/slow-4.jpg', width = 12, height = 8)

p5 = p4 +  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) 

ggsave('../results/slow-5.jpg', width = 12, height = 8)


# plot FDR 5 --------------------------------------------------------------

p6 = p3 +
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) 

ggsave('../results/slow-6.jpg', width = 12, height = 8)
  

p7 = p6 +
  geom_vline(xintercept = hurdle_05, size = linesize, color = MATYELLOW, linetype = 'longdash') +  
  annotate(geom="text", 
           label="FDR = 5%", 
           x=hurdle_05, y=texty, vjust=-0.3, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  ) 

ggsave('../results/slow-7.jpg', width = 12, height = 8)

p8 = p7 + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )

ggsave('../results/slow-8.jpg', width = 12, height = 8)
  

# Plot big data -------------------------------------------------------

p1 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,300)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(0,500,50)) +  
  chen_theme +
  xlab(TeX("|t-statistic|")) +
  ylab(TeX("Expected Return (bps p.m.)"))

ggsave('../results/hlz-1.jpg', width = 12, height = 8)

# classical hurdle only
p5 = p1 + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )

ggsave('../results/hlz-classical-only.jpg', width = 12, height = 8)

# Plot post-truth -------------------------------------------------------

lab_hlz_false = '"False": Insignificant'
lab_hlz_true = '"True": Significant'

# post-truth labels
small = small %>% 
  mutate(
    v_hlz = case_when(
      tabs > hurdle_bonf05 ~ TRUE
      , tabs <= hurdle_bonf05 ~ FALSE
    )
  ) %>% mutate(
    v_hlz = factor(
      v_hlz, levels = c(TRUE,FALSE), labels = c(lab_hlz_true, lab_hlz_false)
    )
  ) %>% 
  mutate(
    v_hlz_alt = case_when(
      tabs > hurdle_05 ~ TRUE
      , tabs <= hurdle_05 ~ FALSE
    )
  ) %>% mutate(
    v_hlz_alt = factor(
      v_hlz_alt, levels = c(TRUE,FALSE), labels = c(lab_hlz_true, lab_hlz_false)
    )
  ) 


# plot classical only
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
  ylab(TeX("Expected Return (bps p.m.)")) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )

ggsave('../results/post-truth-1.jpg', width = 12, height = 8)

# add bonf
p2 = p1 + geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) 
ggsave('../results/post-truth-2.jpg', width = 12, height = 8)  
  

# relabel
p3 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v_hlz, color = v_hlz, shape = v_hlz), size = 2.5) +
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
  ylab(TeX("Expected Return (bps p.m.)")) + 
  geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )
 
ggsave('../results/post-truth-3.jpg', width = 12, height = 8)  

# add FDR 5%
p4 = p3  +
  geom_vline(xintercept = hurdle_05, size = linesize, color = MATYELLOW, linetype = 'longdash') +  
  annotate(geom="text", 
           label="FDR = 5%", 
           x=hurdle_05, y=texty, vjust=-0.3, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  )  
ggsave('../results/post-truth-4.jpg', width = 12, height = 8)  
  
# alternative relabelling
p5 = ggplot(small, aes(x=tselect,y=mu_scatter)) +
  geom_point(aes(group = v_hlz_alt, color = v_hlz_alt, shape = v_hlz_alt), size = 2.5) +
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
  ylab(TeX("Expected Return (bps p.m.)")) + geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  )  +
  geom_vline(xintercept = hurdle_05, size = linesize, color = MATYELLOW, linetype = 'longdash') +  
  annotate(geom="text", 
           label="FDR = 5%", 
           x=hurdle_05, y=texty, vjust=-0.3, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATYELLOW
  )  + 
  geom_vline(xintercept = hurdle_bonf05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Bonferroni 5\\%"), 
           x=hurdle_bonf05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) + 
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  )

ggsave('../results/post-truth-5.jpg', width = 12, height = 8)  

