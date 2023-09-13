# Doing FDR methods on YZ correctly
# for Panel on the p-hacking debate at FMA 2023
# Andrew Chen

# Environment ==================================================================
rm(list = ls())
library(data.table)
library(ggplot2)
library(tidyverse)
library(haven)
library(janitor)
library(gridExtra)

# read YZ sas file
yz0 <- read_sas("data-YZ/Yan_Zheng_RFS_Data.sas7bdat")

yz <- clean_names(yz0)
setDT(yz)

# combine transformation and fsvariable
yz[, signalname := paste0(transformation, "/", fsvariable)]

yz <- yz %>%
    pivot_longer(
        cols = c("ddiff_ew", "ddiff_vw"),
        names_to = "sweight",
        values_to = "ret",
        names_prefix = "ddiff_"
    ) %>%
    select(signalname, date, sweight, ret) %>%
    setDT()

yzsum <- yz %>%
    group_by(signalname, sweight) %>%
    summarize(
        rbar = mean(ret), vol = sd(ret), nmonth = n(),
        t = rbar / vol * sqrt(nmonth),
        .groups = "drop"
    ) %>%
    setDT()

# Storey FDR bounds ============================================================

lambda <- 0.8 # called theta in Harvey and Liu 2020, it seems

make_storey_dat <- function(lambda, stock_weight) {
    dat <- yzsum %>%
        filter(sweight == stock_weight) %>%
        mutate(
            t = abs(t),
            pval = 2 * pnorm(-t)
        )

    # Storey 2002 eqn 6
    m <- nrow(dat)
    pi0 <- sum(dat$pval > lambda) / (m * (1 - lambda))

    # loop over hurdles
    t_hurdles <- seq(0, 10, 0.01)
    storey <- tibble()
    for (i in 1:length(t_hurdles)) {
        temp <- tibble(
            t_hurdle = t_hurdles[i],
            Pr_signif = sum(dat$t > t_hurdles[i]) / m,
            p_val = 2 * pnorm(-t_hurdles[i]),
            fdrhat = pi0 * p_val / Pr_signif
        )
        storey <- bind_rows(storey, temp)
    }
    storey <- storey %>% filter(Pr_signif > 0)

    storey$lambda <- lambda
    storey$sweight <- stock_weight
    storey$pi0 <- pi0

    storey <- storey %>% select(lambda, sweight, pi0, everything())

    return(storey)
} # end make_storey_dat

storey <- rbind(
    make_storey_dat(lambda, "vw"),
    make_storey_dat(lambda, "ew")
)

# plot FDR vs hurdles and histogram ============================================
t_limits <- c(0, 6)

t_fdr_05 <- storey %>%
    filter(fdrhat < 0.05) %>%
    summarize(h = min(t_hurdle)) %>%
    pull(h)
pct_signals_signif <- sum(dat$t > t_fdr_05) / m * 100

fdr_t_gt_2 <- storey %>%
    filter(t_hurdle > 1.96) %>%
    filter(t_hurdle == min(t_hurdle)) %>%
    pull(fdrhat)

plt2 <- ggplot(dat, aes(x = t)) +
    geom_histogram(fill = "gray50", breaks = seq(0, 20, 0.5)) +
    geom_vline(xintercept = t_fdr_05, linetype = "twodash", color = "red") +
    labs(title = "Histogram of Yan-Zheng t-stats (raw VW return)", x = "|t-statistic|", y = "Count") +
    theme_bw() +
    annotate(
        "text",
        x = t_fdr_05 + 0.25, y = 3000, label = "5% signif hurdle", color = "red",
        angle = 90, size = 3
    ) +
    coord_cartesian(xlim = t_limits) +
    annotate("text",
        x = t_fdr_05 + 1, y = 500, label = paste0(round(pct_signals_signif, 1), "%"),
        color = "blue"
    ) +
    annotate("text",
        x = 1, y = 500, label = paste0(100 - round(pct_signals_signif, 1), "%"),
        color = "darkred"
    )

plt1 <- ggplot(storey, aes(x = t_hurdle, y = fdrhat)) +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    geom_vline(xintercept = t_fdr_05, linetype = "twodash", color = "red") +
    annotate(
        "text",
        x = t_fdr_05 + 0.25, y = 0.5, label = "5% signif hurdle", color = "red",
        angle = 90, size = 3
    ) +
    labs(title = "Storey 2002 FDR (theta=0.8)", x = "t-hurdle", y = "FDR") +
    theme_bw() +
    coord_cartesian(
        xlim = t_limits, ylim = c(0, 0.9)
    )

plt.both <- arrangeGrob(plt1, plt2, ncol = 1)
ggsave("temp.pdf", plt.both, scale = 0.5)
ggsave("yz-correctly.png", plt.both, scale = 0.5)


# plot histogram + null + signif hurdle ========================================

# reused plotting function
plot_dist <- function(stock_weight, t_edge) {
    Femp <- ecdf(abs(yzsum$t[yzsum$sweight == stock_weight]))
    m <- nrow(yzsum[yzsum$sweight == stock_weight])
    pi0 <- storey %>%
        filter(sweight == stock_weight) %>%
        distinct(pi0) %>%
        pull(pi0)
    t_fdr_05 <- storey %>%
        filter(sweight == stock_weight) %>%
        filter(fdrhat < 0.05) %>%
        summarize(h = min(t_hurdle)) %>%
        pull(h)

    histdat <- tibble(
        t_mid = (t_edge[1:(length(t_edge) - 1)] + t_edge[2:length(t_edge)]) / 2,
        data = diff(Femp(t_edge)) * m,
        null = diff(2 * pnorm(t_edge)) * m * pi0
    )

    plt <- ggplot(histdat, aes(x = t_mid)) +
        geom_bar(aes(y = data, fill = "data"), stat = "identity", width = 0.5) +
        geom_line(aes(y = null, group = 1, color = "null")) +
        scale_color_manual(" ", values = c("data" = "gray70", "null" = "red")) +
        scale_fill_manual("", values = "gray70") +
        theme_bw() +
        theme(
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.box = "vertical",
            legend.position = c(8, 7) / 10,
            legend.spacing.y = unit(-0, "cm")
        ) +
        labs(title = "Histogram of Yan-Zheng t-stats (raw VW return)", x = "|t-statistic|", y = "Num of Strategies") +
        coord_cartesian(xlim = c(0, 6)) +
        scale_y_continuous(breaks = seq(0, 8000, 500)) +
        scale_x_continuous(breaks = seq(0, 6, 1)) +
        geom_vline(xintercept = t_fdr_05, linetype = "twodash", color = "blue") +
        annotate(
            "text",
            x = t_fdr_05 + 0.1, y = max(histdat$data)*0.5, label = "5% signif hurdle", color = "blue",
            angle = 90, size = 3
        )

    dat = list(sweight = stock_weight
        , m=m
        , pi0 = pi0
        , t_fdr_05 = t_fdr_05
        , Pr_signif_05 = 1-Femp(t_fdr_05)
        , plt = plt)

    return(dat)
} # end plot_dist

t_edge <- seq(0, 10, 0.25)
dat_ew = plot_dist("ew", t_edge)
dat_vw = plot_dist("vw", t_edge)
plt_ew <- dat_ew$plt + labs(title = 'Raw EW return')
plt_vw <- dat_vw$plt + labs(title = 'Raw VW return')
plt_both <- arrangeGrob(plt_ew, plt_vw, ncol = 1)

ggsave("temp.pdf", plt_both, scale = 0.5)


# check pi0 as a function of lambda ===========================================



pval_dat = yzsum  %>% 
mutate(pval = 2 * pnorm(-abs(t)))  %>%
select(signalname,sweight,pval)  %>% 
pivot_wider(names_from = sweight, values_from = pval)   %>% 
filter(!is.na(ew+vw))

lambda_list <- seq(0.4, 0.95, 0.05)
pi0dat = tibble()
for (i in 1:length(lambda_list)) {
    temp = tibble(
        lambda = lambda_list[i]
        , pi0_ew = sum(pval_dat$ew > lambda_list[i]) / (nrow(pval_dat) * (1 - lambda_list[i]))
        , pi0_vw = sum(pval_dat$vw > lambda_list[i]) / (nrow(pval_dat) * (1 - lambda_list[i]))
    )
    pi0dat = rbind(pi0dat, temp)
}
pi0dat
i = 1
lambda_list[i]



# test 2 ======================================================================
h <- 0.25
denom <- 2 * (pnorm(h) - 0.5)
numer <- Femp(h)
numer
denom
numer / denom

lambda <- 0.5
numer2 <- sum(dat$pval > lambda) / m
denom2 <- (1 - lambda)
numer2
denom2
numer2 / denom2

histdat  %>% mutate(delta = data - null)  %>% 
    summarize(sum(delta)/m)
