library(tidyverse)
df <- readr::read_csv("data/sample_data.csv") %>%
  mutate(
    sex = factor(sex, levels=c("F","M","U")),
    log_hg_blood   = log(hg_blood_ug_g),
    log_hg_feather = log(hg_feather_ug_g)
  )
df %>% summarise(across(c(hg_blood_ug_g,hg_feather_ug_g,mass_g,wing_length_mm,tarsus_length_mm),
                        list(n=~sum(!is.na(.)), mean=mean, sd=sd, min=min, max=max), .names="{.col}.{.fn}"))
cor.test(df$log_hg_blood, df$log_hg_feather, method="pearson")
cor.test(df$hg_blood_ug_g, df$mass_g, method="spearman")
cor.test(df$hg_feather_ug_g, df$wing_length_mm, method="spearman")

anova_blood  <- aov(log_hg_blood ~ sex, data=df);  summary(anova_blood)
anova_feath  <- aov(log_hg_feather ~ sex, data=df); summary(anova_feath)
TukeyHSD(anova_blood); TukeyHSD(anova_feath)

kruskal.test(hg_blood_ug_g ~ sex, data=df)

m1 <- lm(log_hg_blood ~ sex + mass_g + wing_length_mm + tarsus_length_mm, data=df)
m2 <- lm(log_hg_feather ~ sex + mass_g + wing_length_mm + tarsus_length_mm, data=df)
summary(m1); summary(m2)

par(mfrow=c(1,2)); plot(m1); plot(m2); par(mfrow=c(1,1))

m_krew <- lm(log_hg_blood ~ log_hg_feather + sex + mass_g + wing_length_mm + tarsus_length_mm, data=df)
summary(m_krew)
library(MASS)
rlm_blood <- rlm(log_hg_blood ~ log_hg_feather + sex + mass_g, data=df); summary(rlm_blood)

library(quantreg)
rq_blood  <- rq(log_hg_blood ~ log_hg_feather + sex + mass_g, data=df, tau=0.5)
summary(rq_blood)

vars <- df %>% select(log_hg_blood, log_hg_feather, mass_g, wing_length_mm, tarsus_length_mm) %>% drop_na()
pca <- prcomp(scale(vars), center=TRUE, scale.=TRUE)
summary(pca); pca$rotation

library(ggplot2)
ggplot(df, aes(sex, log_hg_blood)) + geom_boxplot() + geom_jitter(width=0.1, alpha=.6)
ggplot(df, aes(log_hg_feather, log_hg_blood)) + geom_point() + geom_smooth(method="lm", se=TRUE)

