# Code for all analyses for "No evidence for extrafloral nectaries in *Erythranthe angulosa*" Martin-Eberhardt, Smith, and Plunkert. 
# July 30th, 2026

# packages used
library(tidyverse)
library(glmmTMB)
library(performance)
library(DHARMa)


#### Data Cleaning ####

  # read in morphology dataset
d <- read.csv("Phenotype_morphology.csv")

# group individual maternal lines into morphs
d$phenotype <- case_when(d$family %in% c("11G", "12G", "50G", "51G") ~ "Robust Morph", 
                         d$family %in% c("53G", "54G", "77G", "78G") ~ "Thin-stemmed Morph", 
                         d$family == "1-1G2" ~ "Accession 1")
  # order levels for plotting
d$phenotype <- factor(d$phenotype, 
                      levels = c("Robust Morph", "Accession 1", "Thin-stemmed Morph"))


  # read in induction experiment dataset
ind <- read.csv("Induction_exp.csv")

  # assign accessions to morph categories
ind$phenotype <- case_when(ind$family %in% c("11G", "12G", "50G", "51G") ~ "Robust Morph", 
                           ind$family %in% c("53G", "54G", "77G", "78G") ~ "Thin-stemmed Morph", 
                           ind$family == "1-1G2" ~ "Accession 1")

  # order treatments
ind$treatment <- factor(ind$treatment, levels = c("ctrl", "soil", "JA", "clip"))

  # order morphs for plotting
ind$phenotype <- factor(ind$phenotype, levels = c("Robust Morph", "Accession 1", "Thin-stemmed Morph"))


#### Induction Experiment ####

# relative contribution of treatment and genetic family to explaining variation in number of axiliary swellings
mi <- glmmTMB(axillary_swellings ~ 1 + (1| treatment) + (1|family),  
              family = nbinom2(),
              data = ind)
summary(mi)
  # intraclass correlation
performance::icc(mi, by_group = TRUE)


# Phenotype predicts number of axillary swellings
m1 <- glmmTMB(axillary_swellings ~ phenotype, 
              family = nbinom2(),
              data = ind)
summary(m1)
check_model(m1)
simulateResiduals(m1) %>% plot()
emmeans::emmeans(m1, specs = "phenotype") %>% pairs()

# robust vs accession 1
# 59.6% fewer in axillary swellings in accession 1 compared to robust
mean <- -0.9066
se <- 0.2930
(1-exp(mean)) * 100
(1-exp(mean - se * 1.96)) * 100 
(1-exp(mean + se * 1.96)) * 100

# robust vs thin stemmed
# 83.7% fewer in axillary swellings in accession 1 compared to robust
mean <- -1.8124
se <- 0.1540
(1-exp(mean)) * 100
(1-exp(mean - se * 1.96)) * 100 
(1-exp(mean + se * 1.96)) * 100

# thin-stemmed vs accession1
mean <- 0.906 
se <- 0.289
exp(mean) * 100
exp(mean - se * 1.96) * 100 
exp(mean + se * 1.96) * 100


# Figure 2
ind %>% 
  mutate(Treatment = case_when(ind$treatment == "ctrl" ~ "Control", 
                               ind$treatment == "soil" ~ "Soil Amendment", 
                               ind$treatment == "JA" ~ "Jasmonic Acid", 
                               ind$treatment == "clip" ~ "Simluated Herbivory")
  ) %>% 
  ggplot(aes(x = phenotype, y = axillary_swellings)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_bar(stat = "summary", position = "dodge") + 
  # geom_errorbar(stat = "summary", position = "dodge") + 
  geom_jitter(position = position_jitterdodge(jitter.height = 0, 
                                              jitter.width = 0.2),
              aes(color = Treatment),  
              alpha = 0.7) + 
  labs(x = "", y = "Number of Axillary Swellings") + 
  theme_bw(base_size = 15)+ 
  theme(legend.position = c(.5, .8))
#ggsave("Figure2.png")


#### Morphology differences in two morphs ####

## stem thickness  
d$phenord <- factor(d$phenotype, levels = c("Thin-stemmed Morph", "Robust Morph"))
m1 <- glmmTMB(scale(log(stem_thickness_mm)) ~ phenord + (1|family), 
              data = d %>% filter(phenotype != "Accession 1"), 
              family = gaussian())
check_model(m1) 
simulateResiduals(m1) %>% plot() 
summary(m1)

# 95% CIs
2.44088 - 1.96 * 0.13118 
2.44088 + 1.96 * 0.13118 

## stem number 
m2 <- glmmTMB(stem_num ~ phenotype + (1|family), 
              data = d %>% filter(phenotype != "Accession 1"), 
              family = nbinom2()) # using neg. binomial b/c poisson overdispersed
check_overdispersion(m2)
check_model(m2) 
simulateResiduals(m2) %>% plot() 
summary(m2)

# estimate and 95% CIs
exp(0.9129)
exp(0.9129 - 1.96 * 0.1246   )
exp(0.9129 + 1.96 * 0.1246   )


## longest leaf length
m3 <- glmmTMB(scale(log(len_longest_leaf_mm)) ~ phenord + (1|family), 
              data = d %>% filter(phenotype != "Accession 1"), 
              family = gaussian()) 
check_model(m3) 
simulateResiduals(m3) %>% plot() 
summary(m3)

# 95% CIs
2.14846 - 1.96 * 0.20740  
2.14846 + 1.96 * 0.20740  



#### Supplemental Figures ####

# Figure S2: longest leaf length
d %>% 
  ggplot(aes(x = phenotype, y = len_longest_leaf_mm)) + 
  geom_boxplot(outlier.shape = NA) + # prevents double plotting of outliers when also plotting all data points in next layer
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  labs(x = "", y = "Longest Leaf Length (mm)") + 
  theme_bw()
#ggsave("FigS2_LongestLfLen.png", height = 4, width = 6.5)

# Figure S3: stem thickness
d %>% 
  ggplot(aes(x = phenotype, y = stem_thickness_mm)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  labs(x = "", y = "Stem Thickness (mm)") + 
  theme_bw()
#ggsave("FigS3_StemThickness.png", height = 4, width = 6.5)

# Figure S4: stem number
d %>% 
  ggplot(aes(x = phenotype, y = stem_num)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) + 
  labs(x = "", y= "Number of Stems") + 
  theme_bw()
#ggsave("FigS4_StemNumber.png", height = 4, width = 6.5)

#### Do plants with more branches just have more axillary swellings? ####

  # bind in phenotype data
ind <- ind %>% 
  left_join(d %>% select(!c("family", 
                            "phenotype", 
                            "treatment")), 
            by = c("plant_num")) 

  # negative bionomial b/c axillary swelling count previously shown to be overdispersed
m3 <- glmmTMB(axillary_swellings ~ stem_num, 
              family = nbinom2(), 
              data = ind)
summary(m3)
sjPlot::plot_model(m3, type = "pred", terms = "stem_num")

effect <- -0.11636
se <- 0.01729

# 11.0% decrease in axillary swellings for every additional stem
(1-exp(-0.11636)) * 100
(1-exp(-0.11636 - se * 1.96)) * 100 
(1-exp(-0.11636 + se * 1.96)) * 100

