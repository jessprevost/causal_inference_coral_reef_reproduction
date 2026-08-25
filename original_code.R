# Install and required packages 
# Note that the rethinking package requires the installation of RStan, with instructions varying depending on the computer used. Detailed instructions for download can be found here: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 

install.packages(“rethinking”) 
library(rethinking) 
install.packages(“rstanarm”)
library(rstanarm)
install.packages(“ggeffects”)
library(ggeffects)
install.packages(“ggplot2”)
library(ggplot2)
install.packages(“dagitty”)
library(dagitty)

# Test DAG-data consistency (see Supplementary Material for detail) 

# Upload data (see attached data file, ‘SeychellesS’)
SeychellesS

#download DAG from dagitty.net

DAG <- downloadGraph("dagitty.net/mbt1AvL") 

#evaluate the d-separation implications of the DAG

r <- localTests(DAG,SeychellesS)

#perform Holm-Bonferroni correction 

r$p.value <- p.adjust(r$p.value)

r # should show all 1s under p-value (if DAG-data consistency is ensured)

# Run a Bayesian logistic regression model (with weakly informed priors) to determine the causal effect of each predictor variable of interest on regime shift trajectory. Apply the backdoor criterion for model selection. Use ‘SeychellesS’ dataset, where predictor variables are standardized.  

# Depth model 

causal_depth <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bd*Depth,
    a ~ cauchy(0,10),
    bd ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for depth  
post <- extract.samples(causal_depth)
post_depth <- post$bd
PI(post_depth, prob = 0.90)
PI(post_depth, prob = 0.50) 
median(post_depth) 

# MPA model 

causal_MPA <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bm*MPA,
    a ~ cauchy(0,10),
    bm ~ dstudent(4,0,2.5),
    bbc ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for MPA  
post <- extract.samples(causal_MPA)
post_MPA <- post$bm
PI(post_MPA, prob = 0.90)
PI(post_MPA, prob = 0.50) 
median(post_MPA) 

# Nutrient model 

causal_nutrient <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bn*Nutrient + bd*Depth,
    a ~ cauchy(0,10) ,
    bn ~ dstudent(4,0,2.5) ,
    bd ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for nutrient  
post <- extract.samples(causal_nutrient)
post_nutrient <- post$bn
PI(post_nutrient, prob = 0.90)
PI(post_nutrient, prob = 0.50) 
median(post_nutrient) 

# Branching coral model 

causal_branching <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bbc*Branching + bd*Depth + bm*Manage + bw*Wave,
    a ~ cauchy(0,10) ,
    bbc ~ dstudent(4,0,2.5),
    bd ~ dstudent(4,0,2.5),
    bm ~ dstudent(4,0,2.5),
    bw ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for branching coral   
post <- extract.samples(causal_branching)
post_branching <- post$bbc
PI(post_branching, prob = 0.90)
PI(post_branching, prob = 0.50) 
median(post_branching) 
# Structural complexity model

causal_struccomplex <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bsc*Struccomplex + bbc*Branching,
    a ~ cauchy(0,10) ,
    bsc ~ dstudent(4,0,2.5),
    bbc ~ dstudent(4,0,2.5),
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for structural complexity 
post <- extract.samples(causal_struccomplex)
post_struccomplex <- post$bsc
PI(post_struccomplex, prob = 0.90)
PI(post_struccomplex, prob = 0.50) 
median(post_struccomplex) 

# Wave model 

causal_wave <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bw*Wave,
    a ~ cauchy(0,10) ,
    bw ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for wave 
post <- extract.samples(causal_wave)
post_wave <- post$bw
PI(post_wave, prob = 0.90)
PI(post_wave, prob = 0.50) 
median(post_wave) 

# Herbivore biomass model 

causal_herb <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bh*Herb + bs*Struccomplex + bb*Branching,
    a ~ cauchy(0,10) ,
    bh ~ dstudent(4,0,2.5),
    bsc ~ dstudent(4,0,2.5),
    bbc ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for herbivore biomass 
post <- extract.samples(causal_herb)
post_herb <- post$bh
PI(post_herb, prob = 0.90)
PI(post_herb, prob = 0.50) 
median(post_herb) 

# Initial macroalgae cover model 

causal_macro <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bm*Macro + bh*Herb + bd*Depth + bn*Nutrients + bw*Wave,
    a ~ cauchy(0,10) ,
    bm ~ dstudent(4,0,2.5),
    bh ~ dstudent(4,0,2.5),
    bd ~ dstudent(4,0,2.5),
    bn ~ dstudent(4,0,2.5),
    bw ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% and 50% percentile interval (PI) for macroalgae 
post <- extract.samples(causal_macro)
post_macro <- post$bm
PI(post_macro, prob = 0.90)
PI(post_macro, prob = 0.50) 
median(post_macro) 

# Use median and 50% and 90% PI values from above models to create Manuscript Figure 3a. 

label <- c("MPA", "Carbon:Nitrogen", "Herbivore Biomass", "Branching Coral", "Depth", "Structural Complexity", "Wave Exposure", "Initial Macroalgae")
median <- c(0.06, -1.61, 0.96, 0.84, -2.12, -3.77, 1.02, 3.3)

lower90 <- c(-1.20, -3.65, -0.87, -1.21, -4.19, -7.33, -0.45, 0.40)
upper90 <- c(1.39, 0.05, 3.35, 3.09, -0.49, -1.63, 2.54, 10.14)

lower50 <- c(-0.46, -2.39, 0.14, 0.02, -2.89, -5.01, 0.38, 1.80)
upper50 <- c(0.63, -0.90, 1.87, 1.74, -1.38, -2.69, 1.63, 5.29)

df <- data.frame(label, median, lower90, upper90, lower50, upper50)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label))

library(ggplot2)
ggplot(data=df, aes(x=label, y=median, ymin=lower90, ymax=upper90)) +
  geom_errorbar(aes(ymin=lower50, ymax=upper50), width =0, size = 1.3)+ #I just added a new error bar on top of your original error bar with the dummy data that I made (lower50, upper50). width = 0 gets ride of the whiskers, size controls the thickness. you can also do color = "put your favoriate color" if you wan to use color 
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab(" ") + ylab("Standardized Effect Size") +
  theme_bw()  # use a white background

# Plot marginal plots of continuous predictor variables as shown in Figure 3b and 4b. Use ‘Seychelles’ dataset, which has unstandardized data.

# Depth marginal plot 

fit <- stan_glm(Shifted ~ Depth, data = Seychelles, family = binomial("logit"),chains = 1)
depth_MP <- ggpredict(fit)
plot(depth_MP, add.data = TRUE, jitter = 0)

# Nutrient marginal plot 

fit <- stan_glm(Shifted ~ Nutrients + Depth, data = Seychelles, family = binomial("logit"), chains = 1)
nutrient_MP <- ggpredict(fit)
plot(nutrient_MP, add.data = TRUE, jitter = 0.0)

# Branching coral marginal plot 

fit <- stan_glm(Shifted ~ Branching + Depth + MPA + Wave, data = Seychelles, family = binomial("logit"), chains = 1)
branching_MP <- ggpredict(fit)
plot(branching_MP, add.data = TRUE, jitter = 0.0)

# Structural complexity marginal plot 

fit <- stan_glm(Shifted ~ Struccomplex + Branching, data = Seychelles, family = binomial("logit"),chains = 1)
struccomplex_MP <- ggpredict(fit) 
plot(struccomplex_MP, add.data = TRUE, jitter = 0.0)

# Wave exposure marginal plot 

fit <- stan_glm(Shifted ~ Wave, data = Seychelles, family = binomial("logit"), chains = 1)
wave_MP <- ggpredict(fit)
plot(wave_MP, add.data = TRUE, jitter = 0.0)

# Herbivore biomass marginal plot 

fit <- stan_glm(Shifted ~ Herb + Struccomplex + Branching, data = Seychelles, family = binomial("logit"), chains = 1)
herb_MP <- ggpredict(fit)
plot(herb_MP, add.data = TRUE, jitter = 0.0)

# Initial macroalgae cover marginal plot 

fit <- stan_glm(Shifted ~ Macro + Herb + Depth + Nutrients + Wave, data = Seychelles, family = binomial("logit"), chains = 1)
macro_MP <- ggpredict(fit, terms = "Macroalgae[0:20 by=0.1]")
plot(macro_MP, add.data = TRUE, jitter = 0.0)

# Predictive model for predicting which of the 12 recovering reefs are expected to regime shift following the 2016 bleaching event 

# Train model using ‘Seychelles’ dataset 

SeyPre <- map(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a  + bd*Depth + bn*Nutrients + bw*Wave + bic*Struccomplex  + bbc*Branching + bim*Macro,
    a ~ dnorm(0,10),
    bd ~ dnorm(0,10),
    bn ~ dnorm(0,10),
    bw ~ dnorm(0,10),
    bic ~ dnorm(0,10),
    bbc ~ dnorm(0,30),
    bim ~ dnorm(0,10)
  ), 
  data = Seychelles)

# Posterior Validation Check

postcheck(SeyPre)

# Predict using 2014 observational dataset “Seychelles2014” 

predict <- link(SeyPre, data=Seychelles2014)
str(predict)
summary(predict)

# The dataset “Predictions” show the individual probabilities of regime shifting estimated for each recovering reef site. 

# Plot predictions as shown in Figure 4a

ggplot(Predictions, aes(x = Probability.of.Regime.Shift, y = 	Reef.Site)) + geom_density_ridges(aes(fill = Reef.Site)) + 
  scale_fill_manual(values = c("yellowgreen", "yellow", 	"orange", "red2", "deeppink", "purple", "lightsteelblue3", 	"darkblue", "slateblue4", "blue", "deepskyblue3", 	"cadetblue2"))

# Supplementary Materials Table 2: Causal salad and no-covariate models using ‘SeychellesS’ dataset 

# Causal salad model 

causal_salad <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bm*MPA + bn*Nutrients + bh*Herb + bb*Branching + bd*Depth + bs*Struccomplex + bw*Wave + bim*Macro,
    a ~ cauchy(0,10) ,
    bm ~ dstudent(4,0,2.5),
    bn ~ dstudent(4,0,2.5),
    bh ~ dstudent(4,0,2.5),
    bb ~ dstudent(4,0,2.5),
    bd ~ dstudent(4,0,2.5),
    bs ~ dstudent(4,0,2.5),
    bw ~ dstudent(4,0,2.5),
    bim ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS) 

# show median and 90% percentile interval (PI) for predictor variables   
post <- extract.samples(causal_salad)
post_MPA <- post$bm
post_nutrients <- post$bn
post_herb <- post$bh
post_branching <- post$bb
post_depth <- post$bd
post_struccomplex <- post$bs
post_wave <- post$bw
post_macro <- post$bim

median(post_MPA)
median(post_nutrients)
median(post_herb)
median(post_branching)
median(post_depth)
median(post_struccomplex)
median(post_wave)
median(post_macro)

PI(post_MPA, prob = 0.90)
PI(post_nutrients, prob = 0.90)
PI(post_herb, prob = 0.90)
PI(post_branching, prob = 0.90)
PI(post_depth, prob = 0.90)
PI(post_struccomplex, prob = 0.90)
PI(post_wave, prob = 0.90)
PI(post_macro, prob = 0.90)

# No-covariate models 

# Nutrient 

nutrient_nc <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bn*Nutrient,
    a ~ cauchy(0,10) ,
    bn ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% percentile interval (PI) for nutrient  
post <- extract.samples(nutrient_nc)
post_nutrient <- post$bn
PI(post_nutrient, prob = 0.90)
median(post_nutrient) 

# Branching coral 

Branching_nc <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bbc*Branching,
    a ~ cauchy(0,10) ,
    bbc ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% percentile interval (PI) for branching coral   
post <- extract.samples(branching_nc)
post_branching <- post$bbc
PI(post_branching, prob = 0.90)
median(post_branching) 

# Structural complexity 

Struccomplex_nc <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bsc*Struccomplex,
    a ~ cauchy(0,10) ,
    bsc ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% percentile interval (PI) for structural complexity 
post <- extract.samples(struccomplex_nc)
post_struccomplex <- post$bsc
PI(post_struccomplex, prob = 0.90)
median(post_struccomplex) 

# Herbivore biomass 

Herb_nc <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bh*Herb,
    a ~ cauchy(0,10) ,
    bh ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% percentile interval (PI) for herbivore biomass 
post <- extract.samples(herb_nc)
post_herb <- post$bh
PI(post_herb, prob = 0.90)
median(post_herb) 

# Initial macroalgae cover 

macro_nc <- map2stan(
  alist(
    Shifted ~ dbinom( 1 , p ) ,
    logit(p) <- a + bm*Macro,
    a ~ cauchy(0,10) ,
    bm ~ dstudent(4,0,2.5)
  ), 
  data = SeychellesS)

# show median and 90% percentile interval (PI) for macroalgae 
post <- extract.samples(macro_nc)
post_macro <- post$bm
PI(post_macro, prob = 0.90)
median(post_macro) 



