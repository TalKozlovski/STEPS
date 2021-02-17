source('R\\help_functions.r')


### Create basic dataset
## Number of subjects in each group
n.subjects <- 200


levels.tab <- expand.grid("ID" = c(1:n.subjects), 
                          "motor_diff" = c(0, 1, 2), 
                          "cog_diff" = c(0, 1 , 2), 
                          "R" = c(1, 2))

ID_tab <- data.frame("ID" = c(1:n.subjects),  ## obs index
                     "severityM" = rnorm(n.subjects, 10, 1), 
                     "severityC" = rnorm(n.subjects, 10, 1))  %>% 
  mutate("severity" =  severityM * severityC) %>% 
  mutate("severity01" = percentize(severity))


sd.para <- 1 # Noise level to be added to explanatory variables
noise.ratio <- 100
px <- 5
rhos <- c(1) # learning rate for NLDSS
lambdaB <- 10
lambdaW <- 0.5
lambdaR <- 1
n.trees <- c(50)
scenario_tab <- expand.grid("Obs" = n.subjects,
                            "p" = px,
                            "Noise_level" = sd.para,
                            "Noise_ratio" = noise.ratio,
                            "lambdaB" = lambdaB,
                            "lambdaW" = lambdaW,
                            "lambdaR" = lambdaR, 
                            "Rho" = rhos, 
                            "ntrees" = n.trees)





iter.num <- 5
cors.results <- NULL # Correlations between performance and severity
conc.results <- NULL # Concordance values 
class.results <- NULL # classification abilities from planes. 
i <- 1
k <- 1
for(i in 1:nrow(scenario_tab))
  
  n.tmp                           <- scenario_tab[i, "Obs"]
p.tmp                           <- scenario_tab[i, "p"]
noise.level.tmp                 <- scenario_tab[i, "Noise_level"]
noise.ratio.tmp                 <- scenario_tab[i, "Noise_ratio"]
lambdaB.tmp                     <- scenario_tab[i, "lambdaB"]
lambdaW.tmp                     <- scenario_tab[i, "lambdaW"]
lambdaR.tmp                     <- scenario_tab[i, "lambdaR"]
rho.tmp                         <- scenario_tab[i, "Rho"]
ntrees.tmp                      <- scenario_tab[i, "ntrees"]

obs_cov <-  data.frame(ID_tab, 
                       "score1" = ID_tab$severityM, 
                       "score2" = ID_tab$severityC)
sim_dat.tmp <- x_sim_fun(ID_tab, levels.tab, pX = p.tmp)


# sim_dat = sim_dat.tmp
# obs_cov = obs_cov
# n = n.tmp
# p = p.tmp
# Rho = rho.tmp
# noise.level = noise.level.tmp
# noise.ratio = noise.ratio.tmp
# lambdaB = lambdaB.tmp
# lambdaW = lambdaW.tmp
# lambdaR = lambdaR.tmp
# results.tab = scenario_tab[i, , drop = F]
# iter.num = iter.num
# ntrees = ntrees.tmp


temp.result <- sim_fun(
  sim_dat = sim_dat.tmp,
  obs_cov = obs_cov,
  n = n.tmp,
  p = p.tmp,
  Rho = rho.tmp,
  noise.level = noise.level.tmp,
  noise.ratio = noise.ratio.tmp,
  lambdaB = lambdaB.tmp,
  lambdaW = lambdaW.tmp,
  lambdaR = lambdaR.tmp ,
  results.tab = scenario_tab[i, , drop = F],
  iter.num = iter.num,
  ntrees = ntrees.tmp
)




#### oBSERVED SEVERITY ILLUSTRATION 
# New facet label names for dose variable
motor.labs <- c("Difficulty domain #1: 0", "Difficulty domain #1: 1", "Difficulty domain #1: 2")
names(motor.labs) <- c("0", "1", "2")

# New facet label names for supp variable
cog.labs <- c("Difficulty domain #2: 0", "Difficulty domain #2: 1", "Difficulty domain #2: 2")
names(cog.labs) <- c("0", "1", "2")

g_severity <- sim_dat.tmp %>%
  ggplot(aes(severity01, obs_severity)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point() +
  facet_grid(cog_diff ~ motor_diff, labeller = labeller(cog_diff = cog.labs, motor_diff = motor.labs)) + 
  ylab("Observed severity udner stress") + 
  xlab("Simulated percentized severity") + 
  theme_classic() + 
  theme(text = element_text(face = "bold", size = 10)) + 
  theme_bw()






ggsave("fig\\simulation_obs_severity.png", plot = g_severity, width=10.5, height = 6, dpi=500)
ggsave("fig\\tiff\\simulation_obs_severity.tiff", plot = g_severity, dpi=500, device = "tiff")

