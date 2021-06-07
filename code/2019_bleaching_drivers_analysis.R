#### hawaii coral bleaching analysis: investigating drivers ####
## written by: morgan winston

## this script uses cluster-level data described & linked to here: {InPort record}
## code performs the following: identifies best-fit model and creates series of plots to visualize model output and predictions

# initialization ####
# set working directory & load data
setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/For InPort/Environmental Drivers")
hcbc <- read.csv("HCBC_2019_ClusteredData.csv")

# load packages
library(ggmap) # for plotting
library(car) # for our variance inflation test
library(MASS) # for stepAIC()
library(plotrix) # for std.error

# load googlemaps API Key
ggmapAPI = readChar("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/Tom's Mapping Example/morgan.winston_googlemapsAPIkey.txt",nchars = 39)
register_google(key = ggmapAPI)

# load useful functions
g_legend<-function(a.gplot){ # function to extract legend from a plot
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
'%!in%' <- function(x,y)!('%in%'(x,y)) # the opposite of %in%

# subset data ####
hcbc_complete <- hcbc[complete.cases(hcbc[c("kdPAR_mn", "SiteSuscp_mn","PAR_surface_mn","PctBleached_hist_mn","WaveAction_value_mn")]),]
hcbc_complete <- hcbc_complete[ which(hcbc_complete$Depth_ft < 60),] # no deep surveys included

# standardize & center drivers to put all on the same scale ####
# if values include 0, sqrt trans
# if values don't include 0, log trans
hcbc_complete$DHW.MeanMax.YR01_mn_s <- as.vector(scale(hcbc_complete$DHW.MeanMax.YR01_mn))
hcbc_complete$DHW.MeanMax.YR10YR01_mn_s <- as.vector(scale(hcbc_complete$DHW.MeanMax.YR10YR01_mn))
hcbc_complete$Depth_ft_mn_s <- as.vector(scale(hcbc_complete$Depth_ft_mn))
hcbc_complete$PctBleached_hist_mn_s <- as.vector(scale(hcbc_complete$PctBleached_hist_mn))
hcbc_complete$SiteSuscp_mn_log_s <- as.vector(scale(log(hcbc_complete$SiteSuscp_mn))) # log trans
hcbc_complete$PAR_surface_mn_s <- as.vector(scale(hcbc_complete$PAR_surface_mn))
hcbc_complete$WaveAction_value_mn_sqrt_s <- as.vector(scale(sqrt(hcbc_complete$WaveAction_value_mn))) # sqrt trans
hcbc_complete$TotalEffluent_mn_sqrt_s <- as.vector(scale(sqrt(hcbc_complete$TotalEffluent_mn))) # sqrt trans
hcbc_complete$TourRec_10yrAvgPUD_mn_log_s <- as.vector(scale(log(hcbc_complete$TourRec_10yrAvgPUD_mn))) # sqrt trans
hcbc_complete$AgGolf_runoff_mn_sqrt_s <- as.vector(scale(sqrt(hcbc_complete$AgGolf_runoff_mn))) # sqrt
hcbc_complete$Urban_runoff_mn_sqrt_s <- as.vector(scale(sqrt(hcbc_complete$Urban_runoff_mn))) # sqrt
hcbc_complete$kdPAR_mn_s <- as.vector(scale(hcbc_complete$kdPAR_mn))
hcbc_complete$SST_Variability_AllB4_mn_s <- as.vector(scale(hcbc_complete$SST_Variability_AllB4_mn))

# check for correlations ####
full_pred <- hcbc_complete[,c("Depth_ft_mn_s",
                              "SiteSuscp_mn_log_s",
                              "DHW.MeanMax.YR01_mn_s",
                              "DHW.MeanMax.YR10YR01_mn_s",
                              "PctBleached_hist_mn_s",
                              "SST_Variability_AllB4_mn_s",
                              "PAR_surface_mn_s",
                              "kdPAR_mn_s",
                              "WaveAction_value_mn_sqrt_s",
                              "TotalEffluent_mn_sqrt_s",
                              "TourRec_10yrAvgPUD_mn_log_s",
                              "AgGolf_runoff_mn_sqrt_s",
                              "Urban_runoff_mn_sqrt_s")]
pred_cor <- cor(full_pred)
head(round(pred_cor,2))
corrplot::corrplot(pred_cor, method="circle") # looks pretty good (why? no coefficients > 0.6 or < -0.6)
pairs(full_pred)
full_mod_corrplot <- ggcorrplot(pred_cor,
                                 hc.order = TRUE,
                                 type = "lower",
                                 lab = TRUE)
full_mod_pairsplot <- ggpairs(full_pred)

### ### ### ### ### ### ### ### ### ### ### 
# model selection -- using stepAIC() ####

# sqrt transform response variable
hcbc_complete$PctBleached_mean_sqrt <- sqrt(hcbc_complete$CoralBleached_Perc_mn)

# run model first with no interactions to assess VIF
mod_full <- lm(PctBleached_mean_sqrt ~
                     DHW.MeanMax.YR01_mn_s +
                     DHW.MeanMax.YR10YR01_mn_s +
                     Depth_ft_mn_s +
                     PctBleached_hist_mn_s +
                     SiteSuscp_mn_log_s +
                     SST_Variability_AllB4_mn_s +
                     kdPAR_mn_s +
                     PAR_surface_mn_s +
                     WaveAction_value_mn_sqrt_s +
                     TotalEffluent_mn_sqrt_s +
                     TourRec_10yrAvgPUD_mn_log_s +
                     AgGolf_runoff_mn_sqrt_s +
                     Urban_runoff_mn_sqrt_s,
                hcbc_complete,
                weights = Weights_DriversAndSpatialAnalysis)
step_mod_bic <- stepAIC(mod_full, k = log(nrow(hcbc_complete)))
summary(step_mod_bic)

# variance inflaction factors #### 
# checks for collinearity
vif(step_mod_bic) # looks good, all < 2 (if > 4 then that is bad)

# now testing interactions that we have hypotheses about 
mod_full_2 <- lm(PctBleached_mean_sqrt ~ 
                   DHW.MeanMax.YR01_mn_s +
                   DHW.MeanMax.YR10YR01_mn_s + 
                   SST_Variability_AllB4_mn_s +
                   Depth_ft_mn_s + 
                   PctBleached_hist_mn_s + 
                   SiteSuscp_mn_log_s +
                   PAR_surface_mn_s + 
                   kdPAR_mn_s +
                   WaveAction_value_mn_sqrt_s +
                   TotalEffluent_mn_sqrt_s +
                   TourRec_10yrAvgPUD_mn_log_s +
                   AgGolf_runoff_mn_sqrt_s +
                   Urban_runoff_mn_sqrt_s + 
                   DHW.MeanMax.YR01_mn_s:kdPAR_mn_s + ### test interactions b/w DHW01 and every possibility
                   DHW.MeanMax.YR01_mn_s:PAR_surface_mn_s + 
                   DHW.MeanMax.YR01_mn_s:SiteSuscp_mn_log_s +
                   DHW.MeanMax.YR01_mn_s:DHW.MeanMax.YR10YR01_mn_s +
                   DHW.MeanMax.YR01_mn_s:SST_Variability_AllB4_mn_s +
                   DHW.MeanMax.YR01_mn_s:PctBleached_hist_mn_s +
                   DHW.MeanMax.YR01_mn_s:WaveAction_value_mn_sqrt_s +
                   DHW.MeanMax.YR01_mn_s:TotalEffluent_mn_sqrt_s +
                   DHW.MeanMax.YR01_mn_s:TourRec_10yrAvgPUD_mn_log_s +
                   DHW.MeanMax.YR01_mn_s:AgGolf_runoff_mn_sqrt_s +
                   DHW.MeanMax.YR01_mn_s:Urban_runoff_mn_sqrt_s + 
                   DHW.MeanMax.YR01_mn_s:Depth_ft_mn_s +
                   Depth_ft_mn_s:SiteSuscp_mn_log_s + ### test interactions b/w suscp and every possibility
                   PAR_surface_mn_s:SiteSuscp_mn_log_s + 
                   kdPAR_mn_s:SiteSuscp_mn_log_s +
                   WaveAction_value_mn_sqrt_s:SiteSuscp_mn_log_s +
                   PctBleached_hist_mn_s:SiteSuscp_mn_log_s +
                   DHW.MeanMax.YR10YR01_mn_s:SiteSuscp_mn_log_s +
                   SST_Variability_AllB4_mn_s:SiteSuscp_mn_log_s +
                   TourRec_10yrAvgPUD_mn_log_s:SiteSuscp_mn_log_s +
                   Urban_runoff_mn_sqrt_s:SiteSuscp_mn_log_s +
                   AgGolf_runoff_mn_sqrt_s:SiteSuscp_mn_log_s +
                   TotalEffluent_mn_sqrt_s:SiteSuscp_mn_log_s,
                 hcbc_complete, 
                 weights = Weights_DriversAndSpatialAnalysis) 
step_mod_bic_2 <- stepAIC(mod_full_2, k = log(nrow(hcbc_complete)))
summary(step_mod_bic_2)
final_mod <- step_mod_bic_2

# perform visual examination of residuals ####
## are the assumptions of normality, homogeneity, and independence met?
par(mfrow = c(2, 2))
plot(final_mod)


### ### ### ### ### ### ### ### ### ### ### 
# PLOTS ####

## i. PARAMETER ESTIMATES +/- SE ####
sum <- summary(final_mod)
sum.co <- data.frame(sum$coefficients)
sum.co$Variable <- rownames(sum.co)
sum.co <- data.frame(sum.co, row.names = NULL)
sum.co <- sum.co[ order(abs(sum.co$Estimate), decreasing = T),]
var_ord <- sum.co$Variable

sum.co$Variable <- factor(sum.co$Variable, levels = var_ord)
sum.co <- sum.co[order(factor(sum.co$Variable, levels = var_ord)),]
sum.co$Variable_plot <- factor(c("x", 
                                 "Acute Thermal Stress x Historical Thermal Stress",
                                 "Historical Bleaching",
                                 "Susceptibility",
                                 "Acute Thermal Stress",
                                 "Urban Run-off",
                                 "Depth",
                                 "Acute Thermal Stress x Tourism",
                                 "Susceptibility x Urban Run-off",
                                 "Historical Bleaching x Susceptibility",
                                 "Surface PAR",
                                 "Effluent",
                                 "Depth x Susceptibility",
                                 "Historical Thermal Stress",
                                 "Tourism"), 
                               levels = c("x", 
                                          "Acute Thermal Stress x Historical Thermal Stress",
                                          "Historical Bleaching",
                                          "Susceptibility",
                                          "Acute Thermal Stress",
                                          "Urban Run-off",
                                          "Depth",
                                          "Acute Thermal Stress x Tourism",
                                          "Susceptibility x Urban Run-off",
                                          "Historical Bleaching x Susceptibility",
                                          "Surface PAR",
                                          "Effluent",
                                          "Depth x Susceptibility",
                                          "Historical Thermal Stress",
                                          "Tourism"))

sum.co$Sig <- NA
sum.co <- transform(sum.co, 
                    Sig=ifelse(Pr...t..<0.05,"Significant","Non significant"))

sum.co$SigLeg <- NA
for(i in c(1:nrow(sum.co))){
  if(sum.co$Pr...t..[i] > 0.05){
    sum.co$SigLeg[i] <- "NS"
  }
  if(sum.co$Pr...t..[i] <= 0.05 & sum.co$Pr...t..[i] > 0.01){
    sum.co$SigLeg[i] <- "p < 0.05"
  }
  if(sum.co$Pr...t..[i] <= 0.01 & sum.co$Pr...t..[i] > 0.001){
    sum.co$SigLeg[i] <- "p < 0.01"
  }
  if(sum.co$Pr...t..[i] <= 0.001){
    sum.co$SigLeg[i] <- "p < 0.001"
  }
}

sum.co <- sum.co[ which(sum.co$Variable != "(Intercept)"),]
sum.co$Variable_plot <- factor(sum.co$Variable_plot)

var_plot <- ggplot(sum.co, aes(x = Variable_plot, y = Estimate)) + 
  geom_hline(yintercept = 0, color = "lightgray") +
  geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error, color = SigLeg), width=.2,
                position=position_dodge(.9))  +
  geom_point(aes(color = SigLeg), size = 3) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = c(0.85,0.1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "lines"), 
    panel.background = element_rect(colour = "black", fill = "white"),
    text = element_text(size = 18)) +
  xlab("") +
  ylab("\nParameter Estimate") +
  scale_y_continuous(limits = c(-1,2)) +
  scale_x_discrete(limits = rev(levels(sum.co$Variable_plot))) +
  scale_color_discrete(name = "")

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Drivers Analysis/Weighted Results")
png(width = 750, height = 750, filename = "ParameterEstimates.png")
var_plot
dev.off()


## ii. MARGINAL RESIDUALS PER VARIABLE ####

# function to predict bleaching and create a plot based on variable of interest
margResidplot.fun <- function(var, mod, dat){
  m <- mod # set final model
  sum_int_f <- summary(m)$coefficients[ which(row.names(summary(m)$coefficients) == "(Intercept)"),"Estimate"] # intercept of final model
  v = var # set variable of interest
  d = dat
  
  # update model to remove variable of interest & calculate residuals
  var.ints <- grep(v, row.names(summary(m)$coefficients), value = T) # this will give a list elements that the variable appears in and we can use it to update the model
  var.update <- paste(var.ints, collapse = "-") # put the list of elements into one string
  mod.var <- update(m, as.formula(paste0(".~.-",var.update))) # create updated model without the variable of interest & its interaction effects
  d$residuals.var <- residuals(mod.var) # save the residuals from the model
  int_var <- coef(mod.var)[1] # intercept of reduced model without variable of interest
  
  # generate list of predictors we need calculate mean values for
  sum.co <- data.frame(summary(m)$coefficients)
  sum.co$Variable <- rownames(sum.co)
  sum.co <- data.frame(sum.co, row.names = NULL)
  # rem <- sum.co[ which(sum.co$Pr...t.. > 0.05),]$Variable
  sum.co <- sum.co[ which(sum.co$Variable != "(Intercept)"),] # remove intercept
  vars = sum.co$Variable[!grepl(paste0(":", collapse = "|"), sum.co$Variable)] # remove interaction effects
  
  # generate new data frame to hold unique values of variable of interest and mean of other predictors
  new_temp <- setNames(data.frame(matrix(ncol = length(vars), nrow = 0)), vars)
  add <- as.data.frame(unique(d[v])) # get unique values of variable of interest
  new_temp <- merge(new_temp, add, all = TRUE) # add in unique values of variable of interest
  # calculate means of other values of interest
  for(i in c(1:length(vars))){
    if(vars[i] != v){
      col.name <- vars[i]
      new_temp[[col.name]] <- mean(d[[col.name]]) # we want the mean of all other other variables - they stay constant while we vary the variable of interest
    }
  }
  
  # make predictions of bleaching on new data
  new_temp$Predict <- predict(m, newdata = new_temp, interval = "confidence")   # add predicted values of bleaching
  new_temp$Predict.fit <- new_temp$Predict[,1]
  new_temp$Predict.lwr <- new_temp$Predict[,2] # confidence interval upper bound
  new_temp$Predict.upr <- new_temp$Predict[,3] # confidence interval lower bound
  
  print(head(new_temp,10)) # check that newdata generated properly
  
  # set color for line
  pcol <- NA
  pv <- summary(m)$coefficients[,4][[v]]
  
  if(pv > 0.05){pcol <- "#F8766D"}
  if(pv <= 0.05 & pv > 0.01){
    pcol <- "#C77CFF"
  }
  if(pv <= 0.01 & pv > 0.001){
    pcol <- "#00BFC4"
  }
  if(pv <= 0.001){
    pcol <- "#7CAE00"
  }
  
  ggplot() + 
    geom_point(data=d, aes(x = .data[[v]], y =residuals.var, size = Inverse_SE)) + 
    geom_ribbon(data=new_temp, aes(x=.data[[v]], ymin=Predict.lwr-int_var, ymax=Predict.upr-int_var), alpha= 0.2, fill="black") +
    geom_line(data=new_temp, aes(x= .data[[v]], y=Predict.fit-int_var), color=pcol, size = 1) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title = element_text(face = "bold"),
      text = element_text(size = 18),
      panel.grid = element_blank()
    ) +
    scale_size_continuous(guide = F) +
    scale_y_continuous(breaks = c(-7.5,-5,-2.5,0,2.5,5,7.5,10), labels = c(-(7.5^2), -(5^2), -(2.5^2), 0, (2.5^2), (5^2), (7.5^2), (10^2))) 
}

par_mn <- mean(hcbc_complete$PAR_surface_mn)
par_sd <- sd(hcbc_complete$PAR_surface_mn)
dhw1_mn <- mean(hcbc_complete$DHW.MeanMax.YR01_mn)
dhw1_sd <- sd(hcbc_complete$DHW.MeanMax.YR01_mn)
dhw10_mn <- mean(hcbc_complete$DHW.MeanMax.YR10YR01_mn)
dhw10_sd <- sd(hcbc_complete$DHW.MeanMax.YR10YR01_mn)
sus_mn <- mean(log(hcbc_complete$SiteSuscp_mn))
sus_sd <- sd(log(hcbc_complete$SiteSuscp_mn))
dep_mn <- mean(hcbc_complete$Depth_ft_mn)
dep_sd <- sd(hcbc_complete$Depth_ft_mn)
hist_mn <- mean(hcbc_complete$PctBleached_hist_mn)
hist_sd <- sd(hcbc_complete$PctBleached_hist_mn)
eff_mn <- mean(sqrt(hcbc_complete$TotalEffluent_mn))
eff_sd <- sd(sqrt(hcbc_complete$TotalEffluent_mn))
urb_mn <- mean(sqrt(hcbc_complete$Urban_runoff_mn))
urb_sd <- sd(sqrt(hcbc_complete$Urban_runoff_mn))
tou_mn <- mean(log(hcbc_complete$TourRec_10yrAvgPUD_mn))
tou_sd <- sd(log(hcbc_complete$TourRec_10yrAvgPUD_mn))

par_resid_plot <- margResidplot.fun("PAR_surface_mn_s", final_mod, hcbc_complete) + xlab("Surface Light (PAR)") + scale_x_continuous(breaks = c((42.5-par_mn)/par_sd, (45-par_mn)/par_sd, (47.5-par_mn)/par_sd), labels = c(42.5, 45, 47.5))  # unscale x axis
dhw1_resid_plot <-  margResidplot.fun("DHW.MeanMax.YR01_mn_s", final_mod, hcbc_complete) + xlab("Acute Thermal Stress (DHW)") + scale_x_continuous(breaks = c((0-dhw1_mn)/dhw1_sd, (3-dhw1_mn)/dhw1_sd, (6-dhw1_mn)/dhw1_sd, (9-dhw1_mn)/dhw1_sd), labels = c(0,3,6,9)) 
dhw10_resid_plot <-  margResidplot.fun("DHW.MeanMax.YR10YR01_mn_s", final_mod, hcbc_complete) + xlab("Historical Thermal Stress (DHW)") + scale_x_continuous(breaks = c((0-dhw10_mn)/dhw10_sd, (3-dhw10_mn)/dhw10_sd, (6-dhw10_mn)/dhw10_sd, (9-dhw10_mn)/dhw10_sd, (12-dhw10_mn)/dhw10_sd), labels = c(0,3,6,9,12)) 
suscp_resid_plot <- margResidplot.fun("SiteSuscp_mn_log_s", final_mod, hcbc_complete) + xlab("Susceptibiity Score") +  scale_x_continuous(breaks = c((log(2)-sus_mn)/sus_sd, (log(2.5)-sus_mn)/sus_sd, (log(3)-sus_mn)/sus_sd, (log(3.5)-sus_mn)/sus_sd, (log(4)-sus_mn)/sus_sd), labels = seq(2,4,by=.5)) 
depth_resid_plot <- margResidplot.fun("Depth_ft_mn_s", final_mod, hcbc_complete) + xlab("Depth (ft)") + scale_x_continuous(breaks = c((10-dep_mn)/dep_sd, (30-dep_mn)/dep_sd, (50-dep_mn)/dep_sd, (70-dep_mn)/dep_sd), labels = c(10,30,50,70)) 
hist_resid_plot <-  margResidplot.fun("PctBleached_hist_mn_s", final_mod, hcbc_complete) + xlab("Historical Bleaching (%)") +  scale_x_continuous(breaks = c((10-hist_mn)/hist_sd, (20-hist_mn)/hist_sd, (30-hist_mn)/hist_sd, (40-hist_mn)/hist_sd, (50-hist_mn)/hist_sd, (60-hist_mn)/hist_sd), labels = c(10, 20,30,40,50,60))
eff_resid_plot <- margResidplot.fun("TotalEffluent_mn_sqrt_s", final_mod, hcbc_complete) + xlab("Sewage Effluent") + scale_x_continuous(breaks = c((sqrt(0)-eff_mn)/eff_sd, (sqrt(1000)-eff_mn)/eff_sd, (sqrt(10000)-eff_mn)/eff_sd, (sqrt(30000)-eff_mn)/eff_sd), labels = c(0,1000,10000,30000)) 
urb_resid_plot <-  margResidplot.fun("Urban_runoff_mn_sqrt_s", final_mod, hcbc_complete) + xlab("Urban Run-off") +  scale_x_continuous(breaks = c((sqrt(0)-urb_mn)/urb_sd, (sqrt(0.05)-urb_mn)/urb_sd, (sqrt(0.2)-urb_mn)/urb_sd, (sqrt(0.4)-urb_mn)/urb_sd, (sqrt(0.6)-urb_mn)/urb_sd), labels = c(0,0.05,0.2,0.4,0.6))
tour_resid_plot <-  margResidplot.fun("TourRec_10yrAvgPUD_mn_log_s", final_mod, hcbc_complete) + xlab("Tourism & Recreation") +  scale_x_continuous(breaks = c((log(1)-tou_mn)/tou_sd, (log(10)-tou_mn)/tou_sd, (log(100)-tou_mn)/tou_sd, (log(1000)-tou_mn)/tou_sd), labels = c(1,10,100,1000))

# save legend as separate plot (use the legend from parameter estimate plots!)
mylegend<-g_legend(ggplot(sum.co, aes(x = Variable_plot, y = Estimate, color = SigLeg)) + theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18), legend.position = "bottom", legend.direction = "horizontal") +
                     geom_line() + scale_color_discrete(name = "Significance"))

# save full plot
setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Drivers Analysis/Weighted Results")
ytitle <- text_grob("Partial Residual", size = 18, face = "bold", rot = 90)
png(width = 1050, height = 950, filename = "Drivers_PartialResiduals.png")
grid.arrange(arrangeGrob(dhw1_resid_plot + ggtitle("a)"), 
                         par_resid_plot + ggtitle("b)"),
                         depth_resid_plot + ggtitle("c)"),
                         suscp_resid_plot + ggtitle("d)"), 
                         dhw10_resid_plot + ggtitle("e)"),
                         hist_resid_plot + ggtitle("f)"), 
                         eff_resid_plot + ggtitle("g)"), 
                         urb_resid_plot + ggtitle("h)"),
                         tour_resid_plot + ggtitle("i)"),
                         nrow = 3), 
             mylegend, nrow = 2, heights = c(10,1),
             left = ytitle)
dev.off()


## iii. INTERACTION SURFACES ####
intSurfplot.fun <- function(var1, var2, mod, dat){
  m <- mod # set final model
  v1 = var1 # set interaction variable 1
  v2 = var2 # set interaction variable 2
  d = dat # set database of observations
  
  v1_seq <- seq(from = min(d[[v1]]), to = max(d[[v1]]), length = 100)
  v2_seq <- seq(from = min(d[[v2]]), to = max(d[[v2]]), length = 100)
  
  # generate list of predictors we need calculate mean values for (everything but the two variables specified)
  sum.co <- data.frame(summary(m)$coefficients)
  sum.co$Variable <- rownames(sum.co)
  sum.co <- data.frame(sum.co, row.names = NULL)
  sum.co <- sum.co[ which(sum.co$Variable != "(Intercept)"),] # remove intercept
  vars = sum.co$Variable[!grepl(paste0(":", collapse = "|"), sum.co$Variable)] # remove interaction effects
  vars2 = vars[!(vars %in% c(v1,v2))] # remove v1 and v2
  
  # generate new data frame to hold unique values of variable of interest and mean of other predictors
  new_temp <- setNames(data.frame(matrix(ncol = length(vars), nrow = 100*100)), vars)
  new_temp[[v1]] <- rep(v1_seq, each = 100)
  new_temp[[v2]] <- rep(v2_seq, 100)
  # calculate means of other values of interest
  for(i in c(1:length(vars2))){
    col.name <- vars2[i]
    new_temp[[col.name]] <- mean(d[[col.name]]) # we want the mean of all other other variables - they stay constant while we vary the variable of interest
  }
  
  # make predictions of bleaching on new data
  new_temp$pred_ble <- predict(m, newdata = new_temp)   # add predicted values of bleaching
  
  new_temp <- transform(new_temp,
                        pred_ble=ifelse(pred_ble>10,10,pred_ble))
  
  new_temp <- transform(new_temp,
                        pred_ble=ifelse(pred_ble<0,0,pred_ble))
  
  ggplot(new_temp, aes(x = .data[[v1]], y = .data[[v2]])) + 
    geom_raster(aes(fill=pred_ble)) + 
    scale_fill_viridis_c(name = "Predicted Bleaching (%)", breaks = c(0,2.5,5,7.5,10), labels = c(0, (2.5^2), 25, (7.5^2), 100), limits=c(-0.1,10.1), na.value = "black") +
    geom_point(data = d, aes(x = .data[[v1]], y = .data[[v2]]),
               shape = 21, size = 3, stroke = 1, color = "black") +
    theme(
      text = element_text(size = 20),
      axis.text = element_text(color = "black")
    )

}

## save the plots per interactions specified ## one for int w/ DHW, one for int w/ suscp
### interactions with acute stress:
dhw_int_plot <- intSurfplot.fun("DHW.MeanMax.YR01_mn_s", "DHW.MeanMax.YR10YR01_mn_s", final_mod, hcbc_complete) +
  ggtitle("a)") +
  xlab("Acute Thermal Stress (DHW)") + 
  scale_x_continuous(expand = c(0,0), breaks = c((0-dhw1_mn)/dhw1_sd, (3-dhw1_mn)/dhw1_sd, (6-dhw1_mn)/dhw1_sd, (9-dhw1_mn)/dhw1_sd), labels = c(0,3,6,9)) +
  ylab("Historical Thermal Stress (DHW)\n") +
  scale_y_continuous(expand = c(0,0), breaks = c((0-dhw10_mn)/dhw10_sd, (3-dhw10_mn)/dhw10_sd, (6-dhw10_mn)/dhw10_sd, (9-dhw10_mn)/dhw10_sd, (12-dhw10_mn)/dhw10_sd), labels = c(0,3,6,9,12)) 

dhw_tour_int_plot <- intSurfplot.fun("DHW.MeanMax.YR01_mn_s", "TourRec_10yrAvgPUD_mn_log_s", final_mod, hcbc_complete) +
  ggtitle("b)") +
  xlab("\nAcute Thermal Stress (DHW)") + 
  scale_x_continuous(expand = c(0,0), breaks = c((0-dhw1_mn)/dhw1_sd, (3-dhw1_mn)/dhw1_sd, (6-dhw1_mn)/dhw1_sd, (9-dhw1_mn)/dhw1_sd), labels = c(0,3,6,9)) +
  ylab("Tourism & Recreation") +
  scale_y_continuous(expand = c(0,0), breaks = c((log(1)-tou_mn)/tou_sd, (log(10)-tou_mn)/tou_sd, (log(100)-tou_mn)/tou_sd, (log(1000)-tou_mn)/tou_sd), labels = c(1,10,100,1000))

# combine plots and save --
dhw_legend <- g_legend(dhw_int_plot +  
                         theme(legend.direction = "horizontal",
                               legend.spacing.y = unit(0, "mm"), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
                               legend.background = element_blank()) +
                         guides(size = guide_legend(nrow = 1)) + 
                         guides(fill=guide_colourbar(barwidth=15))) ## save legend as separate object

dhw <- grid.arrange(dhw_int_plot + xlab("") + theme(legend.position = "none"), 
                    dhw_tour_int_plot + xlab("") + theme(legend.position = "none"), nrow = 2,
                    bottom = text_grob("Acute Thermal Stress (DHW)", size = 20))

### interactions with susceptibility  
sus_dep_int_plot <- intSurfplot.fun("SiteSuscp_mn_log_s", "Depth_ft_mn_s", final_mod, hcbc_complete) +
  ggtitle("c)") +
  xlab("Susceptibility Score") + 
  scale_x_continuous(expand = c(0,0), breaks = c((log(2)-sus_mn)/sus_sd, (log(2.5)-sus_mn)/sus_sd, (log(3)-sus_mn)/sus_sd, (log(3.5)-sus_mn)/sus_sd, (log(4)-sus_mn)/sus_sd), labels = seq(2,4,by=.5)) +  ylab("Depth (ft)") +
  ylab("\nDepth (ft)\n") +
  scale_y_continuous(expand = c(0,0), breaks = c((10-dep_mn)/dep_sd, (30-dep_mn)/dep_sd, (50-dep_mn)/dep_sd, (70-dep_mn)/dep_sd), labels = c(10,30,50,70)) 

sus_hist_int_plot <- intSurfplot.fun("SiteSuscp_mn_log_s", "PctBleached_hist_mn_s", final_mod, hcbc_complete) +
  ggtitle("d)") +
  xlab("Susceptibility Score") + 
  scale_x_continuous(expand = c(0,0), breaks = c((log(2)-sus_mn)/sus_sd, (log(2.5)-sus_mn)/sus_sd, (log(3)-sus_mn)/sus_sd, (log(3.5)-sus_mn)/sus_sd, (log(4)-sus_mn)/sus_sd), labels = seq(2,4,by=.5)) +  
  ylab("\nHistorical Bleaching (%)\n") +
  scale_y_continuous(expand = c(0,0), breaks = c((10-hist_mn)/hist_sd, (20-hist_mn)/hist_sd, (30-hist_mn)/hist_sd, (40-hist_mn)/hist_sd, (50-hist_mn)/hist_sd, (60-hist_mn)/hist_sd), labels = c(10, 20,30,40,50,60))

sus_urb_int_plot <- intSurfplot.fun("SiteSuscp_mn_log_s", "Urban_runoff_mn_sqrt_s", final_mod, hcbc_complete) +
  ggtitle("e)") +
  xlab("Susceptibility Score") + 
  scale_x_continuous(expand = c(0,0), breaks = c((log(2)-sus_mn)/sus_sd, (log(2.5)-sus_mn)/sus_sd, (log(3)-sus_mn)/sus_sd, (log(3.5)-sus_mn)/sus_sd, (log(4)-sus_mn)/sus_sd), labels = seq(2,4,by=.5)) +  
  ylab("\nUrban Run-off") +
  scale_y_continuous(expand = c(0,0), breaks = c((sqrt(0)-urb_mn)/urb_sd, (sqrt(0.05)-urb_mn)/urb_sd, (sqrt(0.2)-urb_mn)/urb_sd, (sqrt(0.4)-urb_mn)/urb_sd, (sqrt(0.6)-urb_mn)/urb_sd), labels = c(0,0.05,0.2,0.4,0.6))

# combine plots and save --
sus_legend <- g_legend(sus_dep_int_plot + 
                         theme(legend.spacing.y = unit(0, "mm"), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
                               legend.background = element_blank()) + 
                         guides(fill=guide_colourbar(barheight=10)) + 
                         guides(size=guide_legend(keywidth=0.35, keyheight=0.35, default.unit="inch")))
btitle <- text_grob("Taxonomic Susceptibility", size = 20)

sus <- grid.arrange(sus_dep_int_plot + xlab("") + theme(legend.position = "none"), 
                    sus_hist_int_plot + xlab("") + theme(legend.position = "none"), 
                    sus_urb_int_plot + xlab("") + theme(legend.position = "none"), 
                    nrow = 2, 
                    bottom = btitle)

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Drivers Analysis/Weighted Results")
png(width = 1450, height = 950, filename = "Drivers_Int_Wtd_all.png")
grid.arrange(arrangeGrob(dhw, sus, nrow = 1, widths = c(2,4)), dhw_legend, heights = c(10,1), nrow = 2)
dev.off()


## iv. MODEL PERTURBATIONS #### 
# take model preds, compare preds of model with data as is, to preds of model with one of these variables tuned down/up
pertDat_fun <- function(var, var_s, trans = "none", dat, mod, pert = 1, summary = "all"){
  v = var
  v_s = var_s
  d = dat
  t = trans
  m = mod
  p = pert
  s = summary
  
  d$predict_ble <- predict(m, d)
  
  if(t == "none"){
    mn_v = mean(d[[v]])
    sd_v = sd(d[[v]])
  }
  
  if(t == "log"){
    mn_v = mean(log(d[[v]]))
    sd_v = sd(log(d[[v]]))
  }
  
  if(t == "sqrt"){
    mn_v = mean(sqrt(d[[v]]))
    sd_v = sd(sqrt(d[[v]]))
  }
  
  sd_pert <- sd(d[[v]])*p
  
  ## add x% of mean
  te_pert <- d[[v]] + sd_pert
  # sometimes adding this will make the value exceed the observed range, which we do not want
  te_pert[ te_pert > max(d[[v]], na.rm=T)] = max(d[[v]], na.rm=T) # anything greater than observed range gets set to obs range
  # scale the new values accordingly
  if(t == "none"){
    d[[v_s]] <- (te_pert- mn_v)/sd_v
  }
  if(t == "log"){
    d[[v_s]] <- (log(te_pert)- mn_v)/sd_v
  }
  if(t == "sqrt"){
    d[[v_s]] <- (sqrt(te_pert)- mn_v)/sd_v
  }
  # make predictions based on new values
  d$pred.pert.plus <- predict(m, d)
  
  ## subtract x% of mean
  te_pert <- d[[v]] - sd_pert
  # or sometimes subtracting this will make value fall below observed range, also bad
  te_pert[ te_pert < min(d[[v]], na.rm=T)] = min(d[[v]], na.rm=T) # anything less than observed range gets set to obs range
  if(t == "none"){
    d[[v_s]] <- (te_pert- mn_v)/sd_v
  }
  if(t == "log"){
    d[[v_s]] <- (log(te_pert)- mn_v)/sd_v
  }
  if(t == "sqrt"){
    d[[v_s]] <- (sqrt(te_pert)- mn_v)/sd_v
  }
  d$pred.pert.minus <- predict(m, d)
  
  d2 <- d[,c("Cluster_DepthBin_ID", "predict_ble", "pred.pert.plus", "pred.pert.minus")]
  
  d2$ble_inc <- d2$pred.pert.plus - d2$predict_ble
  d2$ble_dec <- d2$predict_ble - d2$pred.pert.minus
  
  colnames(d2)[c(3:6)] <- paste(colnames(d2)[c(3:6)],v)
  
  if(s == "all"){
    return(d2)
  }
  
  se_var <- apply(d2[,c(3:4)],2,std.error)
  mn_var <- apply(d2[,c(3:4)],2,mean)
  
  d3 <- data.frame("var" = v, "pert" = c("plus", "minus"), "pred_mn" = mn_var, "pred_se" = se_var, row.names = NULL)
  
  if(s == "mean"){
    return(d3)
  }
}

# calculate mean predicted bleaching by increasing and decreasing each variable by 1SD
#### average change in bleaching by variable ####
avg_pert <- rbind(pertDat_fun("DHW.MeanMax.YR01_mn", "DHW.MeanMax.YR01_mn_s", "none", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("DHW.MeanMax.YR10YR01_mn", "DHW.MeanMax.YR10YR01_mn_s", "none", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("Depth_ft_mn", "Depth_ft_mn_s", "none", hcbc_complete, final_mod,  1, summary = "mean"),
                  pertDat_fun("PctBleached_hist_mn", "PctBleached_hist_mn_s", "none", hcbc_complete, final_mod, 1, summary = "mean"), 
                  pertDat_fun("SiteSuscp_mn", "SiteSuscp_mn_log_s", "log", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("PAR_surface_mn", "PAR_surface_mn_s", "none", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("TotalEffluent_mn", "TotalEffluent_mn_sqrt_s", "sqrt", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("TourRec_10yrAvgPUD_mn", "TourRec_10yrAvgPUD_mn_log_s", "log", hcbc_complete, final_mod, 1, summary = "mean"),
                  pertDat_fun("Urban_runoff_mn", "Urban_runoff_mn_sqrt_s", "sqrt", hcbc_complete, final_mod, 1, summary = "mean"))          
base_pred <- mean(predict(final_mod, hcbc_complete))

cols <- c("inc" = "orange", "dec" = "cadetblue3")

avg_pert_plot <- ggplot() +
  geom_hline(yintercept = base_pred) +
  geom_point(data = avg_pert[ which(avg_pert$pert == "plus"),], aes(x = var, y = pred_mn, color = "inc")) +
  geom_errorbar(data = avg_pert[ which(avg_pert$pert == "plus"),], aes(x = var, ymin = pred_mn - pred_se, ymax = pred_mn + pred_se, color = "inc"), width=.2,
                position=position_dodge(.9)) +
  geom_point(data = avg_pert[ which(avg_pert$pert == "minus"),], aes(x = var, y = pred_mn, color = "dec")) +
  geom_errorbar(data = avg_pert[ which(avg_pert$pert == "minus"),], aes(x = var, ymin = pred_mn - pred_se, ymax = pred_mn + pred_se, color = "dec"), width=.2,
                position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 48, hjust = 1),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        text = element_text(size = 16)) +
  ylab("Predicted Bleaching (%)\n") +
  xlab("") +
  scale_x_discrete(labels = c("Acute Thermal Stress (DHW)", "Historic Thermal Stress (DHW)", "Depth (ft)", "Historic Bleaching (%)", 
                              "Susceptibility", "Surface Light (PAR)", "Sewage Effluent", "Tourism & Recreation", "Urban Run-off")) +
  scale_y_continuous(breaks = c(3,4,5), labels = c(9,16,25)) +
  scale_color_manual(name = "", 
                     breaks = c("inc", "dec"), 
                     values = cols,
                     labels = c("Increase (+1 SD)", "Decrease (-1 SD)"))


setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Drivers Analysis/Weighted Results")
png(width = 855, height = 650, filename = "Avg_Ble_Pert_1SD.png")
avg_pert_plot
dev.off()

### ### ### ### ### ### ### ### ###
### v. PLOT SITE LEVEL PERTURBS ###
### ### ### ### ### ### ### ### ###

## only look @ variables that can be effectively managed: run-off, effluent, PAR, suscp, & tourism
# we will hold acute thermal stress @ 90-95th percentile instead of at mean to simulate a heating event
mn_dhw1 <- mean(hcbc_complete$DHW.MeanMax.YR01_mn) # calc mn and sd for scaling the 95th percentile
sd_dhw1 <- sd(hcbc_complete$DHW.MeanMax.YR01_mn)
# calculate 95th percentile of acute thermal stress
dhw1_95 <- quantile(hcbc_complete$DHW.MeanMax.YR01_mn, 0.95)[[1]]
hcbc_complete2 <- hcbc_complete
hcbc_complete2$DHW_YR01_mn_s <- (dhw1_95 - mn_dhw1)/sd_dhw1
sus.th <- pertDat_fun("SiteSuscp_mn", "SiteSuscp_mn_log_s", "log", hcbc_complete2, final_mod, 1)
par.th <- pertDat_fun("PAR_surface_mn", "PAR_surface_mn_s", "none", hcbc_complete2, final_mod, 1)
eff.th <- pertDat_fun("TotalEffluent_mn", "TotalEffluent_mn_sqrt_s", "sqrt", hcbc_complete2, final_mod, 1)
urb.th <- pertDat_fun("Urban_runoff_mn", "Urban_runoff_mn_sqrt_s", "sqrt", hcbc_complete2, final_mod, 1)
tour.th <-  pertDat_fun("TourRec_10yrAvgPUD_mn", "TourRec_mn_log_s", "log", hcbc_complete2, final_mod, 1)
all.th <- Reduce(function(...) merge(..., all.x=TRUE), list(tour.th, sus.th, par.th, eff.th, urb.th))
pred_dec.th <- all.th[,c(1,6,10,14,18,22)]
meta <- hcbc_complete[,c("ZoneName", "Island_Name", "Bluster", "Depth_bin", "Longitude_DD_mn", "Latitude_DD_mn")]
pred_dec.th <- merge(pred_dec.th, meta, all.x = T)
pred_dec.th$FirstLargestDec <- colnames(pred_dec.th[,c(2:6),])[max.col(pred_dec.th[,c(2:6),],ties.method="first")] # first highest
pred_dec.th$SecondLargestDec <- colnames(pred_dec.th[,c(2:6),])[apply(pred_dec.th[,c(2:6),],1,function(x)which(x==sort(x,partial=4)[4])[1])] # second highest
pred_dec.th$ThirdLargestDec <- colnames(pred_dec.th[,c(2:6),])[apply(pred_dec.th[,c(2:6),],1,function(x)which(x==sort(x,partial=3)[3])[1])] # third highest
pred_dec.th <- pred_dec.th %>%
  mutate(FirstLargestDec = str_remove_all(FirstLargestDec, "ble_dec "))
pred_dec.th <- pred_dec.th %>%
  mutate(SecondLargestDec = str_remove_all(SecondLargestDec, "ble_dec "))
pred_dec.th <- pred_dec.th %>%
  mutate(ThirdLargestDec = str_remove_all(ThirdLargestDec, "ble_dec "))
pred_pert_dat <- pred_dec.th[,c(1,7:14)]

# function to create maps based upon data and model and specific ranking 
map_p_fun <- function(dat, fill.var, mod){
  var = fill.var
  d = dat
  m = mod
  
  # import zone shapefile - readOGR() converts the shapefile into a Spatial Polygons Data Frame
  zones = readOGR("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/Current Database/For Analysis/Shapefiles/HCBC_Zones.shp")
  
  sum.co <- data.frame(summary(m)$coefficients)
  sum.co$Variable <- rownames(sum.co)
  vars <- sum.co$Variable[c(2:10)]
  vars <- gsub("(mn).*", "\\1", vars)
  
  # set the desired color for each variable that we are altering
  col = setNames(hue_pal()(5), levels(as.factor(vars[c(5:9)])))
  
  # convert the 2019 bleaching observations to a spatial points data frame for mapping
  h.sp = SpatialPointsDataFrame(coords=d[,c("Longitude_DD_mn","Latitude_DD_mn")], data=d)
  crs(h.sp) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  scale = 1.3 # set the scale by which to expand the map extent (generaly 1.3 works for all islands with exceptions)
  
  ## ## ## ## ## ##  
  ### MAUI NUI ####
  ## ## ## ## ## ##
  isl_bb.maui = scale*extent(h.sp[ which(h.sp$Island_Name %in% c("Maui","Lanai")),]) # set the extent (xmin, xmax, ymin, ymax) for the Maui base map
  
  # convert zone shapefile to simple feature object (sf) for stable plotting in ggplot
  ISL_zone.maui <- st_as_sf(zones)
  ISL_zone.maui <- ISL_zone.maui[ which(ISL_zone.maui$ISLAND %in% c("Maui","Lanai")),] # only include Maui zones for this map
  
  isl_base.maui = get_map(location = c(mean(isl_bb.maui[1:2]), mean(isl_bb.maui[3:4])), zoom=10, maptype="satellite") # zoom 10 works for nui
  
  # create a dataframe that has coordinates for the labels of each of Maui's zones
  zone_points.maui <- sf::st_point_on_surface(ISL_zone.maui)
  zone_coords.maui <- as.data.frame(sf::st_coordinates(zone_points.maui))
  zone_coords.maui$NAME <- ISL_zone.maui$ZoneName
  zone_coords.maui$NAME <- c("South", "West", "West\nNorthwest", "North", "Northwest", "South", "Northeast", "Southeast") # i know from looking at zone_points that this is the order of zone names
  zone_coords.maui$X <- c(-156.4, -156.55, -156.735, -156.42, -156.63, -156.9122, -156.8186, -156) # manually set the coordinates for where the zone names are located on the map
  zone_coords.maui$Y <- c(20.5, 20.69615, 20.79, 21, 21.11, 20.675, 20.95, 20.58)
  zone_coords.maui <- zone_coords.maui[ which(zone_coords.maui$NAME %!in% c("North","Southeast")),]
  ISL_zone.maui <- ISL_zone.maui[ which(ISL_zone.maui$ZoneName %!in% c("Maui_N","Maui_SE")),] 
  
  maui_p= ### plot for all variables perturbed
    ggmap(isl_base.maui) + # plots base map
    geom_sf(data = ISL_zone.maui,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones, transparent with white border
    geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = .data[[var]], shape = Depth_bin), color = "black", data = as.data.frame(d), size = 3) +
    scale_color_manual(guide = F) +
    scale_fill_manual(values = col, drop = F) + # sets the palette for coloring the cluster points
    scale_shape_manual(values = c(24,25)) + 
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    geom_text(data = zone_coords.maui, aes(X, Y, label = NAME), colour = "white", size = 6) + # adds the text 
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      title = element_text(size = 18),
      legend.position = "none",
      legend.text = element_text(size = 16),
      legend.box="vertical", 
      legend.margin=margin()
    ) +
    coord_sf(expand = F) + # if you don't include this then gray space remains after you crop the map
    # crop the map to remove extra space and manually set the breaks/labels
    scale_x_continuous(limits = c(-157.05, -156.25), breaks = c(-157, -156.8, -156.6, -156.4), labels = c(-157, -156.8, -156.6, -156.4)) + 
    scale_y_continuous(limits = c(20.47,21.13), breaks = c(20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1), labels = c(20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1)) +
    ggtitle("Maui Nui") 
  
  ## ## ## ## ##
  ## HAWAII ####
  ## ## ## ## ##
  ## since Hawaii is so big and only want to show data from 2 zones, create inset maps for each zone
  isl_bb.haw = scale*extent(h.sp[ which(h.sp$Island_Name == "Hawaii"),])
  ISL_zone.haw <- st_as_sf(zones)
  ISL_zone.haw <- st_as_sf(crop(x=gBuffer(zones,byid = T,width=0),y=isl_bb.haw))
  isl_base.haw = get_map(location = c(-156, mean(isl_bb.haw[3:4])), zoom=7, maptype="satellite")
  
  hawaii=ggmap(isl_base.haw) +
    geom_sf(data = ISL_zone.haw,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
    coord_sf(expand = F) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      legend.position = "none"
    ) +
    scale_x_discrete(breaks = c(-159, -158, -157, -156, -155, -154, -153), labels = c(-159, -158, -157, -156, -155, -154, -153)) +
    scale_y_discrete(breaks = c(17, 18, 19, 20, 21, 22), labels = c(17, 18, 19, 20, 21, 22))
  
  ## HAWAII NW ###
  # create inset map for NW zone
  isl_bb.haw.nw = scale*extent( h.sp[ which(h.sp$ZoneName == "Hawaii_NW"),])
  ISL_zone.haw.nw <- st_as_sf(zones[ which(zones$ZoneName == "Hawaii_NW"),])
  isl_base.haw.nw = get_map(location = c(mean(isl_bb.haw.nw[1:2]), mean(isl_bb.haw.nw[3:4])), zoom=9, maptype="satellite")
  
  hawaii_nw=ggmap(isl_base.haw.nw) +
    geom_sf(data = ISL_zone.haw.nw,fill=NA,color="white", inherit.aes = FALSE) +
    geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = .data[[var]]), data = as.data.frame(pred_pert_dat[ which(pred_pert_dat$ZoneName == "Hawaii_NW"),]), size = 4) +
    scale_fill_manual(values = col, drop = F) +
    scale_shape_manual(values = c(24,25)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = c(-156.1,-155.82), breaks = seq(-156.1, -155.8, len = 4)) +
    scale_y_continuous(limits = c(19.46, 20.18))
  
  ## HAWAII SW ###
  # create inset map for SW zone
  isl_bb.haw.sw = scale*extent(h.sp[ which(h.sp$ZoneName == "Hawaii_SW"),])
  ISL_zone.haw.sw <- st_as_sf(zones[ which(zones$ZoneName == "Hawaii_SW"),])
  isl_base.haw.sw = get_map(location = c(mean(isl_bb.haw.sw[1:2]), mean(isl_bb.haw.sw[3:4])), zoom=9, maptype="satellite")
  
  hawaii_sw=ggmap(isl_base.haw.sw) +
    geom_sf(data = ISL_zone.haw.sw,fill=NA,color="white", inherit.aes = FALSE) +
    geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = .data[[var]]), data = as.data.frame(pred_pert_dat[ which(pred_pert_dat$ZoneName == "Hawaii_SW"),]), size = 4) +
    scale_fill_manual(values = col, drop = F) +
    scale_shape_manual(values = c(24,25)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = c(-156.035,-155.85)) +
    scale_y_continuous(limits = c(19, 19.51)) 
  
  ## COMBINED HAWAII ###
  # now start to combine the overall map with the inset maps
  remove <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(),
                  panel.border = element_rect(colour = "white", fill=NA, size=2))
  h_nw <- hawaii_nw + remove + ggtitle("  Northwest") + theme(plot.title = element_text(size = 16, colour = "white", margin = margin(t = 4, b = -20)))
  h_sw <- hawaii_sw + remove + ggtitle("  Southwest") + theme(plot.title = element_text(size = 16, colour = "white", margin = margin(t = 4, b = -20)))
  nw <- ggplotGrob(h_nw)
  sw <- ggplotGrob(h_sw)
  
  # create final map
  hawaii_all_p <- hawaii + 
    geom_segment(aes(x = -156.14, y = 19.42, xend = -157, yend = 16.75), colour="white") +
    annotation_custom(grob = nw, xmin = -159.8, xmax = -156.5,
                      ymin = 16.5, ymax = 22.7) +
    geom_rect(aes(xmin = -156.14, xmax = -155.79, ymin = 19.42, ymax =  20.25), # nw
              fill = "transparent", color = "white", size = 1) +
    
    geom_segment(aes(x = -155.83, y = 18.96, xend = -154.3, yend = 16.9), colour="white") +
    annotation_custom(grob = sw, xmin = -154.6, xmax = -152.6,
                      ymin = 16, ymax = 22) +
    geom_rect(aes(xmin = -156.035, xmax = -155.83, ymin = 18.96, ymax =  19.54), # sw
              fill = "transparent", color = "white", size = 1) +
    
    ggtitle("Hawaii") + theme(title = element_text(size = 18))
  
  
  ## ## ## ## ##
  #### OAHU ####
  ## ## ## ## ##
  isl_bb.oahu <- scale*extent(h.sp[ which(h.sp$Island_Name == "Oahu"),])
  ISL_zone.oahu <- st_as_sf(zones)
  ISL_zone.oahu <- ISL_zone.oahu[ which(ISL_zone.oahu$ISLAND == "Oahu"),]
  
  isl_base.oahu = get_map(location = c(mean(isl_bb.oahu[1:2]), mean(isl_bb.oahu[3:4])), zoom=10, maptype="satellite")
  zone_points.oahu <- sf::st_point_on_surface(ISL_zone.oahu)
  zone_coords.oahu <- as.data.frame(sf::st_coordinates(zone_points.oahu))
  zone_coords.oahu$NAME <- ISL_zone.oahu$ZoneName
  zone_coords.oahu$NAME <- c("South","West","East","North")
  zone_coords.oahu$X <- c(-157.9699, -158.290, -157.63, -157.9038)
  zone_coords.oahu$Y <- c(21.227, 21.41524, 21.44492, 21.77)
  zone_coords.oahu <- zone_coords.oahu[ which(zone_coords.oahu$NAME %in% c("South", "East")),]
  
  oahu_p=ggmap(isl_base.oahu) +
    geom_sf(data = ISL_zone.oahu[ which(ISL_zone.oahu$ZoneName %in% c("Oahu_S", "Oahu_E")),],fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
    geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = .data[[var]]), data = as.data.frame(pred_pert_dat), color = "black", size = 4) +
    scale_fill_manual(values = col, drop = F) +
    scale_shape_manual(values = c(24,25)) +
    coord_sf(expand = F) +
    geom_text(data = zone_coords.oahu, aes(X, Y, label = NAME), colour = "white", size = 8) + # labels, remove for NW and SW Hawaii
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      legend.position = "none",
      title = element_text(size = 18),
      legend.key.width = unit(1, "inch") 
    ) +
    ggtitle("Oahu") +
    scale_x_continuous(limits = c(-158.2, -157.59), breaks = c(-158.1, -157.9, -157.7), labels = c(-158.1, -157.9, -157.7)) +
    scale_y_continuous(limits = c(21.18, 21.72), breaks = c(21.2, 21.3, 21.4, 21.5, 21.6, 21.7), labels = c(21.2, 21.3, 21.4, 21.5, 21.6, 21.7))
  
  ### create dummy plot for extracting a legend ###
  if(var == "MaxDec_All"){
    leg.title <- "Greatest Simulated Reduction in Bleaching \nAttributed to Decrease in:"
  }
  if(var == "FirstLargestDec"){
    leg.title <- "Greatest Simulated Reduction in Bleaching \nAttributed to Decrease in:"
  }
  if(var == "SecondLargestDec"){
    leg.title <- "Second Greatest Simulated Reduction in Bleaching \nAttributed to Decrease in:"
  }
  if(var == "ThirdLargestDec"){
    leg.title <- "Third Greatest Simulated Reduction in Bleaching \nAttributed to Decrease in:"
  }
  
  vars <- vars[c(5:9)]
  pred_pert_dat$leg.var <- rep(vars, length = nrow(pred_pert_dat))
  key <-
    ggmap(isl_base.oahu) +
    geom_sf(data = ISL_zone.oahu,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
    geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = leg.var), data = as.data.frame(pred_pert_dat), color = "black", size = 4) +
    scale_fill_manual(name = leg.title,
                      values = col, 
                      labels = c("Surface Light (PAR)",
                                 "Taxonomic Susceptibility Score",
                                 "Sewage Effluent",
                                 "Tourism & Recreation",
                                 "Urban Run-off")) +
    scale_shape_manual(name = "", values = c(24,25)) +
    guides(fill=guide_legend(nrow = 6, title.position="top", title.hjust = 0.5, override.aes=list(shape=21))) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 18, colour = "black"),
      title = element_text(size = 18),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.box="vertical", 
      legend.margin=margin()
    )
  
  pert_leg <- g_legend(key)
  
  ### put all three island maps together ###
  grid.arrange(maui_p + theme(plot.margin=unit(c(0,0,0,0),"cm")),
               hawaii_all_p + theme(plot.margin=unit(c(0,0,0,0),"cm")), 
               oahu_p + theme(plot.margin=unit(c(0,0,0,0),"cm")), 
               pert_leg, nrow = 2)
}

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Drivers Analysis/Weighted Results")
png(width = 900, height = 850, filename = "Drivers_Pert_1SD_1st_Map.png")
map_p_fun(pred_pert_dat, "FirstLargestDec", final_mod)
dev.off()

sim_2 <- map_p_fun(pred_pert_dat, "SecondLargestDec", final_mod)
sim_3 <- map_p_fun(pred_pert_dat, "ThirdLargestDec", final_mod)

png(width = 900, height = 850, filename = "Drivers_Pert_1SD_2nd_Map.png")
dev.off()

png(width = 900, height = 850, filename = "Drivers_Pert_1SD_3rd_Map.png")
dev.off()
