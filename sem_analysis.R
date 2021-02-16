rm(list=ls())
cat("\14")
graphics.off()

library(bigleaf)
library(lavaan)
library(semPlot)
library(car)
library(lavaanPlot)
library(lubridate)

wdir <- "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"

setwd(paste0(wdir,"/data/ancillary_data"))

phenoflux_metadata <- read.csv("pheno_flux_sites_to_use.csv")

sites <- phenoflux_metadata$fluxsite
vegtypes <- phenoflux_metadata$vegtype
vegtypes[sites == "US-Ro4"] <- "GR"
vegtypes[sites == "US-Mpj" | sites == "US-Ton"] <- "SA"

##########################
mdl = 
'
VPD~t_air+precip_15day+gcc
gcc~gdd+cdd+precip_15day+VPD
Gs~gcc+VPD+precip_15day
EF~Gs+VPD
'

transformations = data.frame(var=c("VPD","precip_15day","Gs"),transformation=c("log","sqrt","log"))

resp_vars <- c()
pred_vars <- c()

mdl_parts <- strsplit(mdl,"\n")[[1]]
mdl_parts <- mdl_parts[2:length(mdl_parts)]

for (m in 1:length(mdl_parts)){
  
  mdl_m = mdl_parts[m]
  mdl_m_split = strsplit(mdl_m,"\\~")[[1]]
  
  resp_vars <- c(resp_vars,mdl_m_split[1])
  
  pred_vars_m <- strsplit(mdl_m_split[2],"\\+")[[1]]
  
  for (p in 1:length(pred_vars_m)){
    
    pred_vars_mp <- pred_vars_m[p]
    
    if (any(pred_vars == pred_vars_mp)){
      
      next
      
    }
    
    pred_vars <- c(pred_vars,pred_vars_mp)
    
  }
  
}

col_names <- c("site","response","vegtype",pred_vars)

gof_stats <- matrix(NA,nrow = length(sites),2)
r2 <- matrix(NA,nrow = length(sites),length(resp_vars))

rownames(gof_stats) <- sites; colnames(gof_stats) <- c("CFI","SRMR")
colnames(r2) <- resp_vars

params <- matrix(-9999,nrow=length(sites)*length(resp_vars),length(col_names))
colnames(params) <- col_names
params[,1] <- rbind(t(as.matrix(sites)),t(sites),t(as.matrix(sites)),t(as.matrix(sites)))#,t(as.matrix(sites)))
params[,2] <- rep(resp_vars,length(sites))
params[,3] <- rbind(t(as.matrix(vegtypes)),t(vegtypes),t(as.matrix(vegtypes)),t(as.matrix(vegtypes)))#,t(as.matrix(vegtypes)))
params <- as.data.frame(params)

params_se <- matrix(-9999,nrow=length(sites)*length(resp_vars),length(col_names))
colnames(params_se) <- col_names
params_se[,1] <- rbind(t(as.matrix(sites)),t(sites),t(as.matrix(sites)),t(as.matrix(sites)))#,t(as.matrix(sites)))
params_se[,2] <- rep(resp_vars,length(sites))
params_se[,3] <- rbind(t(as.matrix(vegtypes)),t(vegtypes),t(as.matrix(vegtypes)),t(as.matrix(vegtypes)))#,t(as.matrix(vegtypes)))
params_se <- as.data.frame(params_se)

nyr = matrix(NA,nrow=length(sites),ncol=1)

##########################################################
for (i in 1:length(sites)){
  
  setwd(paste0(wdir,"/results/4_daily_flux_data/"))
  data = read.csv(sprintf("%s_daily_values.csv",sites[i]))
  data[data == -9999] = NA
  
  # Shift baseline gcc values for US-UMB due to change in camera FOV
  if (sites[i] == "US-UMB"){
    
    data$gcc[as.Date(data$date)>"2014-05-01"] <- data$gcc[as.Date(data$date)>"2014-05-01"] + 0.06
    
  }
  
  # Remove 2014 data
  if (sites[i] == "US-KFS"){
    
    data$to_remove_id[year(as.Date(data$date)) == 2014] <- 1
    
  }
  
  data$gdd[data$gdd < 0] <- NA
  data$cdd[data$cdd < 0] <- NA
  
  data$t_surf = radiometric.surface.temp(data = data,
                                         LW_up = "lw_out",
                                         LW_down = "lw_in",
                                         emissivity = phenoflux_metadata$emissivity[i])$Trad_degC
  
  data$Gah <- data$H/(air.density(data$t_air,data$pressure) * bigleaf.constants()$cp * (data$t_surf-data$t_air))
  
  data$to_remove_id[data$t_surf - data$t_air < 0] <- 1
  data$to_remove_id[data$Gs < 0] <- 1
  data$to_remove_id[data$Gah > 0.5 | data$Gah <= 0] <- 1
  
  model_matrix = data[,unique(c(resp_vars,pred_vars))]
  
  for (t in 1:length(transformations)){
    
    col_id <- which(colnames(model_matrix) == transformations$var[t])
    
    var_t <- model_matrix[,col_id]
    var_t[var_t <= 0] <- NA
    var_t <- apply(cbind(var_t),1,FUN=transformations$transformation[t])
    
    model_matrix[,col_id] <- var_t
    
  }
  
  # model_matrix[is.infinite(model_matrix)] <- NA
  complete_cases <- data$to_remove_id == 0 # & !is.na(rowSums(model_matrix))
  model_matrix = model_matrix[complete_cases,]
  
  mm_std = scale(model_matrix)
  sem_mdl_std = lavaan::sem(model = mdl_w_defined_params,data=as.data.frame(mm_std))
  sem_results_std = summary(sem_mdl_std,fit=TRUE,rsquare=TRUE)
  
  semPaths(sem_mdl_std,
           'std',
           layout='tree2',
           fade=TRUE)
  
  mlr_models <- list()
  
  for (k in 1:length(resp_vars)){
    
    lm_i = eval(parse(text = sprintf("lm(%s,data=as.data.frame(mm_std))",mdl_parts[k])))
    mlr_models[[k]] <- lm_i
    avp <- avPlots(lm_i,main=sites[i])
    
    resp_var_k <- resp_vars[k]
    pred_vars_k <- strsplit(mdl_parts[k],paste0(resp_var_k,"~"))[[1]][2]
    pred_vars_k <- strsplit(pred_vars_k,"\\+")[[1]]
    
    mtx_row_id <- which(params$site == sites[i] & params[,2] == resp_var_k)
    
    for (v in 1:length(pred_vars_k)){
      
      id_v <- which(sem_results_std$PE$lhs == resp_var_k & 
                      sem_results_std$PE$rhs == pred_vars_k[v])
      
      mtx_col_id <- which(colnames(params) == pred_vars_k[v])
      
      params[mtx_row_id,mtx_col_id] <- round(sem_results_std$PE$est[id_v],3)
      params_se[mtx_row_id,mtx_col_id] <- round(sem_results_std$PE$se[id_v],3)
      
      pp_plot_data <- data.frame(avp[[v]])
      pp_plot_data[,paste0(resp_var_k,"_pred")] <- as.vector(predict(lm(pp_plot_data[,2] ~ pp_plot_data[,1])))
      
      setwd(paste0(wdir,"/results/sem_output/v3"))
      write.csv(pp_plot_data,file = sprintf("%s_partReg_%s-%s.csv",sites[i],resp_var_k,pred_vars_k[v]),
                row.names = FALSE)       
      
    }
    
    r2_id <- which(sem_results_std$PE$op == "r2" &  
                     sem_results_std$PE$lhs == resp_var_k)
    r2[i,k] <- round(sem_results_std$PE$est[r2_id],2)
    
  }
  
  gof_stats[i,1] <- sem_results_std$FIT[9]
  gof_stats[i,2] <- sem_results_std$FIT[21]
  
  nyr[i] = sum(complete_cases)/365
  
  # }
  
}

gof_stats <- data.frame(gof_stats)
r2 <- data.frame(r2)
params <- data.frame(params)
params_se <- data.frame(params_se)

setwd(paste0(wdir,"/results/sem_output/v3"))
write.csv(params,"sem_param_est.csv",row.names = FALSE)
write.csv(params_se,"sem_param_se.csv",row.names = FALSE)
write.csv(gof_stats,"sem_gof_stats.csv")
write.csv(r2,"rsquared_stats.csv",row.names = FALSE)
write.csv(nyr,"n_sites_years.csv",row.names = FALSE)
