###### EPU, monetary policy, and housing price in the US ######
#Author: Marvin Huang #

# load required packages
library(svars)
library(zoo)
library(quantmod)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(tseries)
library(forecast)

# load the dataset
macro <- read.csv(file.choose())
head(macro)

# convert to times series object
pop <- ts(macro$'population', start = c(1985,1,1), frequency = 4)
epu <- ts(macro$'epu', start = c(1985,1,1), frequency = 4)
rgdp <- ts(macro$'rgdp', start = c(1985,1,1), frequency = 4)
cpi <- ts(macro$'cpi', start = c(1985,1,1), frequency = 4)
ffr <- ts(macro$'ffr', start = c(1985,1,1), frequency = 4)
m2 <- ts(macro$'m2', start = c(1985,1,1), frequency = 4)
mspus <- ts(macro$'mspus', start = c(1985,1,1), frequency = 4)
u <- ts(macro$'u', start = c(1985,1,1), frequency = 4)

# adjust the time variable
macro$time <- as.yearqtr(macro$time, format = "%m/%d/%y")
macro$date <- as.Date(macro$time)



### 1. generate EPU on USREC plot
# Fetch the recession data
getSymbols("USREC", src="FRED")

# Convert to data.frame
epu.df <- data.frame(date = as.Date(macro$date), epu = macro$epu)

# Calculate the start and end dates of the recessions
start <- index(USREC[which(diff(USREC$USREC) == 1)])
end   <- index(USREC[which(diff(USREC$USREC) == -1) - 1])
recession.df <- data.frame(start = start, end = end[-1])
recession.df <- subset(recession.df, start >= min(epu.df$date))

p_epu <- ggplot() +
  # Plotting EPU data with thicker line
  geom_line(data = epu.df, aes(x = date, y = epu), linewidth = 1) +
  # Adding recession shades
  geom_rect(data = recession.df,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
            fill = "grey", alpha = 0.4) +
  # Adding a red line at y = 100
  geom_hline(yintercept = 100, color = "red", linetype = "solid", size = 0.5) +
  # Labels and title
  labs(
    x = "Date",
    y = "EPU"
  ) +
  # Adjusting the x-axis to have less frequent date breaks
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +   # Change here: from "2 years" to "5 years"
  # Adjusting the y-axis to have more frequent breaks
  scale_y_continuous(breaks = seq(floor(min(epu.df$epu, na.rm = TRUE)/50)*50, 
                                  ceiling(max(epu.df$epu, na.rm = TRUE)/50)*50, 
                                  by = 50)) +
  # A more professional theme
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "inch")  # Adjust the plot margins to make the main plot area smaller
  )
show(p_epu)



### 2.generate SVAR and IRF plots
# generate first-difference variables
diff_pop <- diff(pop)
diff_rgdp <- diff(rgdp)
diff_cpi <- diff(cpi)
diff_mspus <- diff(mspus)

# SVAR restrictions
amat <- diag(6)
amat[2,1] <- NA
amat[3,1] <- NA
amat[3,2] <- NA
amat[4,1] <- NA
amat[4,2] <- NA
amat[4,3] <- NA
amat[5,1] <- NA
amat[5,2] <- NA
amat[5,3] <- NA
amat[5,4] <- NA
amat[6,1] <- NA
amat[6,2] <- NA
amat[6,3] <- NA
amat[6,4] <- NA
amat[6,5] <- NA

# built the model
sv <- cbind(epu, diff_rgdp, diff_cpi, u, ffr, diff_mspus)
sv <- window(sv, start=c(1985, 2))

# lag order selection
lagselect <- VARselect(sv, lag.max = 8, type = "const")
lagselect$selection

# Estimate the model
Model1 <- VAR(sv, p = 3, season = NULL, exogen = NULL, type = "const")
SVARMod1 <- SVAR(Model1, Amat = amat, Bmat = NULL, hessian = TRUE, estmethod = c("scoring", "direct"))
SVARMod1

# impulse response function
# I. Generate All IRFs:
variables <- colnames(sv)
irf_list <- list()

for (impulse_var in variables) {
  for (response_var in variables) {
    current_irf <- irf(SVARMod1, impulse = impulse_var, response = response_var, n.ahead = 48, ortho = TRUE, 
                       cumulative = FALSE, boot = TRUE, ci = 0.95, runs = 500)
    irf_list[[paste(impulse_var, response_var, sep = "_to_")]] <- current_irf
  }
}

# II. Plotting All IRFs in a Grid:
# Convert IRFs to data.frames for ggplot
df_list <- lapply(names(irf_list), function(name) {
  df <- data.frame(Time = 1:length(irf_list[[name]]$irf[[1]]))
  
  df$Estimate <- irf_list[[name]]$irf[[1]]
  df$Lower <- irf_list[[name]]$Lower[[1]]
  df$Upper <- irf_list[[name]]$Upper[[1]]
  
  df$Impulse <- strsplit(name, "_to_")[[1]][1]
  df$Response <- strsplit(name, "_to_")[[1]][2]
  return(df)
})

# Bind all data frames in the list into one data frame
combined_df <- dplyr::bind_rows(df_list)

# Set the order for Facet
ordering <- c("epu", "diff_rgdp", "diff_cpi", "u", "ffr", "diff_mspus")
combined_df$Impulse <- factor(combined_df$Impulse, levels = ordering)
combined_df$Response <- factor(combined_df$Response, levels = ordering)

facet_ordering <- as.vector(outer(ordering, ordering, FUN = paste, sep=" ~ "))
combined_df$Facet <- factor(combined_df$Facet, levels = facet_ordering)

# Plot using ggplot
ggplot(combined_df, aes(x = Time, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid", size = 0.3) +  # Adding the red line at y=0
  facet_wrap(~ Facet, scales = "free", ncol = length(ordering)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )


# Robustness check 1
# ST method
reduced.form <- vars::VAR(sv, p = 3, type = "const")
structural.form <- id.st(reduced.form, nc = 1, c_lower = 0.3, c_upper = 0.7, c_step = 5, c_fix = NULL,
                         transition_variable = NULL, gamma_lower = -3, gamma_upper = 2,
                         gamma_step = 0.5, gamma_fix = NULL, max.iter = 5, crit = 0.001,
                         restriction_matrix = NULL, lr_test = FALSE)
summary(structural.form)

# plot IRF
cores <- parallel::detectCores() - 1
boot.svar <- wild.boot(structural.form, n.ahead = 48, nboot = 500, nc = cores)
plot(boot.svar)


# Robustness check 2
# generate first-difference variables
diff_epu <- diff(epu)
diff_rgdp <- diff(rgdp)
diff_cpi <- diff(cpi)
diff_mspus <- diff(mspus)
diff_u <- diff(u)
diff_ffr <- diff(ffr)

# SVAR restrictions
amat <- diag(6)
amat[2,1] <- NA
amat[3,1] <- NA
amat[3,2] <- NA
amat[4,1] <- NA
amat[4,2] <- NA
amat[4,3] <- NA
amat[5,1] <- NA
amat[5,2] <- NA
amat[5,3] <- NA
amat[5,4] <- NA
amat[6,1] <- NA
amat[6,2] <- NA
amat[6,3] <- NA
amat[6,4] <- NA
amat[6,5] <- NA

# built the model
sv <- cbind(diff_epu, diff_rgdp, diff_cpi, diff_u, diff_ffr, diff_mspus)
sv <- window(sv, start=c(1985, 2))

# lag order selection
lagselect <- VARselect(sv, lag.max = 8, type = "const")
lagselect$selection

# Estimate the model
Model1 <- VAR(sv, p = 3, season = NULL, exogen = NULL, type = "const")
SVARMod1 <- SVAR(Model1, Amat = amat, Bmat = NULL, hessian = TRUE, estmethod = c("scoring", "direct"))
SVARMod1

# impulse response function
# I. Generate All IRFs:
variables <- colnames(sv)
irf_list <- list()

for (impulse_var in variables) {
  for (response_var in variables) {
    current_irf <- irf(SVARMod1, impulse = impulse_var, response = response_var, n.ahead = 48, ortho = TRUE, 
                       cumulative = FALSE, boot = TRUE, ci = 0.95, runs = 500)
    irf_list[[paste(impulse_var, response_var, sep = "_to_")]] <- current_irf
  }
}

# II. Plotting All IRFs in a Grid:
# Convert IRFs to data.frames for ggplot
df_list <- lapply(names(irf_list), function(name) {
  df <- data.frame(Time = 1:length(irf_list[[name]]$irf[[1]]))
  
  df$Estimate <- irf_list[[name]]$irf[[1]]
  df$Lower <- irf_list[[name]]$Lower[[1]]
  df$Upper <- irf_list[[name]]$Upper[[1]]
  
  df$Impulse <- strsplit(name, "_to_")[[1]][1]
  df$Response <- strsplit(name, "_to_")[[1]][2]
  return(df)
})

# Bind all data frames in the list into one data frame
combined_df <- dplyr::bind_rows(df_list)

# Set the order for Facet
ordering <- c("diff_epu", "diff_rgdp", "diff_cpi", "diff_u", "diff_ffr", "diff_mspus")
combined_df$Impulse <- factor(combined_df$Impulse, levels = ordering)
combined_df$Response <- factor(combined_df$Response, levels = ordering)

combined_df$Facet <- paste(combined_df$Impulse, combined_df$Response, sep=" ~ ")
combined_df$Facet <- factor(combined_df$Facet, levels = facet_ordering)

# Plot using ggplot
ggplot(combined_df, aes(x = Time, y = Estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid", size = 0.3) +  # Adding the red line at y=0
  facet_wrap(~ Facet, scales = "free", ncol = length(ordering)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )