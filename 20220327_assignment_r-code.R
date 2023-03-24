# Load necessary libraries after installation
library(vars)
library(tseries)
library(tidyverse)
library(stargazer)
library(ggfortify)
library(dplyr)
library(forecast)
library(dplyr)
library(data.table)
library(stats)

# Load data into R
y = read.csv("C:/Users/Eric/OneDrive - University of Warwick/Economics (Hons) BSc/Year 4/EC306_Econometrics-2-Time-Series/Assignments/Assignment 1/GDP_consumption_altered.csv")

# Create credit spread var long-term - short-term (Ten_yr - T3M)
y$Cre_Spr <- y$Ten_yr - y$T3M
y$Cre_Spr <- abs(y$Cre_Spr) # Cannot be negative, absolute value of (Ten_yr - T3M)

# Rename vars using base R
names(y)[names(y)=="Consumption"] <- "Con"
names(y)[names(y)=="Investment.PNFI"] <- "Inv"

'The dataset contains 6 variables: GDP (GDP), Consumption (Con), Investment (Inv)
3 month treasury bills (T3M), 10 year government bond rate (Ten_yr), Moodys seasoned BAA bond rate (BAA) and credit spread (Cre_Spr)'

# Inspect TS
head(y) # first few obs
tail(y) # last few obs
str(y) # structure of data frame


# Define vars as TS object with start + end date
TsGDP <- ts(y$GDP, start = c(1959,1), frequency = 4)
TsCon <- ts(y$Con, start = c(1959,1), frequency = 4)
TsInv <- ts(y$Inv, start = c(1959,1), frequency = 4)
TsT3M <- ts(y$T3M, start = c(1959,1), frequency = 4)
TsBAA <- ts(y$BAA, start = c(1959,1), frequency = 4)
TsTen_yr <- ts(y$Ten_yr, start = c(1959,1), frequency = 4)
TsCre_Spr <- ts(y$Cre_Spr, start = c(1959,1), frequency = 4)

# Plot data
autoplot(TsGDP) # trend but no mean, non-stationary
#take log of GDP
#lGDP <- log(GDP)
autoplot(log(TsCon), main = "Log of consumption over time", ylab = "Log of consumption")
autoplot(log(TsInv), main = "Log of investment over time", ylab = "Log of investment") # trend + drift until late 2019, then sudden drop, overall not stationary
autoplot(T3M) 'Slightly volatile until 1980, non-zero mean with trend; from 1980s downward trend + non-zero mean,
looks like random walk, overall non-zero mean + non-stationary'
autoplot(TsBAA, main = "Moody's seasoned bond rate BAA") # Upward trend until 1980, then downward trend, looks like random-walk, non-zero mean overall
autoplot(Ten_yr) # Upward trend until 1985, there spike, then downward trend, random walk with drift?
autoplot(TsCre_Spr, main = "First difference of credit spread over time", ylab = "First difference of absolute credit spread") # Highly volatile, non-zero mean (around 1.5)

autoplot(TsCon, main="Consumption over time")
# Plot the first difference to make it stationary

'Preparing the data:
exchange rates usually follow a random walk (ie unit root). Thus, we will do our estimations using differenced data'



####Question 1: VAR -> credit spreads on I; C, I, and credit spread , 1965####

# Subselect timespan sample 1965q1-2006 q4

# Select subsample from 1965q1 to 2006q4
# Select variables consumption, 

quarters <- seq(ISOdate(1959,1,1), length.out=243, by = "quarter")


y1 <- as_tibble(y)
y1 %<>% add_column(t = quarters) %>% 
  dplyr::select(t, Con, Inv, Cre_Spr) %>%
  filter(t >= ISOdate(1965,1,1), t <= ISOdate(2006,12,31)) %>%
  dplyr::select(Con, Inv, Cre_Spr)

y1 <- as.data.frame(y1)

TsCon1 <- ts(y1$Con, start = c(1965,1), frequency = 4)
TsInv1 <- ts(y1$Inv, start = c(1965,1), frequency = 4)
TsCre_Spr1 <- ts(y1$Cre_Spr, start = c(1965,1), frequency = 4)

#plot(TsCon1, main = "Time Series consumption")



# Plot the first difference to make it stationary
TsCon1D <- diff(log(TsCon1), differences = 1) # Non-zero mean, stationary
plot(TsCon1D)
TsInv1D <- diff(log(TsInv1), differences = 1) # Non-zero mean, stationary
plot(TsInv1D)
plot(TsCre_Spr1)# Non-zero mean, stationary

# Check that series have a unit root, stationary I(0) with ACF and PACf

# ADF test from tseries; rest of functions from vars
# H0: Data is not stationary, Has a unit root
# H1: Otherwise, i.e. stationary
adfCon1 <- adf.test(TsCon1D)
adf.test(TsCon1D) # Reject H0. Stationary.
adfInv1 <- adf.test(TsInv1D)
adf.test(TsInv1D) # Not reject H0. Has a unit root. Lag order 5.
adfCre_Spr1 <- adf.test(TsCre_Spr1)
adf.test(TsCre_Spr1) # Reject H0. Stationary. P-value = 0.03.



y1 <- as_tibble(y)
y1 %<>% add_column(t = quarters) %>% 
  dplyr::select(t, Con, Inv, Cre_Spr) %>%
  filter(t >= ISOdate(1965,1,1), t <= ISOdate(2006,12,31)) %>%
  dplyr::select(Con, Inv, Cre_Spr)

# Overwrite vars in y1

y1[,c(1,2)] <- log(y1[,c(1,2)])
# now convert into datatable for diff() function
dy1 = data.table(y1)
y1 <- dy1[, .(Con1D = diff(Con, differences = 1),Inv1D = diff(Inv, differences = 1), Cre_Spr)]


# Specify the optimal lag length

lag <- VARselect(y1) # Use the dataframe here instead of TS object
lag$selection
# Akaike Information Criterion chooses 3, most lags of all methods HQ, SC, and FPE
# loss of power overall but not relevant with t = 244;

# Check order of co-integration Johanson's co-integration
# Check co-integration between Cre_Spr and Inv
# H0: no cointegration on specific rank 
# Ha: maximum number of cointegrating equations at this rank. 
jo_trace1 <- ca.jo(y1, type = "trace", ecdet = "constant", season = 4, K = 2)
summary(jo_trace1)


#TODO: build linear stationary combination

# Estimate VAR(3)
estim1 <- VAR(y1, type = "const", lag.max = 3, ic = c("AIC"), seasonal) # p=3 is optimal lag length, 3 lags (168-3)
summary(estim1)

stargazer(estim1[["varresult"]], type="text") # Like in a publication, use for report

# problem if different trends etc..., check in theory



# Check that model is stable (eigenvalues < 1), if TS stationary then likely < 1
roots(estim1, modulus = TRUE) # I(1), one unit root

# Check number of differences required to achieve stationarity
ndiffs(y1) #in forecast package

# Series have to be stationary
# Granger causality test
grangerCre_Spr <- causality(estim1, cause = "Cre_Spr")
grangerCre_Spr$Granger # Very small p-value, reject null hypothesis. Cre_Spr granger causes rest

grangerInv <- causality(estim1, cause = "Inv1D") # Very small p-value, reject null hypothesis. Inv granger causes rest
grangerInv$Granger

####Question 2a: IRFs credit spreads on investment####

# IRFs
# effect of a shock to credit spreads on investment.
# n.ahead on how many time periods ahead
# boot := bootstrapped error bands for IRFs
irf1 <- irf(estim1, impulse = "Cre_Spr", response = "Inv",
            n.ahead = 20, boot = TRUE, runs = 200, ci = 0.95)
plot(irf1, ylab = "Investment",
     main = "Percentage change of investment response to credit spreads shock")
# Dying out

irf2 <- irf(estim1, impulse = "Inv", response = "Cre_Spr",
            n.ahead = 20, boot = TRUE, runs = 200, ci = 0.95)
plot(irf2, ylab = "Credit Spreads",
     main = "Credit spreads response to investment shock")
# Not dying out

# Var decomposition
# check on how much the shock is due to variables themselves -> breakdown
# tells percentage of shock from var
vd <- fevd(estim, n.ahead = 5)
plot(vd)

####Question 2b: Additional info for credit spread -> Moody's seasoned BAA bond rate####

# Update credit spread?

cor.test(y$Cre_Spr, y$BAA, method = "pearson", data = y)

####Question 3a: Estimate VAR(4) -> C, I, short rate, credit spread 1965q1 to 2000q4####
###compare with ARMA model for C###

# assume trend stationarity -> difference once

quarters <- seq(ISOdate(1959,1,1), length.out=243, by = "quarter")

y3 <- as_tibble(y)
y3 %<>% add_column(t = quarters) %>% 
  dplyr::select(t, Con, Inv, Cre_Spr, T3M) %>%
  filter(t >= ISOdate(1965,1,1), t <= ISOdate(2000,12,31)) %>%
  dplyr::select(Con, Inv, Cre_Spr, T3M)

y3 <- as.data.frame(y3)

# Overwrite vars in y3

y3[,c(1,2)] <- log(y3[,c(1,2)])
plot(y3$Con, ylab = "Percentage change in consumption")
# now convert into datatable for diff() function
#dy3 = data.table(y3)
#y3 <- dy3[, .(Con1D = diff(Con, differences = 1),Inv1D = diff(Inv, differences = 1), Cre_Spr, T3M)]

# Estimate VAR(4)
estim3 <- VAR(y3, p = 4, type = "trend") # p=3 is optimal lag length, 3 lags (168-3)
summary(estim3)

stargazer(estim3[["varresult"]], type="text") # Like in a publication, use for report

# problem if different trends etc..., check in theory

# Check co-integration between Cre_Spr and Inv
# H0: no cointegration on specific rank 
# Ha: maximum number of cointegrating equations at this rank. 

# Check that model is stable (eigenvalues < 1), if TS stationary then likely < 1
roots(estim, modulus = TRUE) 

# Check number of differences required to achieve stationarity
ndiffs(y1) #in forecast package

# Series have to be stationary
# Granger causality test
grangerCre_Spr <- causality(estim, cause = "Cre_Spr")
grangerCre_Spr$Granger # Very small p-value, reject null hypothesis. Cre_Spr granger causes rest

grangerInv <- causality(estim, cause = "Inv") # Very small p-value, reject null hypothesis. Inv granger causes rest
grangerInv$Granger

##Out of sample 1-step forecasting data: 2001Q1-2006Q4
quarters <- seq(ISOdate(1959,1,1), length.out=243, by = "quarter")

y31 <- as_tibble(y)
y31 %<>% add_column(t = quarters) %>% 
  dplyr::select(t, Con, Inv, Cre_Spr, T3M) %>%
  filter(t >= ISOdate(2001,1,1), t <= ISOdate(2006,12,31)) %>%
  dplyr::select(Con, Inv, Cre_Spr, T3M)

y31 <- as.data.frame(y31)

# Overwrite vars in y31
y31[,c(1,2)] <- log(y31[,c(1,2)])

##ARMA training 1 step ahead: 1965Q1 - 2000Q4 then test on 
fit2VAR <- predict(estim3, n.ahead=1, cl=0.95)
plot(fit2VAR)

# 1. Model identification and selection: Assume stationarity
Acf(y3$Con, ci = 0.95, type = c("correlation"), plot = TRUE)
Pacf(y3$Con, ci = 0.95, type = c("partial"), plot = TRUE)
# Seasonality, differencing -> not required

# 2. MLE estimation
fitARMA <- auto.arima(y3[,1])
fitARMA

# 3. Statistical model checking
checkresiduals(y2, test="LB" )

# Forecast ARIMA
fit2ARMA <- Arima(y31[,1], model=fitARMA)
onestepARMA <- forecast(fit2ARMA, h=1)
fit2ARMA
plot(onestepARMA)
####3b: Compare ARMA in 3a with dynamic start 2001, 2006, 2012, better dynamics conditional on estimating up to 2000####

quarters <- seq(ISOdate(1959,1,1), length.out=243, by = "quarter")

y32 <- as_tibble(y)
y32 %<>% add_column(t = quarters) %>% 
  dplyr::select(t, Con, Inv, Cre_Spr, T3M) %>%
  filter(t >= ISOdate(2001,1,1), t <= ISOdate(2019,12,31)) %>%
  dplyr::select(Con, Inv, Cre_Spr, T3M)

y32 <- as.data.frame(y32)

y32[,c(1,2)] <- log(y32[,c(1,2)])

##2001Q1 - 2019Q4 dynamic forecast



####Question 4: arch effects for I growth rate?####

# Check Cholesky decomposition, test for ARCH effect

arch.test(fit2VAR, arch = c("box"), output = TRUE, alpha = 0.05)
