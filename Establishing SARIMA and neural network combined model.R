# 01整理数据 ------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(openxlsx)
library(TSA)
library(forcats)
library(tseries)

dat <- read.xlsx("1origindata/2009-2015.xlsx") %>% select(发病率) %>% 
  as_vector() %>% c(., dat1) %>% 
  ts(frequency = 12, start = c(2009, 01))

##编写一个评价函数##
pingjia <- function(x, y){
  mape <- mean(abs(x - y) / x * 100)
  mae <- mean(abs(x - y))
  rmse <- sum(abs(x - y) ^ 2/length(x)) ^ 0.5
  mer <- mean(abs(x - y)) / mean(x) * 100
  return(c(mape = mape, mae = mae, rmse = rmse, mer = mer))
}

# 02ARIMA建立 ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
load("01整理数据.RData")

##数据处理##
(tdat <- ts(dat[-c(c(length(dat) - 5):length(dat))], frequency = 12, start = c(2009, 1)))##拟合数据
(pdat <- ts(dat[c(c(length(dat) - 5):length(dat))], frequency = 12, start = c(2016, 7)))##预测数据

##forecast中评估差分次数函数
findfrequency(tdat)#找到时间序列的季节性周期，为12个月
ndiffs(tdat)#为控制长期趋势、周期性差分次数，为0
nsdiffs(tdat)#为序列存在明显季节性，季节性差分1次


#考虑到序列的季节性和周期性，对数据进行Box_Cox变换、一次季节差分并进行单位根检验
adf.test(diff(BoxCox(tdat, lambda = "auto"), lag = 12, differences = 1), k = 1)

##图1-绘制时间序列的自相关图（ACF）和偏自相关图（PACF）
ggtsdisplay(tdat, theme = theme_bw(), xlab = "years", ylab = "Incidence(1/10万)")

##利用自动建模函数建立ARIMA模型
myarima <- auto.arima(tdat, stepwise = T, trace = T, lambda = "auto", biasadj = T)
summary(myarima)#查看总体结果

##图2-残差的自相关和偏自相关图
(LB <- checkresiduals(myarima, xlab = "years", ylab = "Residuals", theme = theme_bw(), title = ""))#检查残差分布及Ljung-Box test

##预测2016年7-12月发病率
arimafor <- forecast(myarima, h = 6)

pingjiat_arima <- pingjia(myarima$x, myarima$fitted)##SARIMA拟合结果
pingjiaf_arima <- pingjia(pdat, arimafor$mean)##SARIMA预测结果

##表1-计算季节指数##
season <- order(decompose(tdat, type = "multiplicative")$seasonal[1:12])##求季节指数
bind_cols(decompose(tdat, type = "multiplicative")$seasonal[1:12], 
          season) %>% 
  set_names(c("季节指数", "赋值")) %>% 
  write.csv("result/季节指数.csv")##对季节指数结果进行整理导出
traindata1 <- tibble(x1 = c(myarima[["fitted"]], arimafor$mean), x2 = c(rep(season, ceiling(length(tdat)/12))), 
                     x3 = rep(rep(c(1, 2), each = 12), 4), x4 = c(tdat, pdat))##对ARIMA拟合及预测结果进行整理，为下一步分析

# 03-1支持向量机 ----------------------------------------------------------------
library(e1071)

##数据预处理
traindata <- traindata1[-c(c(nrow(traindata1) - 5):nrow(traindata1)), ] %>% 
  set_names(paste0("V", 1:4))
forcastdata <- traindata1[c(c(nrow(traindata1) - 5):nrow(traindata1)), ] %>% 
  set_names(paste0("V", 1:4))

##筛选参数##
obj_svm <- tune(svm, V4 ~ ., data = traindata,
                ranges = list(gamma = seq(0.1, 0.9, 0.2), cost = 10^seq(0, 6, 2)),
                tunecontrol = tune.control()
)

plot(obj_svm, col = "black", color.palette = function(n)grey.colors(n, start = 0.1, end = 0.9))##查看参数筛选结果
obj_svm[["best.parameters"]]##查看最优参数


##利用最优参数组合建立模型
train_svm <- svm(V4 ~ ., data = traindata, gamma = obj_svm[["best.parameters"]][, "gamma"], cost = obj_svm[["best.parameters"]][, "cost"])

##查看拟合值及拟合效果
fit_svm <- as.numeric(train_svm$fitted)


(pingjiat_svm <- pingjia(tdat, fit_svm))

##预测##
forcast_svm <- predict(train_svm, forcastdata[, 1:3])
(pingjiaf_svm <- pingjia(pdat, forcast_svm))##查看预测效果


# 03-2RBF -----------------------------------------------------------------
library(RSNNS)

##参数筛选
canshu_rbf <- matrix(NA, 10, 10)
inputs <- normalizeData(as.matrix(traindata[1:3]))
inputsnor <- getNormParameters(inputs)
outputs <- normalizeData(traindata[[4]])
ouputsnor <- getNormParameters(outputs)

for (i in seq(0.01, 0.1, 0.01)) {
  for (j in seq(0.01, 0.1, 0.01)) {
    set.seed(9)
    mod <- rbf(inputs, outputs, size = 30, maxit = 1000, initFuncParams = c(0, 1, 0, i, j))
    fits <- as.numeric(denormalizeData(fitted(mod), getNormParameters(outputs)))
    canshu_rbf[i*100, j*100] <- pingjia(traindata[[4]], fits)["rmse"]
  }
}
which(canshu_rbf == min(canshu_rbf), arr.ind = T)
# row col
# 9  2

##利用最优参数组合建立RBF模型##
set.seed(9)
obj_rbf <- rbf(inputs, outputs, size = 30, maxit = 1000, initFuncParams = c(0, 1, 0, 0.09, 0.02))
##计算拟合值及拟合效果##
fit_rbf <- as.numeric(denormalizeData(fitted(mod), getNormParameters(outputs)))
(pingjiat_rbf <- pingjia(tdat, fit_rbf))

##预测##
forcast_rbf <- denormalizeData(predict(obj_rbf, scale(as.matrix(forcastdata[1:3]), 
                                                      center = inputsnor[["colMeans"]], scale = inputsnor[["colSds"]])), getNormParameters(outputs))
(pingjiaf_rbf <- pingjia(pdat, forcast_rbf))#查看预测效果


# 03-3GRNN建立 ----------------------------------------------------------------
library(yager)

##参数筛选##
inputs <- normalizeData(as.matrix(traindata[1:3]))
inputsnor <- getNormParameters(inputs)
outputs <- normalizeData(traindata[[4]])
ouputsnor <- getNormParameters(outputs)
canshu_grnn <- vector(length = 200)
for (i in seq(0.01, 2, 0.01)) {
  obj_grnn <- grnn.fit(inputs, outputs[, 1], sigma = i)
  fit_grnn <- denormalizeData(grnn.predict(obj_grnn, inputs), getNormParameters(outputs))
  (canshu_grnn[c(i*100)] <- pingjia(tdat, fit_grnn)["rmse"])
}
#查看参数筛选结果
plot(1:200, canshu_grnn, type = "l")
which(canshu_grnn==min(canshu_grnn), arr.ind = T)

##利用最优参数创建GRNN模型##
obj_grnn <- grnn.fit(inputs, outputs[, 1], sigma = .01)
##计算拟合值及拟合效果
fit_grnn <- as.numeric(denormalizeData(grnn.predict(obj_grnn, inputs), getNormParameters(outputs)))
(pingjiat_grnn <- pingjia(tdat, fit_grnn))

##预测##
fdata_grnn <- matrix(as.numeric(scale(as.matrix(forcastdata[1:3]), 
                                      center = inputsnor[["colMeans"]], scale = inputsnor[["colSds"]])), 6, 3)
forcast_grnn <- denormalizeData(grnn.predict(obj_grnn, fdata_grnn), getNormParameters(outputs))
(pingjiaf_grnn <- pingjia(pdat, forcast_grnn))


# 05-结果整理 -----------------------------------------------------------------


##表2-模型拟合及预测结果整理##
model_accuracy <- bind_rows(pingjiat_arima, pingjiat_svm, pingjiat_rbf, pingjiat_grnn) %>% 
  add_column(methods = c("SARIMA", "SVM", "RBF", "GRNN"), .before = "mape") %>% 
  bind_cols(bind_rows(pingjiaf_arima, pingjiaf_svm, pingjiaf_rbf, pingjiaf_grnn))
write_csv(model_accuracy, "result/model_accuracy.csv")

##图3-SVM参数筛选结果整理##
plot(obj_svm, col = "black", color.palette = function(n)grey.colors(n, start = 0.1, end = 0.9))##查看参数筛选结果

##图4-RBF参数筛选结果整理##
plotdat <- canshu_rbf %>% as_tibble() %>% set_names(as.character(1:10)) %>% 
  mutate(i = c(1:10) * 0.01) %>% 
  gather(key = "j", value = "value", -i) %>% mutate(j = as.integer(j) * 0.01, rvalue = cut(log(value), breaks = c(min(log(value))-0.1, min(log(value)), 1:7)))

library(ggThemeAssist)
library(scales)
cols <- gray.colors(6, 1, 0.01, 1)
p <- ggplot(plotdat, aes(i, j, fill = rvalue)) + geom_tile(color = "black") + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, 1) * 0.01, name = "deviation") + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, 1) * 0.01, name = "bias") + 
  scale_fill_manual(values = cols, name = "log(RMSE)")
p

##图5-GRNN参数筛选结果整理##
par(mgp = c(2, 0.7, 0))
plot(seq(0.01, 2, 0.01), canshu_grnn, type = "b", pch = 20, xlab = "sigma", ylab = "RMSE")





































