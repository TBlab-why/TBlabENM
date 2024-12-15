
#' @title Create formatted input data from the maxent results folder
#' @description 本函数从maxent的结果文件中提取物种名、T_MTSSte, T_MTSStr和AUC，作为进行
#' 物种丰富度、避难所等分析的输入文件。注意：重新运行该命令会覆盖上一次的结果。
#' @param resultdir maxent结果文件路径。
#' @param proname 字符型向量。要提取的时期名称。根据自己的时期确定。
#' @return 一个格式规范的.csv文件。
#' @export
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @examples
#' #设置工作路径
#' setwd("D:/Desktop")
#' #读取参数结果
#' maxent_spdata(resultdir="./TBlabENM/result",
#'           proname = "present"
#'  )
maxent_spdata <- function(
    resultdir = "./TBlabENM/result",
    proname = "present" ) {

  #检查输入变量
  if (file.exists(resultdir)==FALSE){
    stop("Resultdir not find.")}
##创建一个包含5列的空数据框
spdata <- data.frame(matrix(ncol = 5, nrow = length(list.files(resultdir))))
##重命名列名，分别
names(spdata) <- c("species", "occurrence", "T_MTSSte", "T_MTSStr", "auc")
###从结果文件读取阈值
for (i in seq_along(list.files(resultdir))) {
  #填充第一列"species"
  spdata[1] <- list.files(resultdir)
  #读取每个模拟结果的‘maxentResults.csv’文件
  data <- utils::read.csv(paste0(resultdir, "/", spdata[i,1], "/", proname, "/maxentResults.csv"), row.names=1)

  #提取"occurrence"
  occurrence <- data['species (average)','X.Training.samples']+data['species (average)','X.Test.samples']
  #将值赋给第二列"occurrence"
  spdata[i,2] <- occurrence

  #提取"MTSSte"
  MTSSte <- data['species (average)','Maximum.test.sensitivity.plus.specificity.Logistic.threshold']
  #将值赋给第3列"MTSSte"
  spdata[i,3] <- MTSSte

  #提取"MTSStr"
  MTSStr <- data['species (average)','Maximum.training.sensitivity.plus.specificity.Logistic.threshold']
  #将值赋给第4列"MTSSte"
  spdata[i,4] <- MTSStr

  #提取"Test.AUC"
  auc <- data['species (average)','Test.AUC']
  #将值赋给第5列"Test.AUC"
  spdata[i,5] <- auc
# Sun Oct  1 15:35:33 2023 ------------------------------
  #提取上级目录
  updir <- unlist(stringr::str_split(resultdir, "/"))[length(unlist(stringr::str_split(resultdir, "/")))-1]
  updir1 <- unlist(stringr::str_split(resultdir, "/"))[1:length(unlist(stringr::str_split(resultdir, "/")))-1]
  d <- c()
  for (i in seq_along(updir1)) {
     d1 <-  updir1[i]
     d <- paste0(d, d1, "/")
     }

  utils::write.csv(spdata, row.names=F, file = paste0(paste0(d, "spdata.csv")))
}}
