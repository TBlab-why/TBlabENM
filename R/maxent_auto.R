
#' Maxent Automation:From occurrence to maps of suitable areas
#' @description
#' 本函数自动执行MaxEnt模型，当设置好参数后，首先进行相关参数的选择，包括环境变量、
#'     正则化乘数和特征函数。接着使用选择好的最佳参数进行建模然后投影至指定的时空生成
#'     适宜性地图。
#'
#' @param spdir 要模拟的一或多个物种文件的路径，包括文件名（.csv）。
#' @param evdir 环境变量路径,格式为.asc
#' @param evlist 用于建模的环境变量，是环境变量路径下所有变量的子集。用下标表示，
#' 数字代表要选用的变量，默认使用所有变量。
#' @param factors 字符型向量。指定哪些变量是分类变量。未指定时，则默认全为连续变量。
#' @param nbg 随机背景点数量，当指定参数mybgfile时忽略。
#' @param mybgfile 自定义的背景数据，包含两列（经度、纬度）
#' @param args 自定义的MaxEnt模型参数，详见\code{\link[TBlabENM]{maxent_args}}。
#' @param myenv 指定的一组环境变量，当提供此参数时，将不会再根据相关性进行变量筛选
#' @param fc 特征函数,取值为L,Q,P,H,T,及它们的组合
#' @param rm 正则化乘数 数值型
#' @param r 相关性系数
#' @param cormethod 计算相关性的方法，取值为"pearson" (default), "kendall", or "spearman"
#' @param vif 逻辑型，当为T时对所有候选变量组合进行方差膨胀因子检验，排除具有强共线性的组合
#' @param vifth 方差膨胀因子阈值
#' @param opt 选择最佳模型的指标，可选择"aicc","seq","cbi".
#' @param prodir 用列表储存的投影文件路径，包含了投影时期的名称和与之对应的环境变量路径
#' @param outdir 结果保存路径
#' @param parallel 是否并行计算
#' @param ncpu 如果并行，使用的cpu数
#'
#' @return 包含maxent相关结果的文件夹
#' @export
#'
#' @examples
#' maxent_auto(sp = system.file("extdata", "species/Phoebe sheareri.csv", package = "TBlabENM"),
#' evdir = system.file("extdata", "envar/asc", package = "TBlabENM"),
#' evlist = 1:7,
#' factors = c("dr", "fao90"),
#' nbg = 1000,
#' args = maxent_args(),
#' myenv = NULL,
#' fc = c("lph", "q", "lq"),
#' rm = 1:2,
#' r = 0.7,
#' vif = T,
#' vifth = 5,
#' opt = "aicc",
#' outdir = NULL)
#'
maxent_auto <- function(spdir, evdir, myenv = NULL, evlist = NULL, factors = NULL,
                        mybgfile = NULL, nbg = 10000, args = maxent_args(),
                        fc, rm, r = 0.7, cormethod = "pearson", vif = T, vifth = 5,
                        opt = NULL, prodir = NULL, outdir = NULL, parallel = F, ncpu = 2){

  fun3 <- function(x){
    cat("*****************************************************\n")
    pa <- TBlabENM::maxent_parameter(x = x,
                    evdir = evdir,
                    evlist = evlist,
                    factors = factors,
                    nbg = nbg,
                    fc = fc,
                    rm = rm,
                    r = r,
                    vif = vif,
                    vifth = vifth,
                    opt = opt,
                    outdir = outdir,
                    parallel = parallel,
                    ncpu = ncpu)
  print("将使用以下参数构建最终模型")
  print(pa)
    #模拟
    ##设置args参数
    args[3] <- "linear=FALSE"
    args[4] <- "quadratic=FALSE"
    args[5] <- "product=FALSE"
    args[6] <- "threshold=FALSE"
    args[7] <- "hinge=FALSE"

    args[2] <- paste0("betamultiplier=", pa$rm)
    ff1 <- stringr::str_split_1(pa$fc, "")
    for (j in ff1) {
      if(j == "L"){args[3] <- "linear=TRUE"}
      if(j == "Q"){args[4] <- "quadratic=TRUE"}
      if(j == "P"){args[5] <- "product=TRUE"}
      if(j == "T"){args[6] <- "threshold=TRUE"}
      if(j == "H"){args[7] <- "hinge=TRUE"}
    }
    bio_name <- stringr::str_split_1(pa$env, ",")
    biolistall <- list.files(evdir, pattern = ".asc$", full.names = TRUE)
    evlist <- c()
    for (i in seq_along(bio_name)) {
      evlist1 <- which( stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i])==T)
      evlist <- c(evlist, evlist1)
    }

    #模拟
    cat("*****************modelling*****************\n")
    ms <- TBlabENM::maxent_single(
      x = x,
      evdir = evdir,
      evlist = evlist,
      factors = factors,
      nbg = nbg,
      mybgfile = mybgfile,
      args = args,
      prodir = prodir,
      outdir = outdir,
      parallel = parallel,
      ncpu = ncpu
    )

    }
  maxent_args <- function(replicates = 10,
                          betamultiplier = 1,
                          l = TRUE,
                          q = TRUE,
                          p = FALSE,
                          t = FALSE,
                          h = FALSE,
                          replicatetype = "crossvalidate",
                          responsecurves = TRUE,
                          jackknife = TRUE,
                          pictures = TRUE,
                          outputgrids = FALSE,
                          outputformat = "logistic"
  ) {
    c(
      paste0("replicates=", replicates),
      paste0("betamultiplier=", betamultiplier) , #重复次数和正则化乘数
      paste0("linear=", l),    # 5种特征函数
      paste0("quadratic=", q),
      paste0("product=", p),
      paste0("threshold=", t),
      paste0("hinge=", h),
      paste0("replicatetype=", replicatetype),  #重复类型
      paste0("responsecurves=", responsecurves),   #响应曲线
      paste0("jackknife=", jackknife),     #折刀分析
      paste0("pictures=", pictures),
      paste0("outputgrids=", outputgrids),
      paste0("outputformat=", outputformat) #输出文件格式
    )
  }
  #####################################################
  #新建文件夹
  if(is.null(outdir)){dir.create("./TBlabENM", showWarnings = FALSE)} else{
    dir.create(paste0(outdir, "/TBlabENM"), showWarnings = FALSE)}

  star_time <- Sys.time() ## 记录程序开始时间
  if(parallel == T){
    # library(snowfall)
    # 开启集成
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    # 注册每个环境变量
    snowfall::sfExport("fun3")
    snowfall::sfExport("maxent_args")
   # snowfall::sfLibrary(TBlabENM)
    snowfall::sfExport("spdir")
    snowfall::sfLapply(spdir, fun3)
    snowfall::sfStop()  # 关闭集群

  } else{

    ## 第一个位置：新建起始进度条
    pb <- utils::txtProgressBar(style=3)
    #新建向量保存失败的模型
    failed_species <- c()
    for (x in spdir) {
      fit <- try(fun3(x))
      if('try-error' %in% class(fit)){
        failed <- paste0(x, " was failed!")
      failed_species <- c(failed_species, failed)
      next}
    }

#   for (i in seq_along(spdir)) {
#   x=spdir[i]
#   fun3(x)
#
# }
    ##第二个位置：实时显示进度
    utils::setTxtProgressBar(pb, which(x==spdir)/length(spdir))

  ## 第三个位置关闭进度条
  if(parallel == F){close(pb)
    if(is.null(failed_species) == FALSE) {print(failed_species)}}
  }
  #读取结果文件返回相关参数
  df <- data.frame(matrix(NA,0,10))
  names(df) <- c("species", "number", "env", "fc", "rm", "replicates", "AUCtrain", "AUCtest", "MTSStrain", "MTSStest")
  #提取物种名
  nm <- c()
  for(i in spdir){
    spname1 <- stringr::str_split_1(i, "/")[length(stringr::str_split_1(i, "/"))]
    sp_name <- stringr::str_split_1(spname1, ".csv$")[1]
    nm <- c(nm, sp_name)
  }
  for (i in nm) {
      df1 <- utils::read.csv(paste0(outdir,"/TBlabENM/maxent/", i, "/parameters.csv"))
    df <- rbind(df,df1)
  }

  utils::write.csv(df, paste0(outdir,"/TBlabENM/maxent/allsp_parameters.csv"),row.names = FALSE)

  end_time <- Sys.time()  ## 记录程序结束时间
  print(end_time - star_time)
  cat("****************completion*****************\n")
    return(df)
}

