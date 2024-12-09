
#' Select optimum parameters for the MaxEnt models
#' @description
#' 本函数为MaxEnt模型选择最适合的环境变量组合、正则化乘数和特征函数组合。
#' @details
#' 工作流程：首先对fc和rm进行所有可能的组合，对于每个组合，使用全部变量进行建模获得变量贡献大小，
#'     依据变量贡献大小排除与最重要变量具有强相关性的变量。然后使用保留的变量继续建模获得该组变量
#'     贡献大小，依据变量贡献大小排除与第二重要的变量具有强相关性的变量。重复上述过程直到保留的
#'     变量之间不具有强相关性。经过此步骤就为每个fc和rm组合获得了一组不具有强相关性且保留了贡献大
#'     的候选变量，当vif=T时，还对这一系列候选变量组合进行方差膨胀因子分析，进一步排除具有强共线性
#'     的组合。最后使用ENMevaluate对每个组合进行测试，获得相应的指标（例如“aicc”），根据指标选择
#'     最佳的一组组合作为模型最终参数。
#'  注意：当物种发生数据过少时，可能无法得出正确的AICC和CBI值，如果opt选择aicc或cbi，则使用
#'     auc.val.avg值进行代替。
#'
#' @param x 要模拟的物种文件的路径，包括文件名（.csv）。
#' @param evdir 环境变量路径,格式为.asc
#' @param evlist 用于建模的环境变量，是环境变量路径下所有变量的子集。用下标表示，
#' 数字代表要选用的变量，默认使用所有变量。
#'
#' @param factors 字符型向量。指定哪些变量是分类变量。未指定时，则默认全为连续变量。
#' @param fc 特征函数,取值为L,Q,P,H,T,及它们的组合
#' @param rm 正则化乘数 数值型
#' @param r 相关性系数
#' @param mybgfile 自定义的背景数据，包含两列（经度、纬度）
#' @param cormethod 计算相关性的方法，取值为"pearson" (default), "kendall", or "spearman"
#' @param vif 逻辑型，当为T时对所有候选变量组合进行方差膨胀因子检验，排除具有强共线性的组合
#' @param vifth 方差膨胀因子阈值
#' @param opt 选择最佳模型的指标，可选择"aicc","seq","cbi"或NULL.当选择NULL时，则仅选择变量。
#' @param outdir 结果保存路径
#' @param nbg 随机背景点数量，当指定参数mybgfile时忽略。
#' @param myenv 指定的一组环境变量，当提供此参数时，将不会再根据相关性进行变量筛选
#' @param parallel 是否并行计算
#' @param ncpu 如果并行，使用的cpu数
#'
#' @importFrom dplyr mutate arrange select
#' @importFrom stringr str_split_1
#' @importFrom utils read.csv tail write.csv
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr map map_dbl
#'
#' @return 生成一个csv文件，包含了每种fc、rm、选择的环境变量组合的评估结果
#' @export
#'
#' @examples
#' maxent_parameter(x = system.file("extdata", "species/Phoebe sheareri.csv", package = "TBlabENM"),
#'        evdir = system.file("extdata", "envar/asc", package = "TBlabENM"),
#'        myenv = NULL,
#'        evlist = 1:7,
#'        factors = c("dr", "fao90"),
#'        nbg = 1000,
#'        fc = c("LPH", "Q", "LQ"),
#'        rm = 1:5,
#'        r = 0.7,
#'        vif = T,
#'        vifth = 5,
#'        opt = "aicc",
#'        parallel = F,
#'        ncpu = 2)

maxent_parameter <- function(x, evdir, myenv = NULL, evlist = NULL, factors = NULL,
                             mybgfile = NULL, nbg = 10000, fc, rm, r = 0.7,
                             cormethod = "pearson", vif = T, vifth = 5, opt = "aicc",
                             outdir = NULL, parallel = F, ncpu = 2){
  # unlink(paste0(outdir, "/TBlabENMtemp*") ,recursive = T)

  star_time <- sample(1:100000,1)

#corse_method功能使用相关性选择变量
  corse_method <- function(correlation, importance, vif, n){ #n最初为1

    n = n
    #这里可能会出现变量重要性不一致问题，之前选出来的变量直接保留占位，其他变量排序
    if(n==1){importance1 <- dplyr::arrange(importance, desc(value))} else {
    importance_a <- importance[1:(n-1), , drop = FALSE]
    importance_b <- dplyr::arrange(importance[-(1:(n-1)), , drop = FALSE], desc(value))
    importance1 <- rbind(importance_a, importance_b)}
    #importance1 <- dplyr::arrange(importance, desc(value))
    if (nrow(importance) < n) {importance <- rownames(importance1)[1]} else{
      importance <- rownames(importance1)[n]}

    ff <- correlation[,importance, drop = FALSE]
    #ff <- dplyr::filter(ff, rowSums(ff < r) == ncol(ff))
    ff <-  ff[which(ff < r), ,drop = FALSE]

    bio_name <- c(importance, rownames(ff)) #要保留的变量名称


    #当执行vif分时
    if(vif==T){

      bio_name <- vifse_method(envdf = mybg, env_name = bio_name,
                                         importance = importance1, vifth = vifth, n = n )
      #importance继承自importance1，importance1是排过序的，因此importance不用再排序
    }

    return(bio_name)
  }
  #vifse_method功能使用vif选择变量
  vifse_method <- function(envdf, env_name, importance, vifth, n){
    #对保留的变量进行重要性排序
    #im <- dplyr::arrange(importance[env_name,, drop =FALSE],desc(value))
    #env_name <- rownames(im) #排过序的变量名
    #importance中取出保留的变量
    im <- importance[rownames(importance) %in% env_name, , drop = FALSE]
    env_name <- rownames(im) #这里的env_name是排过序的变量名
    env_name1 <- env_name
    n = n-1

    while (n < length(env_name1)-1) {

      n = n+1
      #从env_name中取出前n+1个变量，判断最重要的变量和次重要的变量之间的共线性，如果存在，次重要的变量删除
      var <- env_name1[c(1:(n+1))]
      #只保留在env_name1中的变量
      var1 <- var[var%in%env_name]
      #对var1计算vif
      vifvalue <- usdm::vif(envdf[var1])

      #当存在共线性则应该排除var1中的最后一个变量,当不存在共线性说明该变量要保留，进入下次循环
      if(sum(vifvalue$VIF > vifth)>0){env_name <- env_name[!env_name%in%var1[length(var1)]]} else{break}
    }

    return(env_name) }


  fun2 <- function(fc1, rm1){

    #根据组合设置 args参数
    args1[2] <- paste0("betamultiplier=", rm1)
    ff1 <- stringr::str_split_1(fc1, "")
    for (j in ff1) {
      if (j == "L") {
        args1[3] <- "linear=TRUE"
      }
      if (j == "Q") {
        args1[4] <- "quadratic=TRUE"
      }
      if (j == "P") {
        args1[5] <- "product=TRUE"
      }
      if (j == "T") {
        args1[6] <- "threshold=TRUE"
      }
      if (j == "H") {
        args1[7] <- "hinge=TRUE"
      }
    }

    n = 0 #用于fun1(),指示第几次循环来确定次重要变量
    if(is.null(factors) == FALSE) {
      #当存在分类变量，分开计算
      #当factors只有1个变量
      if(length(factors)==1){
        while(sum(correlation1 > r) > nrow(correlation1) | n == 0){
          n <- n + 1

          #进行模拟获得重要性列表
          compare <- suppressMessages(maxent_single(
            x = x,
            evdir = evdir,
            evlist = evlist, #evlist应该是相对于evdir也就是全部变量biolstall的下标
            factors = factors,
            nbg = nbg,
            mybgfile = mybgfile,
            args = args1,
            prodir = NULL,
            outdir = paste0(outdir, "/TBlabENMtemp",star_time,"/", fc1,rm1,n),
            parallel = FALSE))
          #下面根据上面模拟的结果删除相关性强的变量
          #变量重要性
          ev_cb <- read.csv(paste0(outdir, "/TBlabENMtemp", star_time, "/", fc1, rm1, n, "/TBlabENM/maxent/", sp_name,"/maxentResults.csv")) %>%
            dplyr::select(., paste0(bio_name, ".contribution")) %>%
            utils::tail(., n =1) %>% t() %>% as.data.frame()
          names(ev_cb) <- "value"

          ##排除贡献小于0.5的变量
          ev_cb <- ev_cb[ev_cb>0.5,, drop = FALSE]

          #修改ev_cb的变量名
          nn <- c()
          for (i in 1:nrow(ev_cb)) {
            nn1 <- stringr::str_split_1(rownames(ev_cb)[i], ".contribution")[1]
            nn <- c(nn, nn1)
          }
          rownames(ev_cb) <- nn
          #将ev_cb按照变量类型拆开
          factors <- rownames(ev_cb)[rownames(ev_cb) %in% factors] #更新factors
          if (length(factors) == 0) {factors <- NULL}

          if(is.null(factors)){ev_cb1 <- ev_cb} else {
            ev_cb1 <- ev_cb[!rownames(ev_cb) %in% factors, , drop = FALSE]}


          #correlation1按照ev_cb提取
          correlation1 <- correlation1[rownames(ev_cb1),rownames(ev_cb1)]

          #连续变量的选择n=1
          bio_name1 <- corse_method(correlation = correlation1, importance = ev_cb1, vif=vif, n=n) #下次模拟的变量
          #组合两种变量
          bio_name <- c(factors, bio_name1)
          #获取保留的变量存储路径的下标
          evlist <- c()
          for (i in seq_along(bio_name)) {
            evlist1 <- which( stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i])==T)
            evlist <- c(evlist, evlist1)
          }
          #计算两种变量的相关性
          correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
        }

        return(bio_name)

      } else { #分类变量大于1个的情况

        while(sum(correlation1 > r) > nrow(correlation1) | sum(correlation2 > r) > nrow(correlation2)| n==0){

        n <- n+1

          #进行模拟获得重要性列表
          compare <-suppressMessages(maxent_single(
            x = x,
            evdir = evdir,
            evlist = evlist,
            factors = factors,
            nbg = nbg,
            mybgfile = mybgfile,
            args = args1,
            prodir = NULL,
            outdir = paste0(outdir, "/TBlabENMtemp", star_time, "/", fc1, rm1, n),
            parallel = FALSE))
          #下面根据上面模拟的结果删除相关性强的变量
          #变量重要性
          ev_cb <- read.csv(paste0(outdir, "/TBlabENMtemp",star_time,"/",fc1,rm1,n, "/TBlabENM/maxent/", sp_name, "/maxentResults.csv")) %>%
            dplyr::select(., paste0(bio_name, ".contribution")) %>%
            utils::tail(., n =1) %>% t() %>% as.data.frame()
          names(ev_cb) <- "value"

          ##排除贡献小于0.5的变量
          ev_cb <- ev_cb[ev_cb>0.5,, drop = FALSE]

          #修改ev_cb的变量名
          nn <- c()
          for (i in 1:nrow(ev_cb)) {
            nn1 <- stringr::str_split_1(rownames(ev_cb)[i], ".contribution")[1]
            nn <- c(nn, nn1)
          }
          rownames(ev_cb) <- nn

          factors <- rownames(ev_cb)[rownames(ev_cb) %in% factors] #更新factors
          if (length(factors) == 0) {factors <- NULL}

          #将ev_cb按照变量类型拆开
          if (is.null(factors) ) {ev_cb1 <- ev_cb} else {
            ev_cb1 <- ev_cb[!rownames(ev_cb) %in% factors, , drop = FALSE]
            ev_cb2 <- ev_cb[factors, , drop = FALSE]
          }

          #correlation1按照ev_cb提取
          if (is.null(factors)) {correlation1 <- correlation1[rownames(ev_cb1), rownames(ev_cb1)]} else {
            if (nrow(ev_cb1) == 0) {
              correlation2 <- correlation2[rownames(ev_cb2), rownames(ev_cb2), drop = FALSE]
            } else {
              correlation1 <- correlation1[rownames(ev_cb1), rownames(ev_cb1), drop = FALSE]
              correlation2 <- correlation2[rownames(ev_cb2), rownames(ev_cb2), drop = FALSE]
            }
            }

          #连续变量的选择n=1
          if (is.null(factors)) {
            if (nrow(ev_cb1) == 0) {bio_name1 <- NULL} else {
              bio_name1 <- corse_method(correlation = correlation1, importance = ev_cb1, vif = vif, n = n)
            }
            bio_name2 <- NULL
          } else{
            if (length(factors) == 1) {
              if (nrow(ev_cb1) == 0) {bio_name1 <- NULL} else {
                bio_name1 <- corse_method(correlation = correlation1, importance = ev_cb1, vif = vif, n = n)
              }
              bio_name2 <- factors} else {
                if (nrow(ev_cb1) == 0) {bio_name1 <- NULL} else {
                  bio_name1 <- corse_method(correlation = correlation1, importance = ev_cb1, vif = vif, n = n)
                }  #下次模拟的变量
                #分类变量的选择
                bio_name2 <- corse_method(correlation = correlation2, importance = ev_cb2, vif=vif, n=n) }} #下次模拟的变量

          #将factors参数更新为bio_name2
          factors <- bio_name2
          #组合两种变量
          bio_name <- c(bio_name2, bio_name1)

          #获取保留的变量存储路径的下标
          evlist <- c()
          for (i in seq_along(bio_name)) {
            evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i]) == TRUE)
            evlist <- c(evlist, evlist1)
          }

          #计算两种变量的相关性
          if (is.null(factors)) {correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
          correlation2 <- as.data.frame(0)
          } else {
            if (is.null(bio_name1)) {correlation1 <- as.data.frame(0)} else {
              correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
            }
            correlation2 <- abs(as.data.frame(cor(mybg2[bio_name2], method = cormethod)))
          }
         # utils::write.csv(correlation1, paste0(outdir, "/TBlabENMtemp", star_time, "/", fc1, rm1, n, "/correlation1.csv"))
         # utils::write.csv(correlation2, paste0(outdir, "/TBlabENMtemp", star_time, "/", fc1, rm1, n, "/correlation2.csv"))
        }

        return(bio_name)}
    } else {


      while(sum(correlation1 > r) > nrow(correlation1)| n==0){ #bio_name初始值为最开始使用的全部变量
        n=n+1

        #进行模拟获得重要性列表
        compare <-suppressMessages(maxent_single(
          x = x,
          evdir = evdir,
          evlist = evlist,
          factors = factors,
          nbg = nbg,
          mybgfile = mybgfile,
          args = args1,
          prodir = NULL,
          outdir = paste0(outdir, "/TBlabENMtemp",star_time,"/",fc1,rm1,n),
          parallel = FALSE))

        #下面根据上面模拟的结果删除相关性强的变量
        #变量重要性
        ev_cb <- read.csv(paste0(outdir,"/TBlabENMtemp",star_time,"/",fc1,rm1,n,"/TBlabENM/maxent/", sp_name, "/maxentResults.csv")) %>%
          dplyr::select(., paste0(bio_name, ".contribution")) %>%
          utils::tail(., n =1) %>% t() %>% as.data.frame()
        names(ev_cb) <- "value"

        ##排除贡献小于0.5的变量
        ev_cb <- ev_cb[ev_cb > 0.5,, drop = FALSE]

        #修改ev_cb的变量名
        nn <- c()
        for (i in 1:nrow(ev_cb)) {
          nn1 <- str_split_1(rownames(ev_cb)[i], ".contribution")[1]
          nn <- c(nn, nn1)
        }
        rownames(ev_cb) <- nn
        #correlation1按照ev_cb提取
        correlation1 <- correlation1[rownames(ev_cb),rownames(ev_cb)]
        bio_name <- corse_method(correlation = correlation1, importance = ev_cb, vif=vif, n=n)  #下次模拟的变量
        evlist <- c()
        for (i in seq_along(bio_name)) {
          evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i])==T)
          evlist <- c(evlist, evlist1) }
        correlation1 <- abs(as.data.frame(cor(mybg[,bio_name], method = cormethod)))

      }
      return(bio_name)
    }
  }

  fun4 <- function(i){
    #设置好fc和rm后进行模拟保留最重要变量i=1
    var <- fun2(fc1 = combin[i,1] ,rm1 = combin[i,2])
    df[i,3] <- paste0(var, collapse = ",")
  }


  #####################################################################
  if (is.null(myenv) == FALSE) {
    cat("The environment variable has been specified and is no longer selected !")
    evlist <- NULL}
  if (is.null(opt) == FALSE) {
    if (is.null(myenv) == FALSE & length(myenv) < 4) {stop("The number of environment variables must be at least three")}
    if (is.null(evlist) == FALSE & length(evlist) < 4) {stop("The number of environment variables must be at least three")}}

  tzhs <- c("L", "Q", "H", "P", "T")
  if(is.null(outdir)){outdir <- "."}
  #dir.create(paste0(outdir, "/TBlabENM"),recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "/TBlabENM/data"), showWarnings = FALSE, recursive = TRUE)
  factors123 <- factors
  biolistall <- list.files(evdir, pattern = ".asc$", full.names = TRUE)

  #获取物种名 对路径拆分并取倒数第一个字符串
  spname1 <- stringr::str_split_1(x, "/")[length(stringr::str_split_1(x, "/"))]
  sp_name <- stringr::str_split_1(spname1, ".csv$")[1]
  cat(paste("Select optimum parameters for", sp_name, "\n"))
  #################################################################################
  #当指定环境变量时不再进行变量的选择
  if(is.null(myenv)){
  #变量名称
  biolist <- list.files(evdir, pattern = ".asc$", full.names = TRUE)
  if(is.null(evlist) == FALSE){biolist <- biolist[evlist]}
  bio_name <- c()
  for(i in seq_along(biolist)){
    bioname1 <- stringr::str_split_1(biolist[i], "/")[length(stringr::str_split_1(biolist[i], "/"))]
    bio_name0 <- stringr::str_split_1(bioname1, ".asc")[1]
    bio_name <- c(bio_name, bio_name0)
  }
  cat(paste("Initial variable:", paste(bio_name, collapse = ","), "\n"))

  #判断分类变量factors是否包含在给定的环境数据集内
  if(is.null(factors)){cat("All variables are continuous variables \n")}else{
    factors1 <- factors
    factors <- factors[factors %in% bio_name]
    if(length(factors)==0){
      stop("No categorical variable is found in the given list of variables, and parameter 'factors' is invalid. Please check whether parameter 'evdir','evlist','factors' are set correctly.")
    }
    if(length(factors)!= length(factors1)){
      stop(paste0("Only variable '", factors,  "' is identified as a categorical variable. Please check whether other categorical variables are in the given set of environment variables"))}
  cat(paste("Variable:", paste(factors, collapse = ","), "as categorical variable", "\n"))}
  cat("***********Selecting parameters***********\n")
  #提取发生数据的环境值
  biostack <- terra::rast(biolist)

  occ <- utils::read.csv(x) #读取物种坐标数据

  occ <- occ[c(2,3)]

  occdata <- terra::extract(biostack, occ, ID=FALSE)

  #提取背景值并计算变量相关性
  ##随机生成10000个点
  if(is.null(mybgfile)) {
    mybg0 <- terra::spatSample(biostack, nbg, na.rm = T, xy = T)
    mybgfile <- mybg0[1:2]
    write.csv(mybgfile, paste0(outdir, "/TBlabENM/data/", sp_name, "_bg.csv"), row.names = FALSE)
    mybg <- mybg0[-(1:2)]} else{
    cat("The background points has been specified !")
    mybg <- terra::extract(biostack, mybgfile, ID = FALSE)
  }
  ##判断是否存在分类变量，若存在则分开计算
  if(is.null(factors) == FALSE){#mybg分为两组，mybg1为连续变量，mybg2为分类变量
    mybg1 <- mybg[,bio_name[!bio_name %in% factors], drop = FALSE]
    mybg2 <- mybg[,factors, drop = FALSE]
    correlation1 <- abs(as.data.frame(cor(mybg1, method = cormethod)))
    correlation2 <- abs(as.data.frame(cor(mybg2, method = cormethod)))


  } else { correlation1 <- abs(as.data.frame(cor(mybg, method = cormethod)))

  }

  #组合fc和 rm
  fc <- toupper(fc)
  combin <- expand.grid(fc, rm, stringsAsFactors = FALSE) %>%
    dplyr::mutate(fc = purrr::map_chr(.x = Var1, .f = function(x){
      paste0(tzhs[tzhs %in% stringr::str_split_1(x, "")], collapse = "")}
    ))
  combin <- combin[c(3,2)]
  df <- data.frame(matrix(NA,nrow(combin),1))
  df <- cbind(combin, df)
  names(df) <- c("fc", "rm", "env")

  #根据组合设置 args参数
  args1 <- TBlabENM::maxent_args(l=FALSE , q=FALSE, p=FALSE, h=FALSE, t=FALSE,
                       responsecurves=FALSE, jackknife=FALSE, pictures=FALSE)
##并行计算
  if(parallel == T){
    # 开启集成
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    # 注册每个环境变量
    snowfall::sfExport("fun2")
    snowfall::sfExport("star_time")
    snowfall::sfLibrary(purrr)
    snowfall::sfLibrary(TBlabENM)
    snowfall::sfExport("x")
    snowfall::sfExport("evdir")
    snowfall::sfExport("evlist")
    snowfall::sfExport("mybgfile")
    snowfall::sfExport("nbg")
    snowfall::sfExport("combin")
    snowfall::sfExport("args1")
    snowfall::sfExport("factors")
    snowfall::sfExport("correlation1")
    snowfall::sfExport("r")
   # snowfall::sfExport("correlation2")
   # snowfall::sfExport("maxent_single")
    k <- snowfall::sfLapply(1:nrow(combin), fun4)
    snowfall::sfStop()  # 关闭集群

    df$env <- unlist(k)


  }else{
  for(i in 1:nrow(combin)){


      #设置好fc和rm后进行模拟保留最重要变量i=1
      var <- fun2(fc1 = combin[i,1] ,rm1 = combin[i,2])
      df[i,3] <- paste0(var, collapse = ",")
  } }
  } else {
    biolist <- list.files(evdir, pattern = ".asc$", full.names = TRUE)
    #判断分类变量factors是否包含在给定的环境数据集内
    ##全部变量名称
    bio_name_all <- c()
    for(i in seq_along(biolistall)){
      bioname1 <- stringr::str_split_1(biolistall[i], "/")[length(stringr::str_split_1(biolistall[i], "/"))]
      bio_name0 <- stringr::str_split_1(bioname1, ".asc")[1]
      bio_name_all <- c(bio_name_all, bio_name0)
    }

    #指定的变量
    bio_name <- myenv

    if(sum(bio_name %in% bio_name_all)!=length(bio_name)){
      stop("One or more specified variables do not exist! Pleaes check parameter 'evdir' and 'myenv'.")
    } else { cat(paste("Initial variable:", paste(bio_name, collapse = ","), "\n"))}
    if(is.null(factors)){cat("All variables are continuous variables \n")}else{
      factors1 <- factors
      factors <- factors[factors %in% bio_name]
      if(length(factors)==0){
      stop("No categorical variable is found in the given list of variables, and parameter 'factors' is invalid. Please check whether parameter 'evdir','evlist','factors' are set correctly.")
      }
      if(length(factors)!=length(factors1)){
        stop(paste0("Only variable '", factors,  "' is identified as a categorical variable. Please check whether other categorical variables are in the given set of environment variables"))}
      cat(paste("Variable:", paste(factors, collapse = ","), "as categorical variable", "\n"))}

#提取发生数据的环境值
##获取bio_name的下标
xb <- c()
for (i in bio_name) {
  xb1 <- which( stringr::str_detect(biolistall, paste0(i, ".asc"))==T)
  xb <- c(xb, xb1) }
biostack <- terra::rast(biolistall[xb])
occ <- utils::read.csv(x) #读取物种坐标数据
occ <- occ[c(2,3)]
occdata <- terra::extract(biostack, occ, ID=FALSE)

#提取背景值并计算变量相关性
##随机生成10000个点
if(is.null(mybgfile)){
  mybg <- terra::spatSample(biostack, nbg, na.rm = T, xy = T)[-(1:2)]} else{
    mybg <- terra::extract(biostack, mybgfile, ID=FALSE)
  }

    #组合fc和 rm
    fc <- toupper(fc)
    combin <- expand.grid(fc, rm, stringsAsFactors = FALSE) %>%
      mutate(fc = purrr::map_chr(.x = Var1, .f = function(x){
        paste0(tzhs[tzhs %in% stringr::str_split_1(x, "")], collapse = "")}
      ))
    combin <- combin[c(3,2)]
    df <- data.frame(matrix(NA,nrow(combin),1))
    df <- cbind(combin, df)
    df[,3] <- rep(paste0(myenv, collapse = ","), nrow(df))
    names(df) <- c("fc", "rm", "env")

  }

if(is.null(opt)){parameter <- df} else{
#剩余的组合进行evaluate通过AICC确定最佳组合
  df <- df %>%
    dplyr::mutate(occdata = purrr::map(.x = env, .f = function(x){
      #获取变量下标
      xb <- c()
      xb2 <- stringr::str_split_1(x, ",")
      for (i in seq_along(xb2)) {
        xb1 <- which( stringr::str_detect(biolist, paste0(xb2, ".asc")[i])==T)
        xb <- c(xb, xb1) }

      occdata <- occdata[xb2]
    }
    )) %>%
    dplyr::mutate(bgdata = purrr::map(.x = env,.f = function(x){
      #获取变量下标
      xb <- c()
      xb2 <- stringr::str_split_1(x, ",")
      for (i in seq_along(xb2)) {
        xb1 <- which( stringr::str_detect(biolist, paste0(xb2, ".asc")[i])==T)
        xb <- c(xb, xb1) }
      bgdata <- mybg[xb2]
    }
    )) %>% #计算变量个数
    mutate(num = purrr::map_dbl(.x = occdata, .f = function(x){
      ncol(x)
    }))
  ##当保留的变量数只有三个时ENMevaluate会报错，需要移除只保留三个变量的结果
  df1 <- dplyr::filter(df, num<4)
  df <- dplyr::filter(df, num>=4)
  if(nrow(df)==0){stop("All combined environment variables have fewer than four")}

#使用pmap函数并行评估（由于ENMevaluate函数的参数大于2个，所以使用pmap函数）
  #设置参数
  if(nrow(occdata)>=25){
    partitions = "randomkfold"
    partition.settings = list(kfolds = 10)} else {
      partitions = "jackknife"
      partition.settings = NULL}

  gz <- purrr::pmap(df, .f = function(fc,rm,occdata, bgdata, ...){
    e <- ENMeval::ENMevaluate(
      occs = occdata,
      bg = bgdata,
      tune.args = list(fc = fc, rm = rm),
      #partitions = "jackknife", #数据分区方式，有2+6种
      partitions = partitions, #数据分区方法|测试集、训练集划分方法
      partition.settings =  partition.settings, #数据分区设置
      other.settings =  list(abs.auc.diff = TRUE,
                             pred.type = "logistic",
                             validation.bg = "full", #与partitions参数分区类型对应
                             other.args = NULL),#其他额外设置，有默认值
      taxon.name = sp_name,
      n.bg = nbg,
      algorithm = "maxent.jar", #使用的模型，有三种
      overlap = FALSE, #生态位重叠
      categoricals = factors123[factors123%in%names(occdata)], #指定分类变量,"IAWC_CLASS", , "T_USDA_TEX_CLASS"
      doClamp = FALSE)
      res <- ENMeval::eval.results(e)
  })

  #评估的一系列参数
  gz1 <- as.data.frame(data.table::rbindlist(gz))

  #合并df和gz1

  cs <- cbind(df[1:3], gz1[c(1:16,19)])
  cs <- cs[-(4:5)]
  cs <- dplyr::arrange(cs, AICc, auc.diff.avg, or.mtp.avg)

  #依据aicc选择最佳模型
  if(opt=="aicc") {opt1 <- dplyr::filter(cs, AICc == min(AICc))
  if(nrow(opt1)==0){opt1 <- dplyr::filter(cs, auc.val.avg == max(auc.val.avg))
  warning("AICC is NA, use 'auc.val.avg' instead.")}

    }
  #使用顺序法选择最佳模型
  if(opt=="seq") {
    opt1 <- cs %>%
      dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
        dplyr::filter(auc.val.avg == max(auc.val.avg))

    }
  #使用cbi选择最佳模型
  if(opt=="cbi") {opt1 <- dplyr::filter(cs, cbi.val.avg == max(cbi.val.avg))
  if(nrow(opt1)==0){opt1 <- dplyr::filter(cs, auc.val.avg == max(auc.val.avg))
  warning("cbi.val.avg is NA, use 'auc.val.avg' instead.")
  }
  }
  #保存结果
  #删除缓存文件
  unlink(paste0(outdir, "/TBlabENMtemp",star_time) ,recursive = T)
  #dir.create(paste0(outdir, "/TBlabENM"), showWarnings = FALSE)

print(paste0(outdir, "/TBlabENM/tuneparameter_", sp_name, ".csv"))

  utils::write.csv(cs, paste0(outdir, "/TBlabENM/data/tuneparameter_", sp_name, ".csv"), row.names = FALSE)
  parameter <- opt1[1,1:3]
  parameter$species <- sp_name
  parameter$number <- nrow(occ)
  parameter <- parameter[,c(4,5,3,1,2)]
  if(nrow(df1)>=1 ){warning("The following combinations were not evaluated because the number of environment variables was less than four")
    print(df1[1:3])}}
  return(parameter)
    }



