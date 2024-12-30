
#' Select optimum parameters for the MaxEnt models
#' @description
#' 本函数为MaxEnt模型选择最适合的环境变量组合、正则化乘数和特征函数组合。
#' @details
#' 工作流程：首先对fc和rm进行所有可能的组合，对于每个组合，使用全部变量进行建模获得变量贡献大小，依据变量贡献大小排除与最重要变量具有强相关性的变量。然后使用保留的变量继续建模获得该组变量贡献大小，依据变量贡献大小排除与第二重要的变量具有强相关性的变量。重复上述过程直到保留的变量之间不具有强相关性。经过此步骤就为每个fc和rm组合获得了一组不具有强相关性且保留了贡献大的候选变量，当vif=T时，还对这一系列候选变量组合进行方差膨胀因子分析，进一步排除具有强共线性的组合。最后使用ENMevaluate对每个组合进行测试，获得相应的指标（例如“AICc”），根据指标选择最佳的一组组合作为模型最终参数。
#' 注意：当物种发生数据过少时，可能无法得出正确的AICc和cbi值，如果opt选择AICc或cbi，则使用auc.val.avg值进行代替。
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
#' @param opt 选择最佳模型的指标，可选择"auc.train, cbi.train, auc.diff.avg, auc.val.avg, cbi.val.avg, or.10p.avg, or.mtp.avg, AICc"或NULL.当选择NULL时，则仅选择变量。当选择多个指标时，则按照顺序筛选最佳模型，如果所选指标具有NA，则跳过该指标，如果所有指标都含有NA，则使用auc.val.avg代替.
#' @param outdir 结果保存路径
#' @param null_model 逻辑值,是否进行NULL模型检验。
#' @param nbg 随机背景点数量，当指定参数mybgfile时忽略。
#' @param myenv 指定的一组环境变量，当提供此参数时，将不会再根据相关性进行变量筛选
#' @param parallel 是否并行计算
#' @param ncpu 如果并行，使用的cpu数
#' @param bgwidth 数值型(单位m), 以发生点为中心, 以bgwidth为半径创建一个缓冲区, 用来选择背景点.如果为NULL则在整个环境层范围内选择背景点.
#'
#' @importFrom dplyr mutate arrange select
#' @importFrom stringr str_split_1
#' @importFrom utils read.csv tail write.csv
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr map map_dbl
#' @importFrom foreach %dopar%
#' @importFrom ggplot2 ggsave
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
#'        bgwidth = 500000,
#'        fc = c("LPH", "Q", "LQ"),
#'        rm = 1:5,
#'        r = 0.7,
#'        vif = T,
#'        vifth = 5,
#'        opt = "auc.val.avg",
#'        parallel = F,
#'        ncpu = 2)

maxent_parameter <- function(x,
                             evdir,
                             myenv = NULL,
                             evlist = NULL,
                             factors = NULL,
                             mybgfile = NULL,
                             nbg = 10000,
                             bgwidth = NULL,
                             fc,
                             rm,
                             r = 0.7,
                             cormethod = "pearson",
                             vif = T,
                             vifth = 5,
                             opt = "auc.val.avg",
                             null_model = TRUE,
                             outdir = NULL,
                             parallel = F,
                             ncpu = 2) {

  #判断参数格式是否正确
  if (is.null(opt) == FALSE) {
    for (i in 1:length(opt)) {
      if (!opt[i] %in% c("auc.train", "cbi.train", "auc.diff.avg", "auc.val.avg", "cbi.val.avg", "or.10p.avg", "or.mtp.avg", "AICc")) {
        stop("'opt' must be a subset of c('auc.train', 'cbi.train', 'auc.diff.avg', 'auc.val.avg', 'cbi.val.avg', 'or.10p.avg', 'or.mtp.avg', 'AICc') or NULL.")
      }
    }
    }
  random_num <- sample(1:100000, 1)
  #绘制模型调优结果的函数
  p_tun <- function(cs1, opt) {
    data <-
      pivot_longer(
        cs1,
        cols = c("auc.train", "cbi.train", "auc.diff.avg", "auc.val.avg",
                 "cbi.val.avg", "or.10p.avg", "or.mtp.avg", "AICc"),
        names_to = "name",
        values_to = "value"
      ) %>%
      dplyr::filter(., name %in% opt)

    if (length(opt) > 1) {
      p <- ggplot(data, mapping = aes(
        x = rm,
        y = value,
        color = fc,
        group = fc
      )) +
        geom_point(size = 2) +
        geom_line() +
        facet_grid(name ~ ., scales = "free") +
        xlab("rm") +
        labs(color = "fc") +
        theme_bw()

      ggsave(
        filename = "model_tun.jpg",
        plot = p,
        path = paste0(outdir, "/maxent/", sp_name),
        width = 14 + (max(cs1$rm) - 4.5)*3,
        height = 7*length(opt),
        units = c("cm"),
        create.dir = TRUE,
        #limitsize = FALSE
      )
    } else {
      p <- ggplot(data, mapping = aes(
        x = rm,
        y = value,
        color = fc,
        group = fc
      )) +
        geom_point(size = 2) +
        geom_line() +
        xlab("rm") +
        ylab(opt) +
        labs(color = "fc") +
        theme_bw()

      ggsave(
        filename = "model_tun.jpg",
        plot = p,
        path = paste0(outdir, "/maxent/", sp_name),
        width = 14,
        height = 10,
        units = c("cm"),
        create.dir = TRUE,
        #limitsize = FALSE
      )
    }
  }

  #修改ENMnulls函数
  ENMnulls <- function(e,
                       mod.settings,
                       no.iter,
                       eval.stats = c("auc.val", "auc.diff", "cbi.val", "or.mtp", "or.10p"),
                       user.enm = NULL,
                       user.eval.type = NULL,
                       userStats.signs = NULL,
                       removeMxTemp = TRUE,
                       parallel = FALSE,
                       numCores = NULL,
                       parallelType = "doSNOW",
                       quiet = FALSE) {
    # assign evaluation type based on partition method
    if (is.null(user.eval.type)) {
      eval.type <- switch(
        e@partition.method,
        randomkfold = "knonspatial",
        jackknife = "knonspatial",
        block = "kspatial",
        checkerboard1 = "kspatial",
        checkerboard2 = "kspatial",
        testing = "testing",
        none = "none"
      )
    } else{
      eval.type <- user.eval.type
    }

    ## checks
    # model settings are all single entries
    if (!all(sapply(mod.settings, length) == 1))
      stop("Please input a single set of model settings.")
    # model settings are correct for input algorithm and are entered in the right order --
    #if not, put them in the right order, else indexing models later will fail because the model
    # name will be incorrect
    if (e@algorithm %in% c("maxent.jar", "maxnet")) {
      if (length(mod.settings) != 2) {
        stop(
          "Please input two complexity settings (fc [feature classes] and rm [regularization
           multipliers]) for mod.settings for maxent.jar and maxnet models."
        )
      }
      if (all(names(mod.settings) %in% c("fc", "rm"))) {
        if (!all(names(mod.settings) == c("fc", "rm"))) {
          mod.settings <- mod.settings[c("fc", "rm")]
        }
      } else{
        stop(
          'Please input only "fc" (feature classes) and "rm" (regularization multipliers) for
           mod.settings for maxent.jar and maxnet models.'
        )
      }
    } else if (e@algorithm == "bioclim") {
      if (length(mod.settings) != 1) {
        stop(
          "Please input one complexity setting (tails) for mod.settings for BIOCLIM models."
        )
      }
      if (!all(names(mod.settings) == "tails")) {
        stop('Please input only "tails" for mod.settings for BIOCLIM models.')
      }
    }


    # assign directionality of sign for evaluation stats
    signs <- c(
      list(
        "auc.val" = 1,
        "auc.train" = 1,
        "cbi.val" = 1,
        "cbi.train" = 1,
        "auc.diff" = -1,
        "or.10p" = -1,
        "or.mtp" = -1
      ),
      userStats.signs
    )

    # record start time
    start.time <- proc.time()

    # assign the number of cross validation iterations
    nk <- max(as.numeric(as.character(e@occs.grp)))

    # get number of occurrence points by partition
    occs.grp.tbl <- table(e@occs.grp)

    # if more than one background partition exists, assume spatial CV and
    # keep existing partitions
    null.samps <- cbind(rbind(e@occs, e@bg), grp = c(e@occs.grp, e@bg.grp))

    if (e@algorithm == "maxent.jar") {
      # create temp directory to store maxent.jar output, for potential removal later
      tmpdir <- paste(tempdir(), runif(1, 0, 1), sep = "/")
      dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
    }

    # assign user algorithm if provided
    if (!is.null(user.enm)) {
      e@algorithm <- user.enm
    }

    ############################## #
    # specify empirical model statistics ####
    ############################## #

    mod.tune.args  <- paste(names(mod.settings),
                            mod.settings,
                            collapse = "_",
                            sep = ".")
    emp.mod <- e@models[[mod.tune.args]]
    emp.mod.res <- e@results %>% dplyr::filter(tune.args == mod.tune.args)

    ############################## #
    # build null models ####
    ############################## #

    if (quiet == FALSE)
      message(paste(
        "Building and evaluating null ENMs with",
        no.iter,
        "iterations..."
      ))
    if (quiet == FALSE)
      pb <- txtProgressBar(0, no.iter, style = 3)

    # set up parallel processing functionality
    if (parallel == TRUE) {
      allCores <- parallel::detectCores()
      if (is.null(numCores)) {
        numCores <- allCores
      }
      cl <- parallel::makeCluster(numCores, setup_strategy = "sequential")
      if (quiet != TRUE)
        progress <- function(n)
          setTxtProgressBar(pb, n)

      if (parallelType == "doParallel") {
        doParallel::registerDoParallel(cl)
        opts <- NULL
      } else if (parallelType == "doSNOW") {
        doSNOW::registerDoSNOW(cl)
        if (quiet != TRUE)
          opts <- list(progress = progress)
        else
          opts <- NULL
      }
      numCoresUsed <- foreach::getDoParWorkers()
      if (quiet != TRUE)
        message(paste0("\nOf ", allCores, " total cores using ", numCoresUsed, "..."))
      if (quiet != TRUE)
        message(paste0("Running in parallel using ", parallelType, "..."))
    }

    if (length(e@clamp.directions) == 0)
      clamp.directions.i <- NULL
    else
      clamp.directions.i <- e@clamp.directions

    # define function to run null model for iteration i
    null_i <- function(i) {
      null.occs.ik <- list()
      if (eval.type == "kspatial") {
        # randomly sample the same number of training occs over each k partition
        # of envs; if kspatial evaluation, only sample over the current spatial
        # partition of envs.z
        for (k in 1:nk) {
          # sample null occurrences only from
          # the records in partition group k
          null.samps.k <- null.samps %>% dplyr::filter(grp == k)
          # randomly sample n null occurrences, where n equals the number
          # of empirical occurrence in partition group k
          samp.k <- sample(1:nrow(null.samps.k), occs.grp.tbl[k])
          null.occs.ik[[k]] <- null.samps.k[samp.k, ]
        }
      } else if (eval.type == "knonspatial") {
        for (k in 1:nk) {
          # randomly sample n null occurrences, where n equals the number
          # of empirical occurrence in partition group k
          samp.k <- sample(1:nrow(null.samps), occs.grp.tbl[k])
          null.occs.ik[[k]] <- null.samps[samp.k, ]
        }
      } else if (eval.type %in% c("testing", "none")) {
        samp.test <- sample(1:nrow(null.samps), occs.grp.tbl)
        null.occs.ik[[1]] <- null.samps[samp.test, ]
      }

      # bind rows together to make full null occurrence dataset
      null.occs.i.df <- dplyr::bind_rows(null.occs.ik)
      if (eval.type == "knonspatial") {
        if (e@partition.method == "randomkfold")
          null.occs.i.df$grp <- get.randomkfold(null.occs.i.df,
                                                e@bg,
                                                kfolds = e@partition.settings$kfolds)$occs.grp
        if (e@partition.method == "jackknife")
          null.occs.i.df$grp <- get.jackknife(null.occs.i.df, e@bg)$occs.grp
      }
      null.occs.i.z <- null.occs.i.df %>% dplyr::select(-grp)
      # shortcuts
      categoricals <- names(which(sapply(e@occs, is.factor)))
      if (length(categoricals) == 0)
        categoricals <- NULL

      if (eval.type %in% c("testing", "none")) {
        partitions <- eval.type
        user.grp <- NULL
        user.val.grps <- NULL
      } else{
        # assign the null occurrence partitions as user partition settings, but
        # keep the empirical model background partitions
        user.grp <- list(occs.grp = null.occs.i.df$grp, bg.grp = e@bg.grp)
        # assign user validation partitions to those used in the empirical model
        user.val.grps <- cbind(e@occs, grp = e@occs.grp)
        partitions <- "user"
      }

      # check if ecospat is installed, and if not, prevent CBI calculations
      if ((requireNamespace("ecospat", quietly = TRUE))) {
        e@other.settings$ecospat.use <- TRUE
      } else{
        e@other.settings$ecospat.use <- FALSE
      }

      args.i <- list(
        occs = null.occs.i.z,
        bg = e@bg,
        tune.args = mod.settings,
        categoricals = categoricals,
        partitions = partitions,
        algorithm = e@algorithm,
        other.settings = e@other.settings,
        partition.settings = e@partition.settings,
        occs.testing = e@occs.testing,
        user.val.grps = user.val.grps,
        user.grp = user.grp,
        doClamp = e@doClamp,
        clamp.directions = clamp.directions.i,
        quiet = TRUE
      )

      null.e.i <- tryCatch({
        do.call(ENMevaluate, args.i)
      }, error = function(cond) {
        if (quiet != TRUE)
          message(paste0("\n", cond, "\n"))
        # Choose a return value in case of error
        return(NULL)
      })

      if (is.null(null.e.i)) {
        results.na <- e@results[1, ] %>% dplyr::mutate(dplyr::across(auc.train:ncoef, ~
                                                                       NA))
        mod.settings.i <- paste(names(mod.settings),
                                mod.settings,
                                collapse = "_",
                                sep = ".")
        if (nrow(e@results.partitions) > 0) {
          results.partitions.na <- e@results.partitions %>% dplyr::filter(tune.args == mod.settings.i) %>% dplyr::mutate(dplyr::across(3:ncol(.), ~
                                                                                                                                         NA)) %>% dplyr::mutate(iter = i)
        } else{
          results.partitions.na <- e@results.partitions
        }

        out <- list(results = results.na, results.partitions = results.partitions.na)
      } else{
        out <- list(
          results = null.e.i@results,
          results.partitions = null.e.i@results.partitions %>% dplyr::mutate(iter = i) %>% dplyr::select(iter, dplyr::everything())
        )
        # restore NA row if partition evaluation is missing (model was NULL)
        if (eval.type != "testing") {
          allParts <- unique(user.grp$occs.grp) %in% out$results.partitions$fold
          if (!all(allParts)) {
            inds <- which(allParts == FALSE)
            newrow <- out$results.partitions[1, ]
            newrow[, 4:ncol(newrow)] <- NA
            for (ind in inds) {
              out$results.partitions <- dplyr::bind_rows(out$results.partitions,
                                                         newrow %>% dplyr::mutate(fold = ind))
            }
            out$results.partitions <- dplyr::arrange(out$results.partitions, fold)
          }
        }
      }

      return(out)
    }

    if (parallel == TRUE) {
      outs <- foreach::foreach(
        i = 1:no.iter,
        .export = 'e',
        .options.snow = opts,
        .packages = c("dplyr", "ENMeval")
      ) %dopar% {
        null_i(i)
      }
    } else{
      outs <- list()
      for (i in 1:no.iter) {
        outs[[i]] <- null_i(i)
        if (quiet == FALSE)
          setTxtProgressBar(pb, i)
      }
    }

    if (quiet != TRUE)
      close(pb)
    if (parallel == TRUE)
      parallel::stopCluster(cl)

    # assemble null evaluation statistics and take summaries
    nulls.ls <- lapply(outs, function(x)
      x$results)
    nulls.grp.ls <- lapply(outs, function(x)
      x$results.partitions)
    nulls <- dplyr::bind_rows(nulls.ls) %>% dplyr::select(!dplyr::contains("AIC"))
    nulls.grp <- dplyr::bind_rows(nulls.grp.ls)
    if (eval.type %in% c("testing", "none")) {
      nulls.avgs <- nulls %>% dplyr::select(dplyr::ends_with("train"),
                                            dplyr::contains(eval.stats)) %>% dplyr::summarize_all(mean, na.rm = TRUE)
      nulls.sds <- nulls %>% dplyr::select(dplyr::ends_with("train"),
                                           dplyr::contains(eval.stats)) %>% dplyr::summarise_all(sd, na.rm = TRUE)
      # get empirical model evaluation statistics for comparison
      emp.avgs <- emp.mod.res %>% dplyr::select(dplyr::ends_with("train"),
                                                dplyr::contains(eval.stats))
    } else{
      nulls.avgs <- nulls %>% dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".avg")) %>% dplyr::summarize_all(mean, na.rm = TRUE)
      nulls.sds <- nulls %>% dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".sd")) %>% dplyr::summarise_all(sd, na.rm = TRUE)
      # get empirical model evaluation statistics for comparison
      emp.avgs <- emp.mod.res %>% dplyr::select(dplyr::ends_with("train"), paste0(eval.stats, ".avg"))
    }
    if (sum(grepl("sd", names(emp.mod.res))) > 0) {
      emp.sds <- emp.mod.res %>% dplyr::select(paste0(eval.stats, ".sd"))
    } else{
      emp.sds <- NULL
    }

    empNull.stats <- as.data.frame(matrix(nrow = 6, ncol = ncol(emp.avgs) +
                                            1))
    names(empNull.stats)[1] <- "statistic"
    empNull.stats$statistic <- c("emp.mean",
                                 "emp.sd",
                                 "null.mean",
                                 "null.sd",
                                 "zscore",
                                 "pvalue")
    names(empNull.stats)[-1] <- gsub(".avg", replacement = "", names(emp.avgs))

    # populate empirical and null means and standard deviations
    empNull.stats[1, -1] <- emp.avgs
    emp.sds.nameStrip <- gsub(".sd", "", names(emp.sds))
    if (length(emp.sds.nameStrip) > 0)
      empNull.stats[2, which(names(empNull.stats) %in% emp.sds.nameStrip)] <- emp.sds
    empNull.stats[3, -1] <- nulls.avgs
    empNull.stats[4, -1] <- nulls.sds
    # calculate z-scores
    empNull.stats[5, -1] <- (emp.avgs - nulls.avgs) / nulls.sds
    # find statistics that need a positive pnorm, and those that need a negative pnorm
    p.pos <- names(signs[sapply(signs, function(x)
      x == 1)])
    p.neg <- names(signs[sapply(signs, function(x)
      x == -1)])
    p.pos.inds <- which(grepl(paste(p.pos, collapse = "|"), names(empNull.stats)))
    p.neg.inds <- which(grepl(paste(p.neg, collapse = "|"), names(empNull.stats)))
    # calculate p-values
    empNull.stats[6, p.pos.inds] <- sapply(empNull.stats[5, p.pos.inds], function(x)
      pnorm(x, lower.tail = FALSE))
    empNull.stats[6, p.neg.inds] <- sapply(empNull.stats[5, p.neg.inds], pnorm)

    mod.settings.tbl <- expand.grid(mod.settings)
    mod.settings.tbl$tune.args <- apply(mod.settings.tbl, 1, function(x)
      paste(names(x), x, collapse = "_", sep = "."))
    mod.settings.tbl <- dplyr::mutate_all(mod.settings.tbl, as.factor)

    # make cbi.val column NA for jackknife partitions
    if (e@partition.method == "jackknife") {
      empNull.stats$cbi.val <- NA
      nulls$cbi.val.avg <- NA
      nulls$cbi.val.sd <- NA
      nulls.grp$cbi.val <- NA
    }

    # first, add clamp directions to other settings so the object can record it
    e@other.settings$clamp.directions <- e@clamp.directions
    # condense mod.args to named matrix for inserting into class slot
    e.n <- ENMnull(
      null.algorithm = e@algorithm,
      null.mod.settings = mod.settings.tbl,
      null.partition.method = e@partition.method,
      null.partition.settings = e@partition.settings,
      null.doClamp = e@doClamp,
      null.other.settings = e@other.settings,
      null.no.iter = no.iter,
      null.results = nulls,
      null.results.partitions = nulls.grp,
      null.emp.results = empNull.stats,
      emp.occs = e@occs,
      emp.occs.grp = e@occs.grp,
      emp.bg = e@bg,
      emp.bg.grp = e@bg.grp
    )

    # optionally remove temp directory for maxent.jar
    if (e@algorithm == "maxent.jar" &
        removeMxTemp == TRUE)
      unlink(tmpdir, recursive = TRUE, force = TRUE)


    timed <- proc.time() - start.time
    t.min <- floor(timed[3] / 60)
    t.sec <- timed[3] - (t.min * 60)
    if (quiet == FALSE)
      message(paste(
        "\nENMnulls completed in",
        t.min,
        "minutes",
        round(t.sec, 1),
        "seconds."
      ))

    return(e.n)
  }

  #corse_method功能使用相关性选择变量
  corse_method <- function(correlation, importance, vif, n) {
    #n最初为1

    n = n
    #这里可能会出现变量重要性不一致问题，之前选出来的变量直接保留占位，其他变量排序
    if (n == 1) {
      importance1 <- dplyr::arrange(importance, desc(value))
    } else {
      importance_a <- importance[1:(n - 1), , drop = FALSE]
      importance_b <- dplyr::arrange(importance[-(1:(n - 1)), , drop = FALSE], desc(value))
      importance1 <- rbind(importance_a, importance_b)
    }
    #importance1 <- dplyr::arrange(importance, desc(value))
    if (nrow(importance) < n) {
      importance <- rownames(importance1)[1]
    } else{
      importance <- rownames(importance1)[n]
    }

    ff <- correlation[, importance, drop = FALSE]
    #ff <- dplyr::filter(ff, rowSums(ff < r) == ncol(ff))
    ff <-  ff[which(ff < r), , drop = FALSE]

    bio_name <- c(importance, rownames(ff)) #要保留的变量名称


    #当执行vif分时
    if (vif == T) {
      bio_name <- vifse_method(
        envdf = mybg,
        env_name = bio_name,
        importance = importance1,
        vifth = vifth,
        n = n
      )
      #importance继承自importance1，importance1是排过序的，因此importance不用再排序
    }

    return(bio_name)
  }
  #vifse_method功能使用vif选择变量
  vifse_method <- function(envdf, env_name, importance, vifth, n) {
    #对保留的变量进行重要性排序
    #im <- dplyr::arrange(importance[env_name,, drop =FALSE],desc(value))
    #env_name <- rownames(im) #排过序的变量名
    #importance中取出保留的变量
    im <- importance[rownames(importance) %in% env_name, , drop = FALSE]
    env_name <- rownames(im) #这里的env_name是排过序的变量名
    env_name1 <- env_name
    n = n - 1

    while (n < length(env_name1) - 1) {
      n = n + 1
      #从env_name中取出前n+1个变量，判断最重要的变量和次重要的变量之间的共线性，如果存在，次重要的变量删除
      var <- env_name1[c(1:(n + 1))]
      #只保留在env_name1中的变量
      var1 <- var[var %in% env_name]
      #对var1计算vif
      vifvalue <- usdm::vif(envdf[var1])

      #当存在共线性则应该排除var1中的最后一个变量,当不存在共线性说明该变量要保留，进入下次循环
      if (sum(vifvalue$VIF > vifth) > 0) {
        env_name <- env_name[!env_name %in% var1[length(var1)]]
      } else{
        break
      }
    }

    return(env_name)
  }

  fun2 <- function(fc1, rm1) {
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
    if (is.null(factors) == FALSE) {
      #当存在分类变量，分开计算
      #当factors只有1个变量
      if (length(factors) == 1) {
        while (sum(correlation1 > r) > nrow(correlation1) | n == 0) {
          n <- n + 1

          #进行模拟获得重要性列表
          compare <- suppressMessages(
            maxent_single(
              x = x,
              evdir = evdir,
              evlist = evlist,
              #evlist应该是相对于evdir也就是全部变量biolstall的下标
              factors = factors,
              nbg = nbg,
              bgwidth = NULL,
              mybgfile = mybgfile,
              args = args1,
              prodir = NULL,
              outdir = paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n),
              parallel = FALSE
            )
          )
          #下面根据上面模拟的结果删除相关性强的变量
          #变量重要性
          ev_cb <- read.csv(
            paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n, "/maxent/", sp_name, "/maxentResults.csv")
          ) %>%
            dplyr::select(., paste0(bio_name, ".contribution")) %>%
            utils::tail(., n = 1) %>% t() %>% as.data.frame()
          names(ev_cb) <- "value"

          ##排除贡献小于0.5的变量
          ev_cb <- ev_cb[ev_cb > 0.5, , drop = FALSE]

          #修改ev_cb的变量名
          nn <- c()
          for (i in 1:nrow(ev_cb)) {
            nn1 <- stringr::str_split_1(rownames(ev_cb)[i], ".contribution")[1]
            nn <- c(nn, nn1)
          }
          rownames(ev_cb) <- nn
          #将ev_cb按照变量类型拆开
          factors <- rownames(ev_cb)[rownames(ev_cb) %in% factors] #更新factors
          if (length(factors) == 0) {
            factors <- NULL
          }

          if (is.null(factors)) {
            ev_cb1 <- ev_cb
          } else {
            ev_cb1 <- ev_cb[!rownames(ev_cb) %in% factors, , drop = FALSE]
          }


          #correlation1按照ev_cb提取
          correlation1 <- correlation1[rownames(ev_cb1), rownames(ev_cb1)]

          #连续变量的选择n=1
          bio_name1 <- corse_method(
            correlation = correlation1,
            importance = ev_cb1,
            vif = vif,
            n = n
          ) #下次模拟的变量
          #组合两种变量
          bio_name <- c(factors, bio_name1)
          #获取保留的变量存储路径的下标
          evlist <- c()
          for (i in seq_along(bio_name)) {
            evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i]) == T)
            if (length(evlist1) == 0) {
              evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".tif")[i]) == T)
            }
            evlist <- c(evlist, evlist1)
          }
          #计算两种变量的相关性
          correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
        }

        return(bio_name)

      } else {
        #分类变量大于1个的情况

        while (sum(correlation1 > r) > nrow(correlation1) |
               sum(correlation2 > r) > nrow(correlation2) | n == 0) {
          n <- n + 1

          #进行模拟获得重要性列表
          compare <- suppressMessages(
            maxent_single(
              x = x,
              evdir = evdir,
              evlist = evlist,
              factors = factors,
              nbg = nbg,
              bgwidth = NULL,
              mybgfile = mybgfile,
              args = args1,
              prodir = NULL,
              outdir = paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n),
              parallel = FALSE
            )
          )
          #下面根据上面模拟的结果删除相关性强的变量
          #变量重要性
          ev_cb <- read.csv(
            paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n, "/maxent/", sp_name, "/maxentResults.csv")
          ) %>%
            dplyr::select(., paste0(bio_name, ".contribution")) %>%
            utils::tail(., n = 1) %>% t() %>% as.data.frame()
          names(ev_cb) <- "value"

          ##排除贡献小于0.5的变量
          ev_cb <- ev_cb[ev_cb > 0.5, , drop = FALSE]

          #修改ev_cb的变量名
          nn <- c()
          for (i in 1:nrow(ev_cb)) {
            nn1 <- stringr::str_split_1(rownames(ev_cb)[i], ".contribution")[1]
            nn <- c(nn, nn1)
          }
          rownames(ev_cb) <- nn

          factors <- rownames(ev_cb)[rownames(ev_cb) %in% factors] #更新factors
          if (length(factors) == 0) {
            factors <- NULL
          }

          #将ev_cb按照变量类型拆开
          if (is.null(factors)) {
            ev_cb1 <- ev_cb
          } else {
            ev_cb1 <- ev_cb[!rownames(ev_cb) %in% factors, , drop = FALSE]
            ev_cb2 <- ev_cb[factors, , drop = FALSE]
          }

          #correlation1按照ev_cb提取
          if (is.null(factors)) {
            correlation1 <- correlation1[rownames(ev_cb1), rownames(ev_cb1)]
          } else {
            if (nrow(ev_cb1) == 0) {
              correlation2 <- correlation2[rownames(ev_cb2), rownames(ev_cb2), drop = FALSE]
            } else {
              correlation1 <- correlation1[rownames(ev_cb1), rownames(ev_cb1), drop = FALSE]
              correlation2 <- correlation2[rownames(ev_cb2), rownames(ev_cb2), drop = FALSE]
            }
          }

          #连续变量的选择n=1
          if (is.null(factors)) {
            if (nrow(ev_cb1) == 0) {
              bio_name1 <- NULL
            } else {
              bio_name1 <- corse_method(
                correlation = correlation1,
                importance = ev_cb1,
                vif = vif,
                n = n
              )
            }
            bio_name2 <- NULL
          } else{
            if (length(factors) == 1) {
              if (nrow(ev_cb1) == 0) {
                bio_name1 <- NULL
              } else {
                bio_name1 <- corse_method(
                  correlation = correlation1,
                  importance = ev_cb1,
                  vif = vif,
                  n = n
                )
              }
              bio_name2 <- factors
            } else {
              if (nrow(ev_cb1) == 0) {
                bio_name1 <- NULL
              } else {
                bio_name1 <- corse_method(
                  correlation = correlation1,
                  importance = ev_cb1,
                  vif = vif,
                  n = n
                )
              }  #下次模拟的变量
              #分类变量的选择
              bio_name2 <- corse_method(
                correlation = correlation2,
                importance = ev_cb2,
                vif = vif,
                n = n
              )
            }
          } #下次模拟的变量

          #将factors参数更新为bio_name2
          factors <- bio_name2
          #组合两种变量
          bio_name <- c(bio_name2, bio_name1)

          #获取保留的变量存储路径的下标
          evlist <- c()
          for (i in seq_along(bio_name)) {
            evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i]) == TRUE)
            if (length(evlist1) == 0) {
              evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".tif")[i]) == T)
            }
            evlist <- c(evlist, evlist1)
          }

          #计算两种变量的相关性
          if (is.null(factors)) {
            correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
            correlation2 <- as.data.frame(0)
          } else {
            if (is.null(bio_name1)) {
              correlation1 <- as.data.frame(0)
            } else {
              correlation1 <- abs(as.data.frame(cor(mybg1[bio_name1], method = cormethod)))
            }
            correlation2 <- abs(as.data.frame(cor(mybg2[bio_name2], method = cormethod)))
          }
          # utils::write.csv(correlation1, paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n, "/correlation1.csv"))
          # utils::write.csv(correlation2, paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n, "/correlation2.csv"))
        }

        return(bio_name)
      }
    } else {
      while (sum(correlation1 > r) > nrow(correlation1) |
             n == 0) {
        #bio_name初始值为最开始使用的全部变量
        n = n + 1

        #进行模拟获得重要性列表
        compare <- suppressMessages(
          maxent_single(
            x = x,
            evdir = evdir,
            evlist = evlist,
            factors = factors,
            nbg = nbg,
            bgwidth = NULL,
            mybgfile = mybgfile,
            args = args1,
            prodir = NULL,
            outdir = paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n),
            parallel = FALSE
          )
        )

        #下面根据上面模拟的结果删除相关性强的变量
        #变量重要性
        ev_cb <- read.csv(
          paste0(outdir, "/TBlabENMtemp", random_num, "/", fc1, rm1, n, "/maxent/", sp_name, "/maxentResults.csv")
        ) %>%
          dplyr::select(., paste0(bio_name, ".contribution")) %>%
          utils::tail(., n = 1) %>% t() %>% as.data.frame()
        names(ev_cb) <- "value"

        ##排除贡献小于0.5的变量
        ev_cb <- ev_cb[ev_cb > 0.5, , drop = FALSE]

        #修改ev_cb的变量名
        nn <- c()
        for (i in 1:nrow(ev_cb)) {
          nn1 <- str_split_1(rownames(ev_cb)[i], ".contribution")[1]
          nn <- c(nn, nn1)
        }
        rownames(ev_cb) <- nn
        #correlation1按照ev_cb提取
        correlation1 <- correlation1[rownames(ev_cb), rownames(ev_cb)]
        bio_name <- corse_method(
          correlation = correlation1,
          importance = ev_cb,
          vif = vif,
          n = n
        )  #下次模拟的变量
        evlist <- c()
        for (i in seq_along(bio_name)) {
          evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".asc")[i]) == T)
          if (length(evlist1) == 0) {
            evlist1 <- which(stringr::str_detect(biolistall, paste0(bio_name, ".tif")[i]) == T)
          }
          evlist <- c(evlist, evlist1)
        }
        correlation1 <- abs(as.data.frame(cor(mybg[bio_name], method = cormethod)))

      }
      return(bio_name)
    }
  }

  fun4 <- function(i) {
    #设置好fc和rm后进行模拟保留最重要变量i=1
    var <- fun2(fc1 = combin[i, 1] , rm1 = combin[i, 2])
    df[i, 3] <- paste0(var, collapse = ",")
  }

  #####################################################################
  if (is.null(myenv) == FALSE) {
    cat("\nThe environment variable has been specified and is no longer selected!\n")
    evlist <- NULL
  }

  tzhs <- c("L", "Q", "H", "P", "T")
  if (is.null(outdir)) {
    outdir <- "."
  }
  #dir.create(paste0(outdir, "/TBlabENM"),recursive = TRUE, showWarnings = FALSE)
  mybgfile_rmd <- mybgfile
  factors123 <- factors
  biolistall <- list.files(evdir, pattern = ".asc$|.tif$", full.names = TRUE)

  #获取物种名 对路径拆分并取倒数第一个字符串
  spname1 <- stringr::str_split_1(x, "/")[length(stringr::str_split_1(x, "/"))]
  sp_name <- stringr::str_split_1(spname1, ".csv$")[1]
  cat(paste("Select optimum parameters for", sp_name, "\n"))
  dir.create(paste0(outdir, "/maxent/", sp_name),
             showWarnings = FALSE,
             recursive = TRUE)
  #################################################################################
  #当指定环境变量时不再进行变量的选择
  if (is.null(myenv)) {
    #变量名称
    biolist <- list.files(evdir, pattern = ".asc$|.tif$", full.names = TRUE)
    if (length(biolist) == 0) {stop("The environment variable does not exist. Check whether it is in `asc` format.\n")}
    if (is.null(evlist) == FALSE) {
      biolist <- biolist[evlist]
    }
    bio_name <- c()
    for (i in seq_along(biolist)) {
      bioname1 <- stringr::str_split_1(biolist[i], "/")[length(stringr::str_split_1(biolist[i], "/"))]
      bio_name0 <- stringr::str_split_1(bioname1, ".asc|.tif")[1]
      bio_name <- c(bio_name, bio_name0)
    }
    bio_name_all <- bio_name
    cat(paste("Initial variable:", paste(bio_name, collapse = ","), "\n"))

    #判断分类变量factors是否包含在给定的环境数据集内
    if (is.null(factors)) {
      cat("All variables are continuous variables \n")
    } else{
      factors1 <- factors
      factors <- factors[factors %in% bio_name]
      if (length(factors) == 0) {
        stop(
          "No categorical variable is found in the given list of variables, and parameter 'factors' is invalid. Please check whether parameter 'evdir','evlist','factors' are set correctly."
        )
      }
      if (length(factors) != length(factors1)) {
        stop(
          paste0(
            "Only variable '",
            factors,
            "' is identified as a categorical variable. Please check whether other categorical variables are in the given set of environment variables"
          )
        )
      }
      cat(paste(
        "Variable:",
        paste(factors, collapse = ","),
        "as categorical variable\n"
      ))
    }
    cat("***********Selecting parameters***********\n")
    #提取发生数据的环境值
    biostack <- terra::rast(biolist)
    occ <- utils::read.csv(x) #读取物种坐标数据
    occ <- occ[c(2, 3)]
    names(occ) <- c("x", "y")
    if (is.null(bgwidth) == FALSE & is.null(mybgfile)) {
      #将发生点转为sf对象
      occs.sf <- sf::st_as_sf(occ, coords = c("x","y"), crs = terra::crs(biostack))
      #转为等面积投影
      eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
      occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)
      #创建缓冲区
      occs.buf <- sf::st_buffer(occs.sf, dist = bgwidth) |>
        sf::st_union() |>
        sf::st_sf() |>
        sf::st_transform(crs = terra::crs(biostack)) #再转为经纬度投影
      #将环境变量裁剪至缓冲区
      biostack <- terra::crop(biostack, occs.buf)
      biostack <- terra::mask(biostack, occs.buf)
      cat(paste0("Crop the environment variable to a buffer with a radius of", bgwidth/1000, "km centered on the point of occurrence\n") )
      }

    occdata <- terra::extract(biostack, occ, ID = FALSE)
    n_na <- nrow(occdata) - nrow(occ)
    #提取背景值并计算变量相关性
    ##随机生成10000个点
    if (is.null(mybgfile)) {

      mybg0 <- terra::spatSample(biostack, nbg, na.rm = T, xy = T)
      mybgfile <- mybg0[1:2]
      write.csv(mybgfile,
                paste0(outdir, "/maxent/", sp_name, "/bg.csv"),
                row.names = FALSE)
      mybg <- mybg0[-(1:2)]
      cat(paste0(nbg ," background points are randomly generated.\n"))
    } else{
      cat("The background points has been specified!\n")
      mybg <- terra::extract(biostack, mybgfile, ID = FALSE)
    }
    ##判断是否存在分类变量，若存在则分开计算
    if (is.null(factors) == FALSE) {
      #mybg分为两组，mybg1为连续变量，mybg2为分类变量
      mybg1 <- mybg[, bio_name[!bio_name %in% factors], drop = FALSE]
      mybg2 <- mybg[, factors, drop = FALSE]

      #绘制相关性热图
      correlation1 <- cor(mybg1, method = cormethod)
      correlation2 <- cor(mybg2, method = cormethod)
      if (nrow(correlation1) > 1) {
      grDevices::jpeg(filename = paste0(outdir, "/maxent/", sp_name, "/cor_continuous.jpg"),
                      width = 15 + nrow(correlation1), height = 15 + nrow(correlation1), units = "cm", res = 300)
      corrplot::corrplot.mixed(correlation1, tl.pos = c( "lt"), tl.col = "black",
                               diag = c("u"), title = "continuous variable")
      dev.off()}
      #
      if (nrow(correlation1) > 1) {
      grDevices::jpeg(filename = paste0(outdir, "/maxent/", sp_name, "/cor_categorical.jpg"),
                      width = 15 + nrow(correlation2), height = 15 + nrow(correlation2), units = "cm", res = 300)
      corrplot::corrplot.mixed(correlation2, tl.pos = c( "lt"), tl.col = "black",
                               diag = c("u"), title = "categorical variable")
      dev.off()}
  #
      correlation1 <- abs(as.data.frame(cor(mybg1, method = cormethod)))
      correlation2 <- abs(as.data.frame(cor(mybg2, method = cormethod)))
    } else {
      ##绘制相关性热图
      correlation1 <- cor(mybg, method = cormethod)
      grDevices::jpeg(filename = paste0(outdir, "/maxent/", sp_name, "/cor_continuous.jpg"),
                      width = 15 + nrow(correlation1), height = 15 + nrow(correlation1), units = "cm", res = 300)
      corrplot::corrplot.mixed(correlation1, tl.pos = c( "lt"), tl.col = "black", diag = c("u"))
      dev.off()
      correlation1 <- abs(as.data.frame(cor(mybg, method = cormethod)))
    }

    #组合fc和 rm
    fc <- toupper(fc)
    combin <- expand.grid(fc, rm, stringsAsFactors = FALSE) %>%
      dplyr::mutate(fc = purrr::map_chr(
        .x = Var1,
        .f = function(x) {
          paste0(tzhs[tzhs %in% stringr::str_split_1(x, "")], collapse = "")
        }
      ))
    combin <- combin[c(3, 2)]
    df <- data.frame(matrix(NA, nrow(combin), 1))
    df <- cbind(combin, df)
    names(df) <- c("fc", "rm", "env")

    #根据组合设置 args参数

    args1 <- TBlabENM::maxent_args(
      l = FALSE,
      q = FALSE,
      p = FALSE,
      h = FALSE,
      t = FALSE,
      responsecurves = FALSE,
      jackknife = FALSE,
      pictures = FALSE
    )
    ##并行计算
    if (parallel == T) {
      # 开启集成
      snowfall::sfInit(parallel = TRUE, cpus = ncpu)
      # 注册每个环境变量
      snowfall::sfExport("fun2")
      snowfall::sfExport("random_num")
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
fit <- try(  #报错调试
      k <- snowfall::sfLapply(1:nrow(combin), fun4)
      )
if ('try-error' %in% class(fit)) {snowfall::sfStop()
  unlink(paste0(outdir, "/TBlabENMtemp", random_num) , recursive = T)
  }
      snowfall::sfStop()  # 关闭集群

      df$env <- unlist(k)

    } else{
fit <- try(  #报错调试
      for (i in 1:nrow(combin)) {
        #设置好fc和rm后进行模拟保留最重要变量i=1
        var <- fun2(fc1 = combin[i, 1] , rm1 = combin[i, 2])
        df[i, 3] <- paste0(var, collapse = ",")
      }
      )
      if ('try-error' %in% class(fit)) {unlink(paste0(outdir, "/TBlabENMtemp", random_num) , recursive = T)}
    }
    #删除缓存文件
    unlink(paste0(outdir, "/TBlabENMtemp", random_num) , recursive = T)
  } else {
    biolist <- list.files(evdir, pattern = ".asc$|.tif$", full.names = TRUE)
    #判断分类变量factors是否包含在给定的环境数据集内
    ##全部变量名称
    bio_name_all <- c()
    for (i in seq_along(biolistall)) {
      bioname1 <- stringr::str_split_1(biolistall[i], "/")[length(stringr::str_split_1(biolistall[i], "/"))]
      bio_name0 <- stringr::str_split_1(bioname1, ".asc|tif")[1]
      bio_name_all <- c(bio_name_all, bio_name0)
    }

    #指定的变量
    bio_name <- myenv

    if (sum(bio_name %in% bio_name_all) != length(bio_name)) {
      stop(
        "One or more specified variables do not exist! Pleaes check parameter 'evdir' and 'myenv'."
      )
    } else {
      cat(paste("Initial variable:", paste(bio_name, collapse = ","), "\n"))
    }
    if (is.null(factors)) {
      cat("All variables are continuous variables \n")
    } else{
      factors1 <- factors
      factors <- factors[factors %in% bio_name]
      if (length(factors) == 0) {
        stop(
          "No categorical variable is found in the given list of variables, and parameter 'factors' is invalid. Please check whether parameter 'evdir','evlist','factors' are set correctly."
        )
      }
      if (length(factors) != length(factors1)) {
        stop(
          paste0("Only variable '", factors,
            "' is identified as a categorical variable. Please check whether other categorical variables are in the given set of environment variables."
          )
        )
      }
      cat(paste(
        "Variable:",
        paste(factors, collapse = ","),
        "as categorical variable\n"
      ))
    }

    #提取发生数据的环境值
    ##获取bio_name的下标
    xb <- c()
    for (i in bio_name) {
      xb1 <- which(stringr::str_detect(biolistall, paste0(i, ".asc")) == T)
      if (length(xb1) == 0) {
        xb1 <- which(stringr::str_detect(biolistall, paste0(i, ".tif")) == T)}
      xb <- c(xb, xb1)
    }
    biostack <- terra::rast(biolistall[xb])
    occ <- utils::read.csv(x) #读取物种坐标数据
    occ <- occ[c(2, 3)]
    occdata <- terra::extract(biostack, occ, ID = FALSE)

    #提取背景值并计算变量相关性
    ##随机生成10000个点
    if (is.null(mybgfile)) {
      mybg <- terra::spatSample(biostack, nbg, na.rm = T, xy = T)[-(1:2)]
    } else{
      mybg <- terra::extract(biostack, mybgfile, ID = FALSE)
    }

    #组合fc和 rm
    fc <- toupper(fc)
    combin <- expand.grid(fc, rm, stringsAsFactors = FALSE) %>%
      mutate(fc = purrr::map_chr(
        .x = Var1,
        .f = function(x) {
          paste0(tzhs[tzhs %in% stringr::str_split_1(x, "")], collapse = "")
        }
      ))
    combin <- combin[c(3, 2)]
    df <- data.frame(matrix(NA, nrow(combin), 1))
    df <- cbind(combin, df)
    df[, 3] <- rep(paste0(myenv, collapse = ","), nrow(df))
    names(df) <- c("fc", "rm", "env")

  }

  if (is.null(opt)) {
    parameter <- df
  } else{
    #剩余的组合进行evaluate通过AICC确定最佳组合
    df <- df %>%
      dplyr::mutate(occdata = purrr::map(
        .x = env,
        .f = function(x) {
          #获取变量下标
          xb <- c()
          xb2 <- stringr::str_split_1(x, ",")
          for (i in seq_along(xb2)) {
            xb1 <- which(stringr::str_detect(biolist, paste0(xb2, ".asc")[i]) == T)
            if (length(xb1) == 0) {
              xb1 <- which(stringr::str_detect(biolist, paste0(xb2, ".tif")[i]) == T)
            }
            xb <- c(xb, xb1)
          }

          occdata <- cbind(occ, occdata[xb2]) #选择变量后添加xy坐标
        }
      )) %>%
      dplyr::mutate(bgdata = purrr::map(
        .x = env,
        .f = function(x) {
          #获取变量下标
          xb <- c()
          xb2 <- stringr::str_split_1(x, ",")
          for (i in seq_along(xb2)) {
            xb1 <- which(stringr::str_detect(biolist, paste0(xb2, ".asc")[i]) == T)
            if (length(xb1) == 0) {
              xb1 <- which(stringr::str_detect(biolist, paste0(xb2, ".tif")[i]) == T)}
            xb <- c(xb, xb1)
          }
          bgdata <- cbind(mybgfile, mybg[xb2]) #选择变量后添加xy坐标
        }
      )) %>% #计算变量个数
      mutate(num = purrr::map_dbl(
        .x = occdata,
        .f = function(x) {
          ncol(x)
        }
      ))
    #将分类变量转为因子型
    for (i in 1:nrow(df)) {
      dff4 <- df[[i,4]]
      dff5 <- df[[i,5]]
      categoricals <- factors123[factors123 %in% names(dff4)]
      if (length(categoricals) > 0) {
      num <- which(names(dff4) %in% factors123 == T)
      for (j in num) {
        dff4[,j] <- factor(dff4[,j] )
        dff5[,j] <- factor(dff5[,j] )
      }
      df[[i,4]] <- dff4
      df[[i,5]] <- dff5
                                     }
    }

    #判断保留的变量有几个，只有一个则无法调优
    df1 <- dplyr::filter(df, num > 1)
    if (nrow(df1) == 0) {
      #删除缓存文件
      unlink(paste0(outdir, "/TBlabENMtemp", random_num) , recursive = T)
      stop("All parameter combinations end up preserving only one variable.")
    }

    #使用pmap函数并行评估（由于ENMevaluate函数的参数大于2个，所以使用pmap函数）
    #设置参数
    block <- ENMeval::get.randomkfold(occ, mybgfile, kfolds = 5)
    if (nrow(occdata) >= 25) {
      #partitions = "randomkfold"
      #partition.settings = list(kfolds = 10)
      partitions = "user"
      user.grp = block
    } else {
      partitions = "jackknife"
      partition.settings = NULL
    }

    gz <- purrr::pmap(
      df,
      .f = function(fc, rm, occdata, bgdata, ...) {
        e <- ENMeval::ENMevaluate(
          occs = occdata,
          bg = bgdata,
          tune.args = list(fc = fc, rm = rm),
          #partitions = "jackknife", #数据分区方式，有2+6种
          partitions = partitions,
          user.grp = block,
          #数据分区方法|测试集、训练集划分方法
         # partition.settings =  partition.settings,
          #数据分区设置
        #  other.settings =  list(
        #    abs.auc.diff = TRUE,
         #   pred.type = "logistic",
        #    validation.bg = "full",
            #与partitions参数分区类型对应
        #    other.args = NULL
       #   ),
          #其他额外设置，有默认值
          taxon.name = sp_name,
          n.bg = nbg,
          algorithm = "maxent.jar",
          #使用的模型，有三种
          overlap = FALSE,
          #生态位重叠
         # categoricals = factors123[factors123 %in% names(occdata)],
          #指定分类变量,"IAWC_CLASS", , "T_USDA_TEX_CLASS"
          doClamp = FALSE
        )
        res <- ENMeval::eval.results(e)
      }
    )

    #评估的一系列参数
    gz1 <- as.data.frame(data.table::rbindlist(gz))
    #合并df和gz1
    cs <- cbind(df[1:3], gz1[c(1:16, 19)])
    cs <- cs[-(4:5)]
    cs <- dplyr::arrange(cs, AICc, auc.diff.avg, or.mtp.avg)
    cs1 <- cs
    #选择最佳模型auc.train, cbi.train, auc.diff.avg, auc.val.avg, cbi.val.avg, or.10p.avg, or.mtp.avg, AICc
    #单一指标法，顺序指标法，最大指标法。#opt1为最佳模型
    for (i in opt) { #i = opt[1]
      if (is.na(min(cs[i]))) {
        warning("'", paste0(i), "'", " is NA in one or more parameter combinations, skip it.")
        next}

      if (i %in% c("auc.diff.avg", "or.10p.avg", "or.mtp.avg", "AICc")) {
        opt1 <- cs[which(min(cs[i]) == cs[[i]]),]
        cs <- opt1
      } else {
        opt1 <- cs[which(max(cs[i]) == cs[[i]]),]
        cs <- opt1
      }
    }
    if (exists("opt1") == FALSE) {
      opt1 <- dplyr::filter(cs, auc.val.avg == max(auc.val.avg))
      opt <- append(opt, "auc.val.avg")
      warning("'", paste0(opt), "'", " is NA in one or more parameter combinations, use 'auc.val.avg' instead.")
    }

    #保存结果
    #删除缓存文件
    #unlink(paste0(outdir, "/TBlabENMtemp", random_num) , recursive = T)
    cat(paste0(outdir, "/maxent/", sp_name, "/tuneparameter.csv"))
    utils::write.csv(cs1,
                     paste0(outdir, "/maxent/", sp_name, "/tuneparameter.csv"),
                     row.names = FALSE)

    #绘制模型调优图
    p_tun(cs1, opt)

    #最佳模型的参数df_best
    df_best <- dplyr::filter(df, fc == opt1$fc, rm == opt1$rm)
    env_best <- stringr::str_split_1(df_best$env, ",")
    env_best_f <- env_best[env_best %in% factors123] #最优模型保留的分类变量
    env_best_c <- env_best[!env_best %in% env_best_f] #最优模型保留的连续变量
    #绘制最佳模型的相关性热图
    mybg <- as.data.frame(tidyr::unnest(df_best[5], cols = c("bgdata")))
    mybg1 <- mybg[env_best_c]
    mybg2 <- mybg[env_best_f]
    correlation1 <- cor(mybg1, method = cormethod)
    correlation2 <- cor(as.data.frame(lapply(mybg2,as.numeric)), method = cormethod)
    if (nrow(correlation1) > 1) {
    grDevices::jpeg(filename = paste0(outdir, "/maxent/", sp_name, "/cor_continuous_best.jpg"),
                    width = 15 + nrow(correlation1), height = 15 + nrow(correlation1), units = "cm", res = 300)
    corrplot::corrplot.mixed(correlation1, tl.pos = c( "lt"), tl.col = "black",
                             diag = c("u"), title = "continuous variable")
    dev.off()
    }
    if (nrow(correlation2) > 1) {
    grDevices::jpeg(filename = paste0(outdir, "/maxent/", sp_name, "/cor_categorical_best.jpg"),
                    width = 15 + nrow(correlation2), height = 15 + nrow(correlation2), units = "cm", res = 300)
    corrplot::corrplot.mixed(correlation2, tl.pos = c( "lt"), tl.col = "black",
                             diag = c("u"), title = "categorical variable")
    dev.off()}
    #零模型检验
    #使用最佳模型的参数重新构建模型测试获得参数e.mx

    if (null_model == TRUE) {
      e.mx <- ENMeval::ENMevaluate(
        occs = df_best$occdata,
        bg = df_best$bgdata,
        tune.args = list(fc = df_best$fc, rm = df_best$rm),
        #partitions = "jackknife", #数据分区方式，有2+6种
        partitions = "randomkfold",
        #其他额外设置，有默认值
        taxon.name = sp_name,
        n.bg = nbg,
        algorithm = "maxent.jar",
        #使用的模型，有三种
        overlap = FALSE,
        #生态位重叠
       # categoricals = env_best_f,
        doClamp = FALSE
      )
      cat("\n**************null model test*****************\n")
      mod.null <- ENMnulls(
        e = e.mx,
        mod.settings = list(fc = df_best$fc, rm = df_best$rm),
        #最佳模型的参数
        no.iter = 100,
        eval.stats = c("auc.val", "or.10p"),
       # user.eval.type = "knonspatial",
        parallel = parallel,
        numCores = ncpu,
        parallelType = "doSNOW",
        quiet = TRUE
      )
      p <- ENMeval::evalplot.nulls(
        mod.null,
        stats = c("auc.val", "or.10p"),
        plot.type = "histogram",
        return.tbl = FALSE
      )
      ggsave(
        filename = "null_model.jpg",
        plot = p,
        path = paste0(outdir, "/maxent/", sp_name),
        width = 24,
        height = 18,
        units = c("cm"),
        create.dir = TRUE
      )
    }

    parameter <- opt1[1, 1:3]
    parameter$species <- sp_name
    parameter$number <- nrow(occ)
    parameter <- parameter[, c(4, 5, 3, 1, 2)]

  }

  ###rmarkdown file
  wendang <- function(x, n_na, occdata, sp_name, evdir, bio_name_all, factors, mybgfile_rmd, nbg,
                      fc, rm, r, cormethod, vif, vifth, opt, bestpar, outdir,env_best_f,env_best_c, null_model){
    Sys.setenv(RSTUDIO_PANDOC = Sys.getenv("RSTUDIO_PANDOC"))
    rmarkdown::render(
     input = system.file("extdata", "Models_detail.Rmd", package = "TBlabENM"),
     #input = "C:/Users/why/TBlabENM/inst/extdata/Models_detail.Rmd",
                      output_file = paste0(outdir, "/maxent/", sp_name, "/Models_tuning.html"),
                      quiet = TRUE,
                      params = list("sp_name" = sp_name, "x" = x, "occdata" = occdata, "n_na" = n_na,
                                    "evdir" = evdir, "bio_name_all" = bio_name_all,
                                    "factors" = factors, "mybgfile_rmd" = mybgfile_rmd, "nbg" = nbg,
                                    "fc" = fc, "rm" = rm, "r" = r, "cormethod" = cormethod, "vif" = vif,
                                    "vifth" = vifth, "opt" = opt, "bestpar" = bestpar, "outdir" = outdir,
                                    "env_best_f" = env_best_f, "env_best_c" = env_best_c,
                                    "null_model" = null_model))
  }

  wendang(sp_name = sp_name, x = x, occdata = occdata, n_na = n_na,
          evdir = evdir, bio_name_all = bio_name_all,
          factors = factors123,
          mybgfile_rmd = mybgfile_rmd, nbg = nbg, fc = fc,
          rm = rm, r = r, cormethod = cormethod, vif = vif, vifth = vifth,
          opt = opt, bestpar = opt1, outdir = outdir,env_best_f = env_best_f, env_best_c = env_best_c,
          null_model = null_model)

  cat("\n")
  return(parameter)
}
