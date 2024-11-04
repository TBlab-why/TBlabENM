#' @title Batch spatially thin species occurence data
#' @description 批量对物种发生数据进行空间过滤.
#' @details
#' 本函数利用\code{\link[spThin]{thin}}对物种发生数据进行地理空间过滤来降低采样偏差，并且额外
#'     生成一个连接文件(compare开头),指明了与原始数据相比哪些记录被保留了下来.
#'
#'
#' @param spdir 要过滤的物种文件的路径的集合，包括文件名（.csv）.
#' @param thin.par 数值型，地理空间过滤距离（km）.
#' @param outdir 结果保存路径，未指定时在当前工作路径自动生成保存路径.
#' @param spec.col 物种名所在的列名.
#' @param long.col 经度所在的列名.
#' @param lat.col 纬度所在的列名.
#'
#' @return 一个文件夹,包含每个物种过滤后的csv文件和过滤前后的连接文件.
#' @export
#'
#' @examples
#' #读取物种路径列表
#' dir <- system.file("extdata", "species", package = "TBlabENM")
#' splist <- list.files(dir, pattern = ".csv$", full.names = T)
#' ENMspthin(spdir = splist,
#'           spec.col = species,
#'           long.col = longitude,
#'           lat.col = latitude,
#'           thin.par = 10,
#'           outdir = NULL)
ENMspthin <- function(spdir, spec.col, long.col, lat.col, thin.par, outdir = NULL) {
  if(is.null(outdir)){outdir = "."}
  dir.create(paste0(outdir, "/TBlabENM/occthin", thin.par, "km"), recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(spdir)) {
    occdata <- utils::read.csv(spdir[i], fileEncoding = "GB18030")
    name <- stringr::str_split_1(spdir[i], pattern = "/")[length(stringr::str_split_1(spdir[i], pattern = "/"))] %>%
      stringr::str_split_1(., pattern = ".csv$")
    if(file.exists(paste0(outdir, "/TBlabENM/occthin", thin.par, "km/", name[1], "_thin1.csv"))){
      file.remove(paste0(outdir, "/TBlabENM/occthin", thin.par, "km/", name[1], "_thin1.csv"))
    }

      spThin::thin(
      loc.data = occdata,
      spec.col = spec.col,
      long.col = long.col,
      lat.col = lat.col,
      thin.par = thin.par,
      reps = 5,
      locs.thinned.list.return = FALSE,
      write.files = TRUE,
      max.files = 1,
      out.dir= paste0(outdir, "/TBlabENM/occthin", thin.par, "km"),
      out.base =  name[1],
      write.log.file = TRUE,
      log.file = "spatial_thin_log.txt",
      verbose = FALSE )
  thinfile <- utils::read.csv(paste0(outdir, "/TBlabENM/occthin", thin.par, "km/", name[1], "_thin1.csv") )
  thinfile$retain <- rep("T", nrow(thinfile))
  compare <-suppressMessages(dplyr::left_join(occdata, thinfile))

  utils::write.csv(compare, paste0(outdir, "/TBlabENM/occthin", thin.par, "km/compare_", name[1], ".csv"),
            fileEncoding = "GB18030",na = "")
    }

  }
