#' Grid centroid calculation
#' @description
#' 计算适生区栅格的质心. 栅格应为二元栅格, 即1代表适生区, 0代表非适生区.
#'
#' @param radir 栅格所在文件夹的路径, 可以是多个文件夹.
#' @param format 栅格格式
#' @param outdir 输出文件路径
#' @param parallel 是否并行
#' @param ncpu 并行时使用的cpu数
#' @param proname 逻辑值, 栅格所在文件夹的名称是否为投影时期名称. 如果为T, 则生成的数据框包含投影时期.
#'
#' @return
#' 数据框
#' @export
#'
#' @examples
#' df <- ENMcentroid(
#'   radir = "F:/sdm/maxent/32/1/p",
#'   proname = F,
#'   format = "tif",
#'   parallel = T,
#'   ncpu = 2
#' )
#' df
#'
ENMcentroid <- function(radir, proname = F, format = "tif", outdir = NULL, parallel = F, ncpu = 2) {
  # 读取栅格数据
  ralist <- list.files(radir, pattern = paste0(format, "$"), full.names = TRUE)
  if (length(ralist) == 0) {
    stop("No grid file is found./n")
  }

  if (proname == FALSE) {
    radf <- tibble::as_tibble(ralist) %>%
      dplyr::mutate(spname = purrr::map_chr(.x = value, .f = function(x) {
        stringr::str_split_1(basename(x), paste0(".", format))[1]
      }))
  } else {
    radf <- tibble::as_tibble(ralist) %>%
      dplyr::mutate(spname = purrr::map_chr(.x = value, .f = function(x) {
        stringr::str_split_1(basename(x), paste0(".", format))[1]
      })) %>%
      dplyr::mutate(proname = basename(dirname(ralist)))
  }

  # 重心（质心）计算公式
  centroid <- function(x, format) {
    # 读取单个栅格数据
    ra <- terra::rast(x)
    names(ra) <- stringr::str_split_1(basename(x), paste0(".", format))[1]
    # 将栅格转为数据框
    df <- terra::as.data.frame(ra, xy = T)
    names(df) <- c("x", "y", "value")
    gi <- df$value
    xi <- df$x
    yi <- df$y
    x <- sum(gi * xi, na.rm = TRUE) / sum(gi, na.rm = TRUE)
    y <- sum(gi * yi, na.rm = TRUE) / sum(gi, na.rm = TRUE)
    xy <- tibble::tibble(x, y)
    xy$spname <- names(ra)
    return(xy)
  }

  if (parallel == TRUE) {
    snowfall::sfInit(parallel = TRUE, cpus = ncpu)
    snowfall::sfLibrary(purrr)
    results <- tryCatch(
      {
        snowfall::sfLapply(ralist, centroid, format = format)
      },
      error = function(e) {
        message("Error in parallel processing: ", e$message)
        snowfall::sfStop()
        stop(e)
      }
    )
    # 合并所有线程返回的数据框
    final_df <- dplyr::bind_rows(results)
    snowfall::sfStop()
  } else {
    final_df <- tibble::tibble(x = numeric(), y = numeric(), species = character())
    for (i in seq_along(ralist)) {
      df1 <- centroid(ralist[i])
      final_df <- rbind(final_df, df1)
    }
  }
  # final_df与radf合并
  final_df <- dplyr::left_join(radf, final_df)
  final_df <- final_df[2:length(final_df)]
  # 保存结果
  if (is.null(outdir)) {
    return(final_df)
  } else {
    write.csv(final_df, outdir, row.names = FALSE)
    return(final_df)
  }
}
