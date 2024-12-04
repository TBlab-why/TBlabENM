
#' Download species data from GBIF and perform routine cleaning
#' @description
#' 使用rgbif包从GBIF下载物种的完整发生数据，并进行常规数据清洗用于物种分布建模. 如果对数据有更多更精细的要求, 建议使用rgbif包进行数据下载.
#' @details
#' 本函数主要参考[https://docs.ropensci.org/rgbif/articles/getting_occurrence_da.html](https://docs.ropensci.org/rgbif/articles/getting_occurrence_da.html)的流程下载和清洗数据. 主要过程如下：
#'
#' **1. rgbif包下载数据**
#' - pred_default() #删除默认的地理空间问题, 仅保留带有坐标的记录, 仅保留存在的记录, 移除化石和活体标本记录. See \code{\link[rgbif]{pred_default}}.
#' - pred_or(pred_lt("coordinateUncertaintyInMeters", uncertainty), pred_isnull("coordinateUncertaintyInMeters")) #对不确定性范围进行筛选并保留为NULL的记录.
#' - pred_or(pred_not(pred_in("establishmentMeans", c("MANAGED","INTRODUCED"))), pred_isnull("establishmentMeans")) #删除动植物园栽培或引进的记录并保留为NULL的记录.
#'
#' **2. CoordinateCleaner包清洗数据**
#' - cc_cen(buffer = 5000)  # remove country centroids within 5km. See \code{\link[CoordinateCleaner]{cc_cen}}.
#' - cc_cap(buffer = 5000)  # remove capitals centroids within 5km. See \code{\link[CoordinateCleaner]{cc_cap}}.
#' - cc_inst(buffer = 5000)  # remove zoo and herbaria within 5km. See \code{\link[CoordinateCleaner]{cc_inst}}.
#' - cc_sea() #remove from ocean. See \code{\link[CoordinateCleaner]{cc_sea}}.
#' - cc_val()  #removes or flags non-numeric and not available coordinates. See \code{\link[CoordinateCleaner]{cc_val}}.
#' - cc_equ()  #删除经纬度相等的坐标. See \code{\link[CoordinateCleaner]{cc_equ}}.
#' - cc_zero(buffer = 0.5)  #删除经度或纬度为零及其半径0.5°内的记录. See \code{\link[CoordinateCleaner]{cc_zero}}.
#'
#' **3. dplyr包筛选数据**
#' - filter(year > 1950) #只选择1950年后的记录.
#' - filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) #删除坐标精度小于0.01的点并保留为NA的点.
#' - filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) #删除已知有错误的坐标不确定度.
#'
#' **4. rWCVP包过滤非原产地的分布点**
#'
#' 如果rWCVP == TRUE还将使用rWCVP包过滤非原产地的记录. 更多细节参考[https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html](https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html).

#' @param spnames character string or character vector. 要下载物种的学名.
#' @param range character string. 用国家或地区代码(ISO 3166-1 alpha2)或WKT多边形(逆时针)表示的地理范围. 如果为NULL，则在全球范围下载.
#' @param rWCVP logical. 是否使用rWCVP包提供的原产地范围进行过滤.
#' @param outdir character string. 输出文件路径.
#' @param filename character string. 输出文件名. 不包括后缀, 格式默认为.csv.
#' @param user character string. GBIF账号.
#' @param pwd character string. GBIF密码.
#' @param email character string. GBIF注册邮箱.
#' @param dist numerical. 缓冲区范围, 单位°。仅rWCVP == TRUE时生效. 位于缓冲区范围内的非原产地坐标点将被保留, 需要注意当位于不同纬度时1°代表的实际地理距离不同.
#' @param uncertainty numerical. 不确定性范围, 单位m. 当coordinateUncertaintyInMeters列已经评估时允许的最大不确定性范围.

#' @importFrom rgbif occ_download pred_default  pred_in  pred_or  pred pred_lt pred_not
#' @importFrom rgbif pred_isnull occ_download_wait occ_download_get occ_download_import
#' @importFrom purrr map_int
#' @importFrom dplyr filter distinct select rename arrange
#' @importFrom stats setNames
#' @importFrom CoordinateCleaner cc_cen cc_cap cc_inst cc_sea cc_val cc_equ cc_coun cc_zero
#' @importFrom sf st_as_sf st_crs st_intersects st_union st_buffer
#' @importFrom rWCVP wcvp_distribution
#' @importFrom utils write.csv
#'
#' @references [https://docs.ropensci.org/rgbif/articles/getting_occurrence_da.html](https://docs.ropensci.org/rgbif/articles/getting_occurrence_da.html)
#'
#' [https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html](https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html)
#'
#' @return
#' 从gbif下载的原始数据压缩包及经过清洗的csv格式的物种发生数据集.
#' @export
#'
#' @examples
#' \donttest{ENMgbif_download(spnames = c("Phoebe formosana", "Decaisnea insignis", "Ginkgo biloba"),
#'                  range = "POLYGON((72.2 56, 137.5 56, 137.5 0, 72.2 0, 72.2 56))", #按给定的地理范围下载
#'                # range = "CN", #按国家或地区代码下载
#'                # range = NULL, #全球区域下载
#'                  uncertainty = 10000,
#'                  rWCVP = TRUE,
#'                  dist = 1,
#'                  outdir = "./gbif",
#'                  filename = "tidy_gbif",
#'                  user = "Your user",
#'                  pwd = "Your pwd",
#'                  email = "Your email")
#'               }

ENMgbif_download <- function(spnames, range = NULL, uncertainty = 10000, rWCVP = TRUE, dist = 1,
                          outdir = NULL, filename, user, pwd, email){
  if (is.null(outdir)) {outdir <- "./TBlabENM/gbif_data"} else {
    outdir <- paste0(outdir, "/TBlabENM/gbif_data")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(range) == FALSE) {
    if ( range == "TW") {
      stop("Please use 'CN'.")}
    if ( range == "MO") {
      stop("Please use 'CN'.")}
    if ( range == "HK") {
      stop("Please use 'CN'.")}
  }

  ##详细用法见https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
  ####首先从GBIF下载初步筛选的记录。由于做分布区相关研究，记录只选择来自标本馆的，并做相应筛选####
  #0.1获取需要下载物种的识别号(key)
  keylist <- spnames %>%
    as.data.frame()  %>%
    mutate(key = map_int(.x = ., .f = function(x){
      key <- rgbif::name_suggest(q = x, rank = 'species')$data$key[1]
     # print(key)
      if (is.null(key)) {
        key <- NA
      }
      key
  }))
  names(keylist)[1] <- "species"
  if (is.na(sum(keylist$key))) {
    matter <- filter(keylist, is.na(key))
    print(matter)
    stop("No key was found for the species above in GBIF.")}
  #0.2 数据下载准备。并获取下载doi，用于文章引用，用于论文必须要有doi.运行后需要等待几分钟让GBIF准备数据
  #全球范围内下载
  if (is.null(range)) {
    gbif_download <- occ_download(
      pred_default(), #默认过滤规则，见帮助文档
      pred_in("taxonKey", c(keylist$key)), #要下载物种的key
      pred_or(
        pred_lt("coordinateUncertaintyInMeters", uncertainty), pred_isnull("coordinateUncertaintyInMeters"),
        pred_not(pred_in("establishmentMeans", c("MANAGED","INTRODUCED"))), pred_isnull("establishmentMeans")
      ), #排除栽培的或引进的记录
      format = "SIMPLE_CSV", #下载格式设置为csv
      user = user, #GBIF的账号
      pwd = pwd, #GBIF的密码
      email = email #注册邮箱
    )} else {
  #按地理范围下载POLYGON
  if (substr(range, 1, 7) == "POLYGON") {
    gbif_download <- occ_download(
      pred_default(), #默认过滤规则，见帮助文档
      pred_in("taxonKey", c(keylist$key)), #要下载物种的key
      pred_within(range), #设置下载的地理区域,逆时针连城一个地理多边形
      pred_or(
        pred_lt("coordinateUncertaintyInMeters", uncertainty), pred_isnull("coordinateUncertaintyInMeters"),
        pred_not(pred_in("establishmentMeans", c("MANAGED","INTRODUCED"))), pred_isnull("establishmentMeans")
      ), #排除栽培的或引进的记录
      format = "SIMPLE_CSV", #下载格式设置为csv
      user = user, #GBIF的账号
      pwd = pwd, #GBIF的密码
      email = email #注册邮箱
    )
  }
  #按国家下载
  if (is.character(range) & substr(range, 1, 7) != "POLYGON") {
    if (range == "CN") {
      gbif_download <- occ_download(
        pred_default(), #默认过滤规则，见帮助文档
        pred_in("taxonKey", c(keylist$key)), #要下载物种的key
        pred_or(pred("country", "CN"), pred("country", "TW"),
                pred("country", "HK"), pred("country", "MO")), #国家必须单独放在一个pred_or()内
        pred_or(
          pred_lt("coordinateUncertaintyInMeters", uncertainty), pred_isnull("coordinateUncertaintyInMeters"),
          pred_not(pred_in("establishmentMeans", c("MANAGED","INTRODUCED"))), pred_isnull("establishmentMeans")
        ), #排除栽培的或引进的记录
        format = "SIMPLE_CSV", #下载格式设置为csv
        user = user, #GBIF的账号
        pwd = pwd, #GBIF的密码
        email = email #注册邮箱
      )
    } else {
      gbif_download <- occ_download(
        pred_default(), #默认过滤规则，见帮助文档
        pred_in("taxonKey", c(keylist$key)), #要下载物种的key
        pred("country", range),
        pred_or(
          pred_lt("coordinateUncertaintyInMeters", uncertainty), pred_isnull("coordinateUncertaintyInMeters"),
          pred_not(pred_in("establishmentMeans", c("MANAGED","INTRODUCED"))), pred_isnull("establishmentMeans")
        ), #排除栽培的或引进的记录
        format = "SIMPLE_CSV", #下载格式设置为csv
        user = user, #GBIF的账号
        pwd = pwd, #GBIF的密码
        email = email #注册邮箱
      )
    }
  }}
  #0.3查看gbif是否准备完毕
  rgbif::occ_download_wait(gbif_download)
  #查看引文doi,安息香科数据集'0048799-230530130749713'
  #rgbif::gbif_citation('0048799-230530130749713')
  #0.4下载到电脑并读取到R
  d <- occ_download_get(gbif_download, outdir) %>%
    occ_download_import()
  # table(d$species)
  ####数据过滤####
  #将使用"CoordinateCleaner"包(Zizka et al. 2019)。以帮助检测和删除有问题的记录,本教
  #程将检测记录点坐标是否围绕首都、国家的中心，是否落入海洋，为零，或在饲养动物的博
  #物馆(机构)周围。读者也可以使用世界自然保护联盟的物种分布范围来剔除分布范围以外的所有记录。
  d1 <- d %>%
    cc_cen(buffer = 2000) %>% # remove country centroids within 2km
    cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
    cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km
    cc_sea() %>% # remove from ocean
    cc_val() %>% #emoves or flags non-numeric and not available coordinates
    cc_equ() %>% #删除经纬度相等的坐标
   # cc_coun(iso3 = "countryCode") %>% #删除或标记地理坐标与其他国家信息之间的不匹配
    cc_zero(buffer = 0.5) %>% #删除坐标为0周围的数据
    #cc_outl() %>%
    setNames(tolower(names(.))) %>% #将列名转换为小写方便后面CoordinateCleaner包的处理
    filter(year > 1950) %>% #只选择1950年后的记录
    filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>%
    filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% #已知有错误的坐标不确定度
    #filter(!decimallatitude == 0 | !decimallongitude == 0) %>% #剔除坐标为0的记录
    distinct(decimallongitude, decimallatitude, specieskey, datasetkey, .keep_all = TRUE) #去重

  #从d中提取“物种名称”、“经度”、“纬度”、位置，“国家名称”，并重命名
    datasel <- dplyr::select(d1, species, decimallongitude, decimallatitude, elevation,
                  names(d1)[18], locality, countrycode, names(d1)[30], names(d1)[37],
                  names(d1)[38], names(d1)[39], names(d1)[45]) %>%
    rename(species = species, longitude = decimallongitude, latitude = decimallatitude,
           province = stateprovince)

  if ("TW" %in% unique(datasel$countrycode)) {
    datasel[datasel$countrycode %in% "TW",]$countrycode <- "CN(TW)"
  }
    if ("HK" %in% unique(datasel$countrycode)) {
      datasel[datasel$countrycode %in% "TW",]$countrycode <- "CN(HK)"
    }
    if ("MO" %in% unique(datasel$countrycode)) {
      datasel[datasel$countrycode %in% "TW",]$countrycode <- "CN(MO)"
    }
    datasel <- arrange(datasel, species)

  if (rWCVP == TRUE) {
    cat("Use rWCVP package to filter distribution points outside the country of origin.")
  #使用rWCVP过滤原产地外的分布点
  occs <- datasel %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) #这里是wgs84坐标系

  #按种进行过滤
  species <- unique(occs$species)

  #使用rWCVPdata查看物种产地，这是一个约等于国家级的尺度
  occs_i <- data.frame(matrix(NA, nrow = 0, ncol = ncol(occs)))
  names(occs_i) <- names(occs)

  for (i in seq_along(species)) {
    occs_sp <- dplyr::filter(occs, species == species[i])
    if (nrow(occs_sp) == 0) {next}
    #获取物种原产地
    native_range <- wcvp_distribution(species[i],
                                      taxon_rank = "species",
                                      introduced = FALSE, #是否包括引种地
                                      extinct = FALSE, #是否包括灭绝物种
                                      location_doubtful = FALSE)
    #过滤非原产地的点，从occs中提取位于st_union(native_range)中的点
    suppressMessages(occs$native <- sf::st_intersects(occs,
                                     #native_range,
                                     st_union(native_range),  #合并不同的区域，当具有多个区域时使用
                                     sparse = FALSE)[,1])
    #图层文件的坐标系是基于度的，所以1km大约为0.009度。但随着纬度升高，0.009度会比1km越来越长，
    #这就要求根据数据和分布区位置做相应的调整。让我们观察一下有多少分布记录落在原生区外1km内。

suppressWarnings(suppressMessages(buffered_dist <- native_range %>%
      st_union() %>%
      st_buffer(dist)))   #构建一个缓冲区
    suppressMessages(occs$native_buffer <- st_intersects(occs, buffered_dist, sparse = FALSE)[,1])

    #现在，我们可以丢弃在原生范围之外 >100 公里（大约）的记录。
    occs_filtered <- occs %>% filter(native_buffer)
    occs_i <- rbind(occs_i, occs_filtered)
  }
  occs_i <- arrange(occs_i, species) %>%
    as.data.frame()

  occs_i <- cbind(occs_i[1], st_coordinates(occs_i$geometry), occs_i[2:10], occs_i[12:13])
  names(occs_i)[2:3] <- c("longitude", "latitude")
  }

  if (rWCVP == TRUE) {
    write.csv(occs_i, paste0(outdir, "/", filename, ".csv"),
              row.names = FALSE, fileEncoding = "GB18030", na = "")
    return(occs_i)
  } else {
    write.csv(datasel, paste0(outdir, "/", filename, ".csv"),
              row.names = FALSE, fileEncoding = "GB18030", na = "")
    return(datasel)
  }

}













