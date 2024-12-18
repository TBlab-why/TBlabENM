% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ENMrefuges.R
\name{ENMrefuges}
\alias{ENMrefuges}
\title{Individual Species Refuge}
\usage{
ENMrefuges(
  parameters,
  x = NULL,
  booleandir,
  key = NULL,
  overwrite = FALSE,
  outdir = NULL,
  parallel = F,
  ncpu = 2
)
}
\arguments{
\item{parameters}{数据框，第一列为物种名(要与resultdir下的名称保持一致)}

\item{x}{数值向量,表示物种序号,与parameters中的物种序号一致。}

\item{booleandir}{适生区和非适生区结果文件路径。}

\item{key}{要提取的关键词，一般对应投影时期}

\item{overwrite}{logical. If TRUE, filename is overwritten}

\item{outdir}{输出文件夹}

\item{parallel}{是否并行}

\item{ncpu}{并行时cpu个数
本文件夹，请自己指定。}
}
\value{
返回表示避难所区域的Boolean图
}
\description{
根据适生区/非适生区栅格计算物种的避难所。
}
\examples{
ENMrefuges(parameters = read.csv("F:/eblf/TBlabENM/spdata.csv"),
           x = c(3,5,7,9),
           booleandir = "F:/eblf/TBlabENM/BooleanmapT_MTSStr",
           key = list("future"= c("Present", "2050RCP8.5", "2070RCP8.5"),
                      "past"= c("Present", "Mid-Holocene", "Last-Glacial-Maximum")),
           overwrite = T,
           outdir = "F:/4",
           parallel = T,
           ncpu = 2)
}
