wendang(

spdir = c("wz001", "wz002", "wz003"),
evdir = "F:/var2/eblf_proj/present/asc",
myenv = NULL,
evlist = NULL,
factors = c("human_lulc", "soil_class", "soil_drainage", "soil_root_depth"),
mybgfile = 1,
nbg = 10000,
args = maxent_args(),
fc = c("lq"),
rm = seq(0.5,5,0.5),
r = 0.7,
cormethod = "pearson",
vif = T,
vifth = 10,
opt = "aicc",
prodir = list(acc2030ssp126 = "/home/215why/sdm/env_asc/acc2030ssp126"),
outdir = "/home/215why/sdm/3"
)
path <- "C:/Users/why/TBlabENM/inst/extdata"
wendang <- function(spdir, evdir, myenv, evlist,
         factors, mybgfile, nbg, args,
         fc, rm, r, cormethod, vif,
         vifth, opt, prodir, outdir) {

  rmarkdown::render(paste0(path, "/Models detail.Rmd"),
                    output_file = "F:/2/report_body.html",
                    params = list(spdir = spdir,
                                  evdir = evdir,
                                  myenv = myenv,
                                  evlist = evlist,
                                  factors = factors,
                                  mybgfile = mybgfile,
                                  nbg = nbg,
                                  args = args,
                                  fc = fc,
                                  rm = rm,
                                  r = r,
                                  cormethod = cormethod,
                                  vif = vif,
                                  vifth = vifth,
                                  opt = opt,
                                  prodir = prodir,
                                  outdir = outdir))

}


rmarkdown::render("C:/Users/why/TBlabENM/inst/extdata/1.Rmd", output_file = "F:/2/1.html",
                  params = list(spdir = "Q",
                                evdir = "M"))
fc <- toupper(fc)
fc <- c("q", "l")
combin <- expand.grid(fc, rm, stringsAsFactors = FALSE) %>%
  dplyr::mutate(fc = purrr::map_chr(.x = Var1, .f = function(x){
    paste0(tzhs[tzhs %in% stringr::str_split_1(x, "")], collapse = "")}
  ))
