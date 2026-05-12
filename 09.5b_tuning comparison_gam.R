### comparing model outputs

### this is a hordcoded comparison between the model AIC's etc. 
##. between the first established models on the wholecovariate set and then on the set established in tuning




files_new<- list.files(here::here("outputs","Sumava","models_gam"),
                    pattern = "^model_gam_.*\\.rds$",
                    full.names = TRUE)

files_old<- list.files(here::here("outputs","Sumava","old_untuned"),
                    pattern = "^model_gam_.*\\.rds$",
                    full.names = TRUE)

get_model_info <- function(path) {
  m <- readRDS(path)
  data.frame(
    file = basename(path),
    response = as.character(stats::formula(m))[2],
    formula = paste(deparse(stats::formula(m)), collapse=" "),
    AIC = AIC(m),
    dev_explained = tryCatch(summary(m)$dev.expl, error=function(e) NA_real_),
    r_sq = tryCatch(summary(m)$r.sq, error=function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}


model_info_new <- do.call(rbind, lapply(files_new, get_model_info))


model_info_old <- do.call(rbind, lapply(files_old, get_model_info))
