
get_results <- function(dir, truncate = NULL) {
  dir.reps <- dir(dir, pattern = "rep=", full.names = TRUE,
    recursive = TRUE, include.dirs = TRUE)
  if (is.null(truncate)) truncate <- seq_along(dir.reps)

  get_OM <- function(file) {
    if (!file.exists(file)) return(NULL)
    load(file)
    index <- Sim$B_tl[,1]/1000
    names(index) <- paste0("om_index_", seq_along(index))
    linear <- lm(log(index) ~ 1 + I(seq_along(index)))$coef
    names(linear) <- NULL
    sets <- Sim$Sim_Settings[c("Range1", "Range2", "SigmaO1",
      "SigmaO2", "SigmaE1", "SigmaE2", "Depth1_km", "Depth2_km")]
    names(sets) <- paste0("om_", tolower(names(sets)))
    om <- unlist(c(index, sets))
    return(c("om_type" = gsub(".+OM=|/.+$", "", file), om,
      "om_linear" = linear[2],
      "n_t" = length(index),
      "om_depth" = ifelse(grepl("nodepth", file), FALSE, TRUE)))
  }
  get_EM <- function(file) {
    if (!file.exists(file)) return(NULL)
    load(file)
    em <- unlist(c(
      setNames(EmSave$Report[c("Range_raw1", "Range_raw2")],
        c("em_range1", "em_range2")),
      setNames(EmSave$ParHat[
        c("L_omega1_z", "L_omega2_z", "L_epsilon1_z", "L_epsilon2_z")],
        c("em_sigmao1", "em_sigmao2", "em_sigmae1", "em_sigmae2")),
      setNames(
        c(mean(EmSave$ParHat[['gamma1_ctp']]),
          mean(EmSave$ParHat[['gamma2_ctp']])),
        c("em_depth1_km", "em_depth2_km"))))
    em[grep("sigma", names(em))] <- abs(em[grep("sigma", names(em))])
    em_type <- gsub(".+rep=[0-9]+/EM=|/EmSave.+$", "", file)
    if (length(dir(dirname(file),
      full.names = TRUE, pattern = "Table_for_SS3")) == 0) {
      index <- rep(NA, times = EmSave$Data_EM$n_t)
      linear <- rep(NA, 2)
    } else {
      index <- read.table(dir(dirname(file),
        full.names = TRUE, pattern = "Table_for_SS3"), sep = ",", header = TRUE)
      index <- index[, "Estimate_metric_tons"]
      linear <- lm(log(index) ~ 1 + I(seq_along(index)))$coef
      names(linear) <- NULL
    }
    names(index) <- paste0("em_index_", seq_along(index))
    return(c("em_type" = gsub("TRUE|FALSE|/", "", em_type),
      "em_depth" = grepl("TRUE", em_type),
      em, index, "em_linear" = linear[2],
      "gradient" = ifelse(is.null(EmSave$Sdreport), NA, max(abs(EmSave$Sdreport$gradient.fixed))),
      "hessian" = ifelse(is.null(EmSave$Sdreport), NA, EmSave$Sdreport$pdHess))) #<1e-6
  }
  info <- lapply(dir.reps[truncate], function(x) {
    omom <- get_OM(file.path(x, "Sim.RData"))
    if (is.null(omom[1])) return(NULL)
    dir.ems <- dir(x, pattern = "EmSave", recursive = TRUE,
      full.names = TRUE, include.dirs = TRUE)
    emem <- t(mapply(get_EM, dir.ems))
    rownames(emem) <- NULL
    stuff <- data.frame(t(omom), emem, stringsAsFactors = FALSE)
    stuff$rep <- gsub(".+rep=", "", x)
    return(stuff)
  })
  info <- do.call("rbind", info)

  numbers <- apply(info, 2,
    function(x) stringr::str_extract(x, "\\-*\\d+\\.*\\d*"))
  info <- cbind(
    info[, apply(numbers, 2, function(x) all(is.na(x)))],
    apply(numbers[, !apply(numbers, 2, function(x) all(is.na(x)))], 2,
      as.numeric))

  info$om_logratio <-
    log(as.numeric(info[, paste0("om_index_", info[1, "n_t"])])) -
    log(as.numeric(info[, paste0("om_index_1")]))
  info$em_logratio <-
    log(as.numeric(info[, paste0("em_index_", info[1, "n_t"])])) -
    log(as.numeric(info[, paste0("em_index_1")]))

  return(info)
}

re_all <- function(data, bind = TRUE) {
  re <- function(text = "index_1", data) {
    omd <- as.numeric(data[, paste0("om_", text)])
    emd <- as.numeric(data[, paste0("em_", text)])

    return((emd - omd) / omd)
  }
  names <- gsub("om_", "", grep("om_", colnames(data), value = TRUE))
  names <- names[!(names %in% c("type"))]
  names <- names[!(names %in% c("logratio"))]
  names <- names[!grepl("depth$|linear|logratio", names)]

  res <- mapply(re, names, MoreArgs = list(data = data))
  colnames(res) <- paste0("re_", colnames(res))
  res <- data.frame(res,
    "re_logratiodiff" = data$om_logratio - data$em_logratio)

  if (bind) {
    return(cbind(data, res))
  } else {return(res)}

}

ggplotre <- function(data, x, y, print = TRUE, gradient = FALSE,
  facetx = c("emname", "em_depth"), facety = "omname",
  labels = NULL, type = "box", scales = c("fixed", "free"),
  dir = NULL) {

  scales <- match.arg(scales)
  facetx <- facetx[!facetx %in% x]

  data$em_depth <- paste("depth in EM =", data$em_depth)

  data$em_type <- gsub("process ", "", data$em_type)
  data$om_type <- gsub("process ", "", data$om_type)
  data$em_type <- gsub("Conventional delta", "Delta-model", data$em_type)
  data$om_type <- gsub("Conventional delta", "Delta-model", data$om_type)
  data$omname <- paste("OM =", data$om_type)
  data$emname <- paste("EM =", data$em_type)

  if (x != "em_depth") data[, x] <- as.numeric(data[, x])
  data[, y] <- as.numeric(data[, y])

  gg <- ggplot(data)
  if (gradient) gg <- ggplot(data,
    aes(col = as.numeric(gradient)))

  if (type == "points"){
    gg <- gg +
      geom_point(aes_string(x = x, y = y),
      alpha = 0.5, size = 2.5) +
      scale_shape_manual(name = "", values = c(15, 19), guide = FALSE) +
      scale_color_gradient(name = "gradient")
    if (grepl("^re", y)) {
      gg <- gg + geom_vline(xintercept = 0, lty = 2, col = "red")
    }
  }
  if (type == "box") {
    gg <- gg +
      geom_boxplot(aes_string(x = x, y = y))
  }
  if (facetx == "" & facety == "") {
    facet <- "."
  } else {
    gg <- gg + facet_grid(as.formula(paste(paste(facety, collapse = "+"), "~",
      paste(facetx, collapse = "+"))), scales = scales)
  }

  gg <- gg +
    theme_bw() + theme(
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1))
  if (grepl("^re", y)) {
    gg <- gg + geom_hline(yintercept = 0, lty = 2, col = "red")
  }


  if(!is.null(labels[1])) gg <- gg +
    xlab(labels[1]) + ylab(labels[2])

  if (print) print(gg)
  if (!is.null(dir)) {
    jpeg(
      filename = file.path(dir,
        paste0("VAST_simulation_depth_", x, "VS", y, ".jpeg")),
      res = 600, units = "in", width = 8, height = 8)
    print(gg)
    dev.off()
  }
  invisible(gg)
}

getrdata <- function(file){
  ne <- new.env()
  load(file, env = ne)

  out <- data.frame(
    "par" = names(
      get("parameter_estimates", env = ne)$SD[c("par.fixed")][[1]]),
    "val" = get("parameter_estimates", env = ne)$SD[c("par.fixed")],
    "se" = sqrt(diag(
      get("parameter_estimates", env = ne)$SD[c("cov.fixed")][[1]])))
  colnames(out)[2] <- "val"
  report <- data.frame(
    "par" = names(get("parameter_estimates", env = ne)$SD$value),
    "val" = get("parameter_estimates", env = ne)$SD$value,
    "se" = sqrt(diag(get("parameter_estimates", env = ne)$SD$cov)))
out$par <- make.unique(as.character(out$par))
report$par <- make.unique(as.character(report$par))

  # ne2 <- new.env()
  # load(file.path(dirname(file), "Save.RData"), env = ne2)
  # depth <- data.frame(
  #   "par" = c("gamma1_ctp", "gamma2_ctp"),
  #   "val" = sapply(get("Save",
  #     env = ne2)$ParHat[c("gamma1_ctp", "gamma2_ctp")], mean),
  #   "se" = NA)
  # report <- rbind(report, depth)

  out <- rbind(out, report)
  out$om_name <- basename(dirname(file))
  return(out)
}
