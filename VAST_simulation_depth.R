###############################################################################
###############################################################################
## Purpose:    VAST simulation for covariates
## Author:     Kelli Faye Johnson
## Contact:    kelli.johnson@noaa.gov
## Date:       2018-03-26
## Comments:   Adapted from VAST_simulation_2018-03-25.R (JTT)
###############################################################################
###############################################################################

#### Libraries
if( !"FishData" %in% installed.packages()) devtools::install_github("james-thorson/FishData")
if( !"TMBhelper" %in% installed.packages()) devtools::install_github("kaskr/TMB_contrib_R/TMBhelper") # Optimize
if( !"VAST" %in% installed.packages()) devtools::install_github("james-thorson/VAST")

# Load libraries
library(ggplot2)
library(TMB)
library(TMBhelper)
library(VAST)

# Database directory
RootDir = file.path(paste0(letters, ":"), "StockAssessment", "VAST_simulation")
RootDir <- RootDir[file.exists(RootDir)]
source(dir(RootDir, pattern = "functions.R", full.names = TRUE))

# DateDir
Date = Sys.Date()
  DateDir = file.path(RootDir, paste0("VAST_simulation_depth_", Date), "/")
  DownloadDir = paste0(DateDir, "downloads/")
  dir.create(DownloadDir, recursive = TRUE, showWarnings = FALSE)

###############
# Settings
###############
if( file.exists(paste0(DateDir,"Record.RData")) ){
  load( file=paste0(DateDir,"Record.RData"))
  attach(Record)
  message("Loaded Record")
}else{
  # West Coast species (lingcod)
  # Species <- "WCGBTS_Ophiodon_elongatus"
  Species <- "EBSBTS_Pleuronectes_quadrituberculatus"
  Version = "VAST_v2_5_0"
  Method = c("Grid", "Mesh")[2]
  grid_size_km = 25
  n_x = 100  # Number of stations
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  OverdispersionConfig = c("eta1"=0, "eta2"=0)
  Kmeans_Config = list(
    "randomseed" = 1,
    "nstart" = 100,
    "iter.max" = 1e3)     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
  Options = c(
    SD_site_density = 0, SD_site_logdensity = 0,
    Calculate_Range = 1, Calculate_evenness = 0, Calculate_effective_area = 1,
    Calculate_Cov_SE = 0, Calculate_Synchrony = 0, Calculate_Coherence = 0)
  Sim_Settings = list(
    "beta1_mean" = 0, "beta2_mean" = 0,
    "beta1_slope" = 0, "beta2_slope" = 0,
    "beta1_sd" = 0, "beta2_sd" = 0,
    "Nyears" = 10, "Nsamp_per_year" = 600,
    "Depth_km" = 0, "Depth_km2" = 0, "Dist_sqrtkm" = 0,
    "SigmaO1" = 0, "SigmaO2" = 0,
    "SigmaE1" = 0, "SigmaE2" = 0,
    "SigmaV1" = 0, "SigmaV2" = 0,
    "SigmaVY1" = 0, "SigmaVY2" = 0,
    "Range1" = 1000, "Range2" = 500, "SigmaM" = 1,
    "ObsModel" = c(2, 0))

  # Gamma p_1
  # "logit-link" == 0 vs. "log-link" == 1
  ObsModel_Set = list( "Conventional delta"=c(2,0), "Poisson-process link"=c(2,1) )
  # 1=Presence-absence; 2=Density given presence;
  # Epsilon=Spatio-temporal; Omega=Spatial
  FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
  n_rep = 100

  # Decide on case-specific settings for use when calculating indices
  # strata.limits <- data.frame(
  #   "STRATA" = c("Coastwide","CA","OR","WA"),
  #   "north_border" = c(49.0, 42.0, 46.0, 49.0),
  #   "south_border" = c(32.0, 32.0, 42.0, 46.0),
  #   "shallow_border" = c(55, 55, 55, 55),
  #   "deep_border" = c(1280, 1280, 1280, 1280)
  #   )
  strata.limits <- data.frame('STRATA'="All_areas")

  # Save
  Record = list("Species"=Species, "Version"=Version, "Method"=Method, "grid_size_km"=grid_size_km, "n_x"=n_x, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig, "Kmeans_Config"=Kmeans_Config, "Options"=Options, "ObsModel_Set"=ObsModel_Set, "n_rep"=n_rep, "strata.limits"=strata.limits, "Sim_Settings"=Sim_Settings)
  save( Record, file=paste0(DateDir,"Record.RData"))
}

###############
# Run multiple operating models
###############

repI = omI = emI = 1
omdepth = "_nodepth"; observed = "_observed"
for(observed in c("_observed", "_5")) {
for(omdepth in c("", "_nodepth")) {
for(repI in 1:n_rep) {
for(omI in 1:length(ObsModel_Set[1])) {
for(emI in 1:length(ObsModel_Set[1])) {

  # Create directory for OM
  OmDir = paste0(DateDir,
    "OM=", names(ObsModel_Set)[omI], omdepth, observed, "/")
  RepDir = paste0(OmDir,"rep=",repI,"/")
  EmDir = paste0(RepDir,"EM=",names(ObsModel_Set)[emI],"/TRUE/")
    dir.create(EmDir, recursive=TRUE)

  if( file.exists(paste0(DateDir,"DatabaseSave.RData")) ){
    load( file=paste0(DateDir,"DatabaseSave.RData"))
  }else{
    # Surveys with a public API  #   # FishData::
    Database = FishData::download_catch_rates(
      survey = strsplit(Species, "_")[[1]][1],
      species_set = gsub("_", " ", gsub("[A-Z]{3}BTS_", "", Species)),
      error_tol = 0.01, localdir = DownloadDir )
    Database = ThorsonUtilities::rename_columns( Database[,c('Sci','Wt','Year','Long','Lat')], newname=c('Sci','Catch_KG','Year','Lon','Lat') )
    Database = cbind(Database,
      'AreaSwept_km2' = 0.01)  # WCGBTS and all AFSC surveys are in KG/Hectare
    if( !("Vessel" %in% names(Database)) ) Database = cbind( Database, 'Vessel'=1 )
    Database = na.omit( Database )
    # tapply( ifelse(Database$Catch_KG>0,1,0), INDEX=Database$Year, FUN=mean )

    # Save
    DatabaseSave = list("Database"=Database)
    save( DatabaseSave, file=paste0(DateDir,"DatabaseSave.RData"))
  }
  Data_Geostat = DatabaseSave$Database

  # Get extrapolation data
  if(grepl("EBSBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="eastern_bering_sea", strata.limits=strata.limits )
    Extrapolation_List$Data_Extrap = ThorsonUtilities::rename_columns( Extrapolation_List$Data_Extrap, origname="Mid_Depth", newname="Depth_km")
  }
  if(grepl("WCGBTS", Species)){
    Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="California_current", strata.limits=strata.limits )
  }
  # Standardize the depth data
  Extrapolation_List$Data_Extrap[,'Depth_km'] <-
    (Extrapolation_List$Data_Extrap[,'Depth_km'] -
    mean(Extrapolation_List$Data_Extrap[,'Depth_km'])) /
    sd(Extrapolation_List$Data_Extrap[,'Depth_km'])


  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateDir )
  Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

  # Generate covariate
  Depth_x = tapply( Extrapolation_List$Data_Extrap[,'Depth_km'], INDEX=Spatial_List$PolygonList$NN_Extrap$nn.idx, FUN=mean )
  X_xtp = Depth_x %o% rep(1,diff(range(Data_Geostat[,'Year']))+1) %o% 1
  # Change depth in Extrapolation_List so that it matches EM
  Extrapolation_List$Data_Extrap[,'Depth_km'] = Depth_x[Spatial_List$PolygonList$NN_Extrap$nn.idx]

  # Plot settings
  Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
  Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

  ##############################
  # Calculate parameters for OM
  ##############################

  if( file.exists(paste0(OmDir,"Save.RData")) ){
    load( file=paste0(OmDir,"Save.RData") )
  }else{
    # Make TMB data list
    TmbData = VAST::Data_Fn(
      "Version" = Version,
      "X_xtp" = X_xtp,
      "FieldConfig" = FieldConfig,
      "RhoConfig" = RhoConfig,
      "ObsModel" = ObsModel_Set[[omI]],
      "c_i" = rep(0,nrow(Data_Geostat)),
      "b_i" = Data_Geostat[,'Catch_KG'],
      "a_i" = Data_Geostat[,'AreaSwept_km2'],
      "v_i" = as.numeric(Data_Geostat[,'Vessel'])-1,
      "s_i" = Data_Geostat[,'knot_i']-1,
      "t_i" = Data_Geostat[,'Year'],
      "a_xl" = Spatial_List$a_xl,
      "MeshList" = Spatial_List$MeshList,
      "GridList" = Spatial_List$GridList,
      "Method" = Spatial_List$Method,
      "Options" = Options)
    # Make TMB data list
    if (omdepth != ""){
        TmbData = VAST::Data_Fn(
          "Version" = Version,
          "X_xtp" = NULL,
          "FieldConfig" = FieldConfig,
          "RhoConfig" = RhoConfig,
          "ObsModel" = ObsModel_Set[[omI]],
          "c_i" = rep(0,nrow(Data_Geostat)),
          "b_i" = Data_Geostat[,'Catch_KG'],
          "a_i" = Data_Geostat[,'AreaSwept_km2'],
          "v_i" = as.numeric(Data_Geostat[,'Vessel'])-1,
          "s_i" = Data_Geostat[,'knot_i']-1,
          "t_i" = Data_Geostat[,'Year'],
          "a_xl" = Spatial_List$a_xl,
          "MeshList" = Spatial_List$MeshList,
          "GridList" = Spatial_List$GridList,
          "Method" = Spatial_List$Method,
          "Options" = Options)
    }
    # Make TMB object
    TmbList = VAST::Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateDir,
      "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)

    # Run model
    Obj = TmbList[["Obj"]]
    OptOM = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=OmDir, newtonsteps=1, bias.correct = TRUE)
    Report = Obj$report()
    Save = list("Report" = Report, "ParHat" = Obj$env$parList(OptOM$par),
      "OptOM" = OptOM, "Data_raw" = TmbData, "Spatial_raw" = Spatial_List,
      "Version" = Version, "RhoConfig" = RhoConfig)
    save(Save, file=paste0(OmDir,"Save.RData"))

    # Plot index
    Index = SpatialDeltaGLMM::PlotIndex_Fn(DirName=OmDir, TmbData=TmbData,
      Sdreport=OptOM[["SD"]], Year_Set=Year_Set,
      strata_names=strata.limits[,1], use_biascorr=TRUE)
  }

  ##############################
  # Simulate new data
  ##############################

  if( file.exists(paste0(RepDir,"Sim.RData")) ){
    load( file=paste0(RepDir,"Sim.RData") )
  }else{
    # Define parameters
    Sim_Settings[["ObsModel"]] = ObsModel_Set[[omI]]
    Sim_Settings[["beta1_mean"]] = mean(Save$ParHat$beta1_ct)
    Sim_Settings[["beta1_sd"]] = sd(Save$ParHat$beta1_ct)
    Sim_Settings[["beta2_mean"]] = mean(Save$ParHat$beta2_ct)
    Sim_Settings[["beta2_sd"]] = sd(Save$ParHat$beta2_ct)
    Sim_Settings[["SigmaM"]] = exp(Save$ParHat$logSigmaM[1,1])
    Sim_Settings[["Depth1_km"]] = mean(Save$ParHat$gamma1_ctp)
    Sim_Settings[["Depth2_km"]] = mean(Save$ParHat$gamma2_ctp)

    # Define spatial and spatio-temporal variation
    # I'm using lower values than observed so that its less likely to have replicates with 0% or 100% encounter rates
    Sim_Settings[["SigmaO1"]] = 0.5 # abs(Save$ParHat$L_omega1_z)
    Sim_Settings[["SigmaO2"]] = 0.5 # abs(Save$ParHat$L_omega2_z)
    Sim_Settings[["SigmaE1"]] = 0.5 # abs(Save$ParHat$L_epsilon1_z)
    Sim_Settings[["SigmaE2"]] = 0.5 # abs(Save$ParHat$L_epsilon1_z)

    if (observed == "_observed") {
      Sim_Settings[["SigmaO1"]] = abs(Save$ParHat$L_omega1_z)
      Sim_Settings[["SigmaO2"]] = abs(Save$ParHat$L_omega2_z)
      Sim_Settings[["SigmaE1"]] = abs(Save$ParHat$L_epsilon1_z)
      Sim_Settings[["SigmaE2"]] = abs(Save$ParHat$L_epsilon1_z)
    }

    # Check encounter rate
    Prop_t = 0
    counter <- 0
    while(any(any(Prop_t==0) | any(Prop_t==1)) & counter < 10) {
      counter <- counter + 1
      Sim = SpatialDeltaGLMM::Geostat_Sim(Sim_Settings=Sim_Settings, Extrapolation_List=Extrapolation_List, Data_Geostat=Data_Geostat )
      Prop_t = tapply( Sim$Data_Geostat[,'Catch_KG'], INDEX=Sim$Data_Geostat[,'Year'], FUN=function(vec){mean(vec>0)} )
    }
    save(Sim, file=paste0(RepDir,"Sim.RData"))
    # Plot
    ThorsonUtilities::save_fig(file=paste0(RepDir,"Index-Sim"),
      width=5, height=5, res=200, units='in')
    plot( x=Year_Set, y=Sim$B_tl[,1]/1000, type="b",
      ylim=c(0,max(Sim$B_tl[,1]/1000)) )
    dev.off()
  }

  ##############################
  # Fit using different models
  ##############################

  # Run model if necessary
  if( !file.exists(paste0(EmDir,"parameter_estimates.RData")) ){

    mats <- Sim$B_tl[,1]/1000
    # Compare true and estimated spatial and spatio-temporal SD
    sds <- data.frame("OM"=unlist(Save$ParHat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')]))
    # Compare true and estimated covariate effects
    covars <- data.frame("OM"=c(mean(Save$ParHat[['gamma1_ctp']]),mean(Save$ParHat[['gamma2_ctp']])))

    for(emCov in c(TRUE, FALSE)) {
    # Make TMB data list
    if (emCov) {
      TmbDataEM = Data_Fn(
        "Version"=Version,
        "X_xtp"=X_xtp,
        "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
        "ObsModel"=ObsModel_Set[[emI]],
        "c_i"=rep(0,nrow(Sim$Data_Geostat)),
        "b_i"=Sim$Data_Geostat[,'Catch_KG'],
        "a_i"=Sim$Data_Geostat[,'AreaSwept_km2'],
        "v_i"=as.numeric(Sim$Data_Geostat[,'Vessel'])-1,
        "s_i"=Spatial_List$knot_i-1,
        "t_i"=Sim$Data_Geostat[,'Year'],
        "a_xl"=Spatial_List$a_xl,
        "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList,
        "Method"=Spatial_List$Method, "Options"=Options)
    }
    if (!emCov) {
      TmbDataEM <- Data_Fn(
        "Version"=Version,
        "X_xtp"= NULL,
        "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig,
        "ObsModel"=ObsModel_Set[[emI]],
        "c_i"=rep(0,nrow(Sim$Data_Geostat)),
        "b_i"=Sim$Data_Geostat[,'Catch_KG'],
        "a_i"=Sim$Data_Geostat[,'AreaSwept_km2'],
        "v_i"=as.numeric(Sim$Data_Geostat[,'Vessel'])-1,
        "s_i"=Spatial_List$knot_i-1,
        "t_i"=Sim$Data_Geostat[,'Year'],
        "a_xl"=Spatial_List$a_xl,
        "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList,
        "Method"=Spatial_List$Method, "Options"=Options )
        EmDir <- gsub("TRUE", "FALSE", EmDir)
        dir.create(EmDir, recursive = TRUE)
    }

    # Make TMB object
    # TmbDir <- c("D:/spatial/VAST/inst/executables",
    #   "D:/StockAssessment/VAST/inst/executables")
    # TmbDir <- TmbDir[file.exists(TmbDir)]
    counter <- 0
    OptEM <- list("SD" = NULL)
    while (is.null(OptEM[["SD"]]) & counter < 2) {
      TmbListEM  = Build_TMB_Fn("TmbData"=TmbDataEM, "RunDir"=DateDir,
        "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)
      # Run model
      ObjEM = TmbListEM[["Obj"]]
      OptEM = try(TMBhelper::Optimize(obj=ObjEM,
        lower=TmbListEM[["Lower"]], upper=TmbListEM[["Upper"]],
        getsd=TRUE, savedir=EmDir, newtonsteps=1, bias.correct=TRUE),
        silent = TRUE)
      if (class(OptEM) == "try-error") OptEM <- list("SD" = NULL)
      counter <- ifelse(is.null(OptEM[["SD"]]), counter + 1, 3)
    }
    ReportEM = ObjEM$report()

    # Plot index
    if (!is.null(OptEM[["SD"]])) {
      IndexEM = SpatialDeltaGLMM::PlotIndex_Fn(DirName=EmDir,
        TmbData=TmbDataEM, Sdreport=OptEM[["SD"]],
        Year_Set=Year_Set, strata_names=strata.limits[,1], use_biascorr=TRUE)
      # save parameters
      ParHat = ObjEM$env$parList( OptEM$par )
    } else {
      IndexEM$Index_ctl <- array(NA, dim = dim(Index$Index_ctl),
        dimnames = list(NULL, NULL, NULL, c("Estimate", "Std. Error")))
      ParHat <- list('L_omega1_z' = NA,'L_epsilon1_z' = NA,'L_omega2_z' = NA,'L_epsilon2_z' = NA, 'gamma1_ctp' = NA, 'gamma2_ctp' = NA)
    }
    EmSave = list("Report"=ReportEM, "ParHat"=ParHat, "Data_EM" = TmbDataEM,
      "Sdreport" = OptEM[["SD"]])
    save(EmSave, file=paste0(EmDir,"EmSave.RData"))
    mats <- cbind(mats, IndexEM$Index_ctl[,,1,'Estimate'])
    sds <- cbind(sds, "EM"=unlist(EmSave$ParHat[c('L_omega1_z','L_epsilon1_z','L_omega2_z','L_epsilon2_z')]))
    covars <- cbind(covars, "EM"=c(mean(EmSave$ParHat[['gamma1_ctp']]),
      mean(EmSave$ParHat[['gamma2_ctp']])))

    # Unlink
    dyn.unload(paste0(DateDir,"/",TMB::dynlib(Version)))
  }

    write.table(sds,
      file.path(gsub("FALSE|TRUE", "", EmDir), "sds.csv"),
      sep = ",")
    write.table(covars,
      file.path(gsub("FALSE|TRUE", "", EmDir), "covars.csv"),
      sep = ",")
    jpeg(file.path(gsub("FALSE|TRUE", "", EmDir), "matplot.jpeg"))
    matplot(mats, col=1:NCOL(mats), lty="solid", type="l",
      ylab = "index", xlab = "year")
    legend("topright", bty = "n", lty = "solid", col = 1:NCOL(mats),
      legend = c("OM", "EM_covariate", "EM_0")[1:NCOL(mats)])
    dev.off()
    rm(covars, mats, sds)
  }
}}}}} # omI

keep <- get_results(gsub("/$", "", DateDir))
keep_re <- re_all(data = keep)
keep_re <- subset(keep_re, gradient < 15)
dim(keep);dim(keep_re)

pdf(file.path(DateDir, "VAST_simulation_depth_plots.pdf"))
for (ii_name in grep("re_", colnames(keep_re), value = TRUE)) {
  for (iii_name in grep("re_", colnames(keep_re), value = TRUE)) {
    if (ii_name == iii_name) next
    if (grepl("index", ii_name) | grepl("index", iii_name)) next
    ggplotre(keep_re, ii_name, iii_name, type = "points")
  }
}
dev.off()

#### Results

## Tables
rawests <- Reduce(function(x,y) merge(x, y, all = TRUE),
  lapply(dir(dir(DateDir, full.names = TRUE, pattern = "observed"),
    recursive = FALSE, "parameter_estimates.RData", full.names = TRUE),
  getrdata))
rawests <- reshape(rawests,
  direction = "wide", idvar = "par", timevar = "om_name")
colnames(rawests) <- gsub(" delta_observed", "", colnames(rawests))

## Plots
results <- keep_re[grepl("Conventional", keep_re$om_type) &
  grepl("Conventional", keep_re$em_type), ]
gg <- ggplotre(results[results$em_depth == TRUE, ],
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "mean depth effect for the first model component",
    "mean depth effect for the second model component"),
  type = "points", facetx = c("om_sigmao1"), facety = c("om_depth"),
  dir = NULL, scales = "fixed")
gg <- gg  +
  geom_vline(aes(xintercept = om_depth1_km),
    col = "red") +
  geom_hline(aes(yintercept = om_depth2_km),
    col = "red")
ggsave(file.path(DateDir,
  "VAST_simulation_depth_em_depth1kmVSem_depth2_km.jpeg"),
  gg, dev = "jpeg")


gg <- ggplotre(results,
  "em_depth", "re_logratiodiff",
  labels = c(
    "", "error in logged ratio of the first and last year"),
  type = "box", facetx = c("om_depth"), facety = c("."),
  dir = DateDir, scales = "fixed")

gg <- ggplotre(results,
  # [results$em_depth == FALSE & results$om_depth == TRUE, ]
  "re_range1", "re_sigmao1",
  labels = c(
    "re in range the spatial field",
    "re in marginal variance of the spatial field for the first component"),
  type = "points", facetx = c("em_depth"), facety = c("om_depth"),
  dir = NULL, scales = "fixed", print = FALSE, gradient = TRUE)
ggg <- ggplotre(
  results[results$em_depth == FALSE & results$om_depth == TRUE, ],
  "re_range1", "re_sigmao1",
  labels = c("", ""),
  type = "points", facetx = c("em_depth"), facety = c("om_depth"),
  dir = NULL, scales = "fixed", print = FALSE, gradient = TRUE)
ggg <- ggg + xlim(c(-0.9, -0.25)) + ylim(c(0.15, 0.8)) + ylab("") + xlab("")+
  scale_colour_gradient(guide = FALSE)
vp <- grid::viewport(width = 0.35, height = 0.35,
  x = 0.26, y = 0.312)
jpeg(file.path(DateDir, "VAST_simulation_depth_re_range1VSre_sigmao1.jpeg"),
  res = 600, units = "in", height = 8, width = 8)
  print(gg)
  theme_set(theme_bw(base_size = 8))
  print(ggg, vp = vp)
dev.off()

#### Appendix A
matched <- keep_re[keep_re$om_depth == TRUE & keep_re$om_sigmao1 == 0.5, ]
matched <- matched[ifelse(matched$om_type != matched$em_type, FALSE, TRUE), ]

# Error in log ratio of last and first years
gg <- ggplotre(matched, "re_logratiodiff", "re_sigmao1",
  labels = c(
    "error in logged ratio of the first and last year",
    "RE in marginal variance of the spatial field for the first component"),
  type = "points",
  facety = c("om_type"), facetx = c("em_depth"),
  dir = NULL, gradient = TRUE, scales = "free")
ggsave(file.path(DateDir,
  "VAST_simulation_depth_AppendixA_re_sigmao1VSre_logratiodiff.jpeg"),
  gg, dev = "jpeg")

gg <- ggplotre(matched[matched$em_depth == TRUE, ],
  "em_depth1_km", "em_depth2_km",
  labels = c(
    "mean depth effect for the first model component",
    "mean depth effect for the second model component"),
  type = "points", facetx = c(""), facety = c(""),
  dir = NULL)
gg <- gg + geom_line(data = data.frame(
    x = c(-0.25, 0.5), y = as.numeric(unique(matched$om_depth2_km)[2])),
    aes(x = x, y = y), col = 2) +
  geom_line(data = data.frame(
    x = c(1.0, 1.4), y = as.numeric(unique(matched$om_depth2_km)[1])),
    aes(x = x, y = y), col = 1) +
  geom_vline(aes(xintercept = as.numeric(om_depth1_km),
    col = as.factor(om_type))) +
  scale_colour_manual(guide = FALSE, values = 1:2)
ggsave(file.path(DateDir,
  "VAST_simulation_depth_AppendixA_em_depth1_kmVSem_depth2_km.jpeg"),
  gg, dev = "jpeg")
dev.off()


# Linear trend in log-abundance
# We next compare the trend in true abundance (fitting a log-linear model to true abundance) to the trend in estimated abundance (Fig. 2).  This again shows that the model is unbiased in each scenario.
# caption -- Error in estimates of trend for scenarios with stable, increasing, or decreasing trends in abundance, where blue line is true expected trend, and the number in the top-left is the average estimated trend for each scenario (true expected trend is 0, -0.1, and 0.1 for the 1st, 2nd, and 3rd panels, respectively)
# with(keep, plot(as.numeric(em_linear) - as.numeric(om_linear),
#   gradient))

# Check density-dependent catchability
# Finally, we calculate density-dependent catchability (whether estimated indices are hyperstable or hypersensitive to changes in the true index).  This shows a delta of approx. 1.0 for every scenario, as also shown in Thorson et al. 2015 ICESJMS.
# caption -- Test for hypersensitive or hyperstable indices (Estimate = 1.0 implies well-calibrated sensitivity)


# TODO: all the below needs to be altered.
# DF = data.frame("True"=as.vector(Results_index[,,"True","Value",]), "Est"=as.vector(Results_index[,,"Est","Value",]), "Scenario_"=dimnames(Results_index)[[1]], "RepNum"=rep(1:dim(Results_index)[2],each=dim(Results_index)[1]), "Year"=rep(1:dim(Results_index)[5],each=prod(dim(Results_index)[1:2])))
# Lm = lme4::lmer(log(Est) ~ 0 + Scenario_ + log(True):Scenario_ + (1 | RepNum:Scenario_), data=DF)
# Coef = summary(Lm)$coef
# Coef[grep("True", rownames(Coef)), ]
