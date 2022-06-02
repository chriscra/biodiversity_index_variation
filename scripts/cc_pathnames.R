## Pathnames --------------------------------------------------------------------------------

# Zambia specific file path names
p_proj <- "/Volumes/GoogleDrive/My Drive/Zambia/agroEcoTradeoff" # new, 4/27/2022
# p_proj <- "/Users/christophercrawford/Google_Drive/_Projects/Zambia/agroEcoTradeoff"
p_dat <-     fp(p_proj,"external/data") # this has the "input_key" folders: ZA
p_ZA <-      fp(p_dat,"ZA")
p_datnew <-  fp(p_proj,"external/data_new")


p_basedat <- fp(p_proj,"external/base_data")
p_dev <-     fp(p_proj,"external/input_devel")
p_bd <-      fp(p_proj,"external/input_devel/biodiversity")
p_output <-  fp(p_proj,"external/output")
p_mod_output <-  fp(p_proj,"external/mod_output")
p_final_inputs <-  fp(p_proj,"external/data_new/final_inputs")
p_mod_inputs <- fp(p_proj, "external/data_new/final_mod_input_rasters")

# for General data sources
p_dat_main <- "/Volumes/GoogleDrive/My Drive/data" # new, 4/27/2022
# p_dat_main <- "/Users/christophercrawford/Google_Drive/_Projects/data"

p_dat_ag <-  fp(p_dat_main,"Ag")
p_dat_bd <-  fp(p_dat_main,"Bd")
p_dat_clim <-  fp(p_dat_main,"Climate")
p_dat_geopol <- fp(p_dat_main,"Geopolitical")
p_dat_hab <- fp(p_dat_main,"Habitats")
p_dat_lu <- fp(p_dat_main,"LU")
p_dat_socio <- fp(p_dat_main,"Socio")
p_dat_water <- fp(p_dat_main,"Water")
p_mam <-  fp(p_dat_bd,"IUCN_RedList/MAMMALS")
p_rep <-  fp(p_dat_bd,"IUCN_RedList/REPTILES")
p_amp <-  fp(p_dat_bd,"IUCN_RedList/AMPHIBIANS")
p_bird <- fp(p_dat_bd,"IUCN_RedList/Birds")

# for data layers created in this notebook
p_iucn_dev <-  fp(p_datnew,"1_IUCN_dev")
p_plants_dev <-  fp(p_datnew,"2_plants_dev")
p_hotspot_dev <-  fp(p_datnew,"4_bd_hotspots_dev")
p_IBA_dev <-  fp(p_datnew,"5_IBA_EBA_dev")
p_ecoreg_dev <-  fp(p_datnew,"6_ecoregions_dev")
p_wilder_dev <-  fp(p_datnew,"7_wilderness_dev")
p_basemaps <- fp(p_datnew,"basemaps")

p_temp <- fp(p_datnew, "temp")

p_plots <- "/Volumes/GoogleDrive/My Drive/Zambia/agroEcoTradeoff/external/plots" # new, 4/27/2022
# p_plots <- "/Users/christophercrawford/Google_Drive/_Projects/Zambia/agroEcoTradeoff/external/plots"
