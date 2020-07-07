# Model steps:
# 1. Prep to one single bd raster (where math takes place)
# 2. Write and reload
# 3. load raster as a data.table
# 4. assign name "cons.priorities"
# 5. fwrite data.table as a csv, load back in as a csv.
# 6. set production targets
# 7. set input folder for specific run.

# prepping for tradeoff_mod run
# 1. do raster combination.
# 2. write_reload raster


# last updated 0ct 7, 2019. See "cc_functions.R"
cc_tradeoff_mod_outputs <- function(mod_toff, input_key, cbetas) {

  # x y coordinates to reconstruct rasters.
  cp <- fread(paste0(p_dat, .Platform$file.sep,
           input_key, .Platform$file.sep,
           input_key, "-cons-priorities.csv"))
  msk_csv <- fread(paste0(p_dat, .Platform$file.sep,
                          input_key, .Platform$file.sep,
                          input_key, "-mask.csv"))

  bd_input <- cp %>%
    cbind(msk_csv, .) %>%
    dt_to_raster(., CRSobj) %>%
    .$cons.priorities

  zambia_xy <- mod_toff$inputs$mask[, .(x, y)]
  CRS_obj <- mod_toff$inputs$sp$crs

  # standard model outputs
  cons_priorities_r <- mod_toff$inputs$cons$cons.priorities %>%
    cbind(zambia_xy, .) %>%
    dt_to_raster(., CRSobj)


  conv_r <- dt_to_raster(mod_toff$conv, CRSobj) # raster brick of areas that get converted, to maize and to soy
  # mod_conv_r has two layers:
  # mod_conv_r$maize
  # mod_conv_r$soy
  conv_r_ms <- calc(dt_to_raster(mod_toff$conv, CRSobj) * c(1, 2), sum) # produces conversion raster with value of 1 for maize, 2 for soy
  conv_r_all <- calc(conv_r, sum) # one raster, with 1 for conversion

  # for visualizations:
  # mod_conv_r_maize[mod_conv_r_maize == 0] <- NA
  # mod_conv_r_soy[mod_conv_r_soy == 0] <- NA # to show only areas converted to soy

  ## create conversion probabilities raster
  convprob_r <- mod_toff %>%
    list(Y = .$inputs$y_std,
         C = .$inputs$carbon_p,
         BD = .$inputs$cons_p,
         COST = .$inputs$cost_p) %>%
    constraints(inlist = ., cbetas = cbetas, silent = FALSE) %>%
    cbind(zambia_xy, .) %>%
    dt_to_raster(., CRSobj)

  ## just converted maize areas gradient
  conv_grad_maize <- conv_r$maize * convprob_r$maize
  conv_grad_soy <- conv_r$soy * convprob_r$soy

  list(
    #name = object,
    bd_input = bd_input,
    bd_input_no_pas = cons_priorities_r,
    conv_r = conv_r,
    conv_r_ms = conv_r_ms,
    conv_r_all = conv_r_all,
    convprob_r = convprob_r,
    conv_grad_maize = conv_grad_maize,
    conv_grad_soy = conv_grad_soy
  )
}

