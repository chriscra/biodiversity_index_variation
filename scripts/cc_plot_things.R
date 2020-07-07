# cc_plot_things
load(file = fp(p_ZA,"parks_roads.rda"), verbose = TRUE) # includes roads (a SpatialLinesDataFrame),
# pas (SpatialPolygonsDataFrame, a shapefile that includes both national parks and GMAs), and
# zambia (SpatialPolygonsDataFrame, outline of Zambia)
msk_shp <- readOGR(fp(p_datnew,"msk.shp")) %>%
  spTransform(CRS("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
rm(zambia, roads)

# plot colors

pacols <- c("grey80", "grey70") # colors for NPs and GMAs
col_pas <- c("grey80", "grey70")
col_pas_all <- col_pas[pas$type]

col_overlap <- c(
  "grey90", # grey # no conversion
#  "#1F78B4", # blue # Mod 1 only
  "#FF7F00", # orange # Mod 1 only
  "#FDBF6F", # yellow # Mod 2 only
#  "#E31A1C", # red # Mod 2 only
  "#33A02C" # green # Both Models
) # light purple # both soy

#display.brewer.all(n=9, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=FALSE)

## main color scheme for plots (0 = green / high bd, 1 = pink / convert)
# col_main <- colorRampPalette(brewer.pal(n=11,name='PiYG'))(100) #pink-yellow-green
# col_main2 <- col_main[15:95] # slightly muted on both ends, for s4,5,6,7
# col_main3 <- col_main[16:100] # toned down pink, for s8
# col_main_legend <- c(col_main[15],col_main[85])
#
# col_grad0 <- colorRampPalette(brewer.pal(n=9,name='YlOrRd'))(100) # for graduated plot of tradeoff_mod results rasters
# col_grad <- col_grad0[1:75] # more mellow red at the top end
# col_grad1 <- col_grad0[30:100] # darker yellow at the bottom
# col_grad2 <- col_grad
# col_grad2[76] <- "#662506"



#col_overlap <- rev(brewer.pal(n=9,name='RdYlGn')) # for diverging plot of overlapping results rasters
#col_overlap[1] <- "grey90" # removing the dark green color from the front

#col_overlap <- brewer.pal(n=9,name='Paired')
# col_overlap <- c(
#   "grey90", # grey # no conversion
#   "#1F78B4", # blue # estes maize only
#   "#B2DF8A", # light green # estes soy only
#   "#E31A1C", # red # laurance maize only
#   "#33A02C", # green # both maize
#   "#FF7F00", # orange # L maize, E soy
#   "#FDBF6F", # yellow # laurance soy only
#   "#FB9A99", # pink # Laurance soy, Estes maize
#   "#CAB2D6") # light purple # both soy
#
# col_overlap_maize <- c("grey90","#fec44f","#e31a1c","#41ab5d") # yellow (E), red (L), green (overlap)
# #col_overlap_maize <- c("grey90",brewer.pal(n=3,name='Paired'))
