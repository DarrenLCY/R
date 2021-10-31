library(geosphere) # for distm()
library(gtools) # for combinations()
library(ggmap) # for register_google()

# Please type your API key here
# You can register one through Google. A credit card is needed, but
# you will not be charged unless you're using your API key excessively. 
register_google(key = "Please type your API key here")

setwd("C:/Users/Owner/Documents/HOTS Database/Final Paper")

br_diff <- read.csv("metadata_with_variables - br diff - SAMPLE.csv")

####################################################
## group plot IDs that are close together         ##
####################################################

## choosing Tasmania as a reference ##

# find furthest distance between loggers deployed in Tasmania

tas_df <-
  br_diff[which(substr(br_diff$loc_name, 1, 8) == "Tasmania"), ]
list_of_combinations <-
  as.data.frame(combinations(length(tas_df$plot_id), 2,
                             tas_df$plot_id))
list_of_combinations$distance <- NA
for (i in 1:length(list_of_combinations$distance)) {
  list_of_combinations$distance[i] <-
    distm(
      c(br_diff$long_in_degrees[which(br_diff$plot_id == list_of_combinations$V1[i])], br_diff$lat_in_degrees[which(br_diff$plot_id == list_of_combinations$V1[i])]),
      c(br_diff$long_in_degrees[which(br_diff$plot_id == list_of_combinations$V2[i])], br_diff$lat_in_degrees[which(br_diff$plot_id == list_of_combinations$V2[i])]),
      fun = distHaversine
    )
}
max(list_of_combinations$distance) 

# choose 500 km as furthest distance for points within same spatial region
diameter <- 500000

list_of_sites <-
  unique(ifelse(
    grepl(", ", br_diff$loc_name, fixed = TRUE),
    gsub(",.*$", "", br_diff$loc_name),
    br_diff$loc_name
  ))

# in this df, each row represents the reference coordinates, around which all
# plot ids within a 250 km radius will fall in the same spatial region
ref_coordinates_df <- data.frame(site = vector(), lat = vector(), lon = vector())

br_diff$spatial_blocks <- NA
# first, set up spatial regions with reference coordinates
for (i in list_of_sites) {
  assign(gsub(" ", "_", paste0(i, "_df")),
         br_diff[grep(i, br_diff$loc_name), ])
  
  while (length(get(gsub(" ", "_", paste0(i, "_df")))$plot_id) != 0) {
    if (length(get(gsub(" ", "_", paste0(i, "_df")))$plot_id) == 1) {
      ref_coordinates_df <- rbind(
        ref_coordinates_df,
        data.frame(
          site = i,
          lat = get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[1],
          lon = get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[1]
        )
      )
      assign(gsub(" ", "_", paste0(i, "_df")),
             get(gsub(" ", "_", paste0(i, "_df")))[-1,])
    } else {
      list_of_combinations <- as.data.frame(combinations(length(get(
        gsub(" ", "_", paste0(i, "_df"))
      )$plot_id), 2,
      get(gsub(
        " ", "_", paste0(i, "_df")
      ))$plot_id))
      list_of_combinations$distance <- NA
      for (j in 1:length(list_of_combinations$distance)) {
        list_of_combinations$distance[j] <-
          distm(
            c(br_diff$long_in_degrees[which(br_diff$plot_id == list_of_combinations$V1[j])], br_diff$lat_in_degrees[which(br_diff$plot_id == list_of_combinations$V1[j])]),
            c(br_diff$long_in_degrees[which(br_diff$plot_id == list_of_combinations$V2[j])], br_diff$lat_in_degrees[which(br_diff$plot_id == list_of_combinations$V2[j])]),
            fun = distHaversine
          )
      }
      
      # format of coordinate is lat then lon
      xtreme_coordinate_1 <- c(get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[which(get(gsub(" ", "_", paste0(i, "_df")))$plot_id
                                                                                          == list_of_combinations$V1[which(list_of_combinations$distance
                                                                                                                           == max(list_of_combinations$distance))[1]])],
                               get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[which(get(gsub(" ", "_", paste0(i, "_df")))$plot_id
                                                                                           == list_of_combinations$V1[which(list_of_combinations$distance
                                                                                                                            == max(list_of_combinations$distance))[1]])])
      
      xtreme_coordinate_2 <- c(get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[which(get(gsub(" ", "_", paste0(i, "_df")))$plot_id
                                                                                          == list_of_combinations$V2[which(list_of_combinations$distance
                                                                                                                           == max(list_of_combinations$distance))[1]])],
                               get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[which(get(gsub(" ", "_", paste0(i, "_df")))$plot_id
                                                                                           == list_of_combinations$V2[which(list_of_combinations$distance
                                                                                                                            == max(list_of_combinations$distance))[1]])])
      
      if (max(list_of_combinations$distance) > diameter) {
        remove_indices <- vector()
        for (k in 1:length(get(gsub(" ", "_", paste0(i, "_df")))$plot_id)) {
          if (distm(
            c(
              get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[k],
              get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[k]
            ),
            c(xtreme_coordinate_1[2], xtreme_coordinate_1[1]),
            fun = distHaversine
          ) <
          diameter / 2) {
            remove_indices <- append(remove_indices, k)
          }
        }
        assign(gsub(" ", "_", paste0(i, "_df")),
               get(gsub(" ", "_", paste0(i, "_df")))[-remove_indices,])
        
        #can use k again because same level
        remove_indices <- vector()
        for (k in 1:length(get(gsub(" ", "_", paste0(i, "_df")))$plot_id)) {
          if (distm(
            c(
              get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[k],
              get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[k]
            ),
            c(xtreme_coordinate_2[2], xtreme_coordinate_2[1]),
            fun = distHaversine
          ) <
          diameter / 2) {
            remove_indices <- append(remove_indices, k)
          }
        }
        assign(gsub(" ", "_", paste0(i, "_df")),
               get(gsub(" ", "_", paste0(i, "_df")))[-remove_indices,])
        
        ref_coordinates_df <- rbind(
          ref_coordinates_df,
          data.frame(
            site = i,
            lat = xtreme_coordinate_1[1],
            lon = xtreme_coordinate_1[2]
          )
        )
        
        ref_coordinates_df <- rbind(
          ref_coordinates_df,
          data.frame(
            site = i,
            lat = xtreme_coordinate_2[1],
            lon = xtreme_coordinate_2[2]
          )
        )
      } else {
        mid_coordinate <-
          c((xtreme_coordinate_1[1] + xtreme_coordinate_2[1]) / 2,
            (xtreme_coordinate_1[2] + xtreme_coordinate_2[2]) /
              2
          )
        
        remove_indices <- vector()
        for (k in 1:length(get(gsub(" ", "_", paste0(i, "_df")))$plot_id)) {
          if (distm(
            c(
              get(gsub(" ", "_", paste0(i, "_df")))$long_in_degrees[k],
              get(gsub(" ", "_", paste0(i, "_df")))$lat_in_degrees[k]
            ),
            c(mid_coordinate[2], mid_coordinate[1]),
            fun = distHaversine
          ) <
          diameter / 2) {
            remove_indices <- append(remove_indices, k)
          }
        }
        assign(gsub(" ", "_", paste0(i, "_df")),
               get(gsub(" ", "_", paste0(i, "_df")))[-remove_indices,])
        
        ref_coordinates_df <- rbind(
          ref_coordinates_df,
          data.frame(
            site = i,
            lat = mid_coordinate[1],
            lon = mid_coordinate[2]
          )
        )
      }
    }
    
  }
}

# then, classify each plot id in each spatial region
# if a plot id falls within 250 km of 2 or more spatial regions, choose closer(est) one

for (i in 1:length(br_diff$plot_id)) {
  region <-
    ifelse(
      grepl(", ", br_diff$loc_name[i], fixed = TRUE),
      gsub(",.*$", "", br_diff$loc_name[i]),
      br_diff$loc_name[i]
    )
  temp <-
    ref_coordinates_df[which(ref_coordinates_df$site == region), ]
  temp$distance <- NA
  for (j in 1:length(temp$site)) {
    temp$distance[j] <-
      distm(
        c(temp$lon[j], temp$lat[j]),
        c(br_diff$long_in_degrees[i], br_diff$lat_in_degrees[i]),
        fun = distHaversine
      )
  }
  br_diff$spatial_blocks[i] <-
    paste0(region, "_", which(temp$distance == min(temp$distance))[1])
}

## draw some maps ## 

# create circles data frame from the centers data frame
# adapted from:
# https://stackoverflow.com/questions/34183049/plot-circle-with-a-certain-radius-around-point-on-a-map-in-ggplot2
make_circles <- function(centers, radius, nPoints = 100) {
  # centers: the data frame of centers with site
  # radius: radius measured in kilometer
  # length per longitude changes with lattitude, so need correction
  
  angle <- seq(0, 2 * pi, length.out = nPoints)
  radiusLat <- radius / 111
  circleDF <-
    data.frame(site = rep(centers$site, each = nPoints))
  circleDF$lat <-
    unlist(lapply(centers$lat, function(x)
      x + radiusLat * sin(angle)))
  circleDF$radiusLon <-
    unlist(lapply(circleDF$lat, function(x)
      radius / 111 / cos(x / 57.3)))
  circleDF$refLon <-
    unlist(lapply(centers$lon, function(x)
      rep(x, nPoints)))
  circleDF$angles <- rep(angle, length(centers$site))
  circleDF$lon <-
    circleDF$refLon + circleDF$radiusLon * cos(circleDF$angles)
  return(circleDF)
}

#calculate centre
calculate_centre <- function(region_name) {
  temp <-
    ref_coordinates_df[which(ref_coordinates_df$site == region_name), ]
  x1 <- min(temp$lon)
  y1 <- min(temp$lat)
  x2 <- max(temp$lon)
  y2 <- max(temp$lat)
  centre <- c((y1 + y2) / 2, (x1 + x2) / 2)
  return(centre)
}

produce_map <- function(region_name) {
  myCoordinates <-
    ref_coordinates_df[which(ref_coordinates_df$site == region_name),]
  myCoordinates$site <- seq(1, length(myCoordinates$site))
  myCircles <- make_circles(myCoordinates, diameter / 1000 / 2)
  if (region_name == "Australia") {
    zoom_var <- 4
  } else if (region_name == "United States") {
    zoom_var <- 3
  } else if (region_name == "Chile") {
    zoom_var <- 3
  } else if (region_name == "Belize") {
    zoom_var <- 7
  }
  ggmap_object <-
    get_map(
      location = c(
        lon = calculate_centre(region_name)[2],
        lat = calculate_centre(region_name)[1]
      ),
      zoom = zoom_var,
      maptype = "satellite"
    )
  ggplot_object <-
    ggmap(ggmap_object, extent = "panel", legend = "bottomright")
  RL <-
    geom_point(aes(x = long_in_degrees, y = lat_in_degrees),
               data = br_diff[which(ifelse(
                 grepl(", ", br_diff$loc_name, fixed = TRUE),
                 gsub(",.*$", "", br_diff$loc_name),
                 br_diff$loc_name
               ) == region_name),],
               color = "#ff0000")
  ggplot_object + RL + geom_polygon(
    data = myCircles,
    aes(lon, lat, group = site),
    color = "red",
    alpha = 0
  )
}

produce_map("Australia")
produce_map("United States")
produce_map("Chile")
produce_map("Belize")
