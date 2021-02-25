# ECMWF data
library(tidyverse)
library(ncdf4)
library(lubridate)
library(ecmwfr)

# Set api login for CDS
# need a key
wf_set_key(user = "21600", service="cds")


dat_levels <- tibble()
for (yy in 2010:2020) {
  for (mm in 1:12) {
# Request the pressure level temperature data
request <- list("dataset_short_name" = "reanalysis-era5-pressure-levels",
                "product_type"   = "reanalysis",
                "variable"       = "temperature",
                "pressure_level" = c("1000","975","950","925","900","875","850"),
                "year"           = as.character(yy),
                "month"          = str_pad(as.character(mm), width=2, pad="0"),
                "day"            = c("01","02","03","04","05","06","07","08","09","10","11",
                                     "12","13","14","15","16","17","18","19","20","21","22",
                                     "23","24","25","26","27","28","29","30","31"),
                "time"           = c("00:00","06:00","12:00","18:00"),
                "area"           = "45.5/-75.75/45.4/-75.75",
                "format"         = "netcdf",
                "target"         = "era5_temps.nc")

# Download the file
n_retries <- 0
while(TRUE) {
  err <- try(wf_request(user = "21600",
                        request = request,   
                        transfer = TRUE,  
                        path = paste("./data"),
                        verbose = FALSE)
  )
  if (class(err) != "try-error") break
  print(paste("Download error; retrying ",n_retries," times"))
  flush.console()
  n_retries <- n_retries + 1
}

nc_file <- nc_open("./data/era5_temps.nc")
lev_ <- ncvar_get(nc_file, "level")
time_ <- ncvar_get(nc_file, "time")
temp_ <- ncvar_get(nc_file, "t")
nc_close(nc_file)

# Time is in hours since 1900-01-01 00:00
# Convert to date
time_f <- as.POSIXct("1900-01-01 00:00", tz="GMT") + as.difftime(time_, units="hours")

# pressure level data
dat_pl <- tibble(
  lev = rep(lev_, length(time_)),
  time = rep(time_f, each = length(lev_)),
  temp = as.numeric(temp_)
)

# Now get surface pressure and wind and boundary layer height
request <- list("dataset_short_name"        = "reanalysis-era5-single-levels",
                "product_type"   = "reanalysis",
                "variable"       = c("surface_pressure","10m_u_component_of_wind","10m_v_component_of_wind","2m_temperature","boundary_layer_height"),
                "year"           = as.character(yy),
                "month"          = str_pad(as.character(mm), width=2, pad="0"),
                "day"            = c("01","02","03","04","05","06","07","08","09","10","11",
                                     "12","13","14","15","16","17","18","19","20","21","22",
                                     "23","24","25","26","27","28","29","30","31"),
                "time"           = c("00:00","06:00","12:00","18:00"),
                "area"           = "45.5/-75.75/45.4/-75.75",
                "format"         = "netcdf",
                "target"         = "ecmwf_surface.nc")

# Download the file
while(TRUE) {
  err <- try(wf_request(user = "21600",
                        request = request,   
                        transfer = TRUE,  
                        path = paste("./data"),
                        verbose = FALSE)
  )
  if (class(err) != "try-error") break
  print(paste("Download error; retrying ",n_retries," times"))
  flush.console()
  n_retries <- n_retries + 1
}

nc_file <- nc_open("./data/ecmwf_surface.nc")
time_ <- ncvar_get(nc_file, "time")
sp_ <- ncvar_get(nc_file, "sp")
u10_ <- ncvar_get(nc_file, "u10")
v10_ <- ncvar_get(nc_file, "v10")
blh_ <- ncvar_get(nc_file, "blh")
t2m_ <- ncvar_get(nc_file, "t2m")
nc_close(nc_file)

# Time is in hours since 1900-01-01 00:00
# Convert to date
time_f <- as.POSIXct("1900-01-01 00:00",tz="GMT") + as.difftime(time_, units="hours")

dat_sp <- tibble(
  time = time_f,
  sp = as.numeric(sp_)/100,
  u10 = as.numeric(u10_),
  v10 = as.numeric(v10_),
  blh = as.numeric(blh_),
  t2m = as.numeric(t2m_)
)


# measure inversions 
# this first measure includes surface temperature in the definition of inversions
window = 5 # Window to drop model temperature around surface in hPa
inversions <- dat_pl %>% 
  bind_rows(dat_sp %>% select(time,lev=sp,temp=t2m)) %>%
  inner_join(dat_sp) %>%
  filter(lev > sp - 100, # only look at bottom 100 hPa
         lev <= sp, # only look at above surface layers
         # only look at model levels outside window
         (lev < sp - window | lev == sp)) %>%
  arrange(time, desc(lev)) %>%
  group_by(time) %>%
  mutate(
    lapse = lag(temp) - temp,
    pd = lag(lev) - lev,
    surf_diff_temp = ifelse(lev != first(lev), first(temp) - temp, NA),
    surf_diff_pres = ifelse(lev != first(lev), first(lev) - lev, NA)) %>%
  summarise(
    # Inversions have negative lapse rates
    lowest_inv = (first(temp) - nth(temp,2)) / (first(lev) - nth(lev,2)),
    any_inv = min(lapse / pd, na.rm=TRUE),
    surf_inv = min(surf_diff_temp / surf_diff_pres, na.rm=TRUE),
    u10 = mean(u10),
    v10 = mean(v10),
    blh = mean(blh),
    t2m = mean(t2m)) 

# this second approach does not include surface temperature in defnition of inversions
inversions_nosurf <- dat_pl %>% 
  inner_join(dat_sp) %>%
  filter(lev > sp - 100, # only look at bottom 100 hPa
         lev <= sp # only look at above surface layers
         ) %>%
  arrange(time, desc(lev)) %>%
  group_by(time) %>%
  mutate(
    lapse = lag(temp) - temp,
    pd = lag(lev) - lev,
    surf_diff_temp = ifelse(lev != first(lev), first(temp) - temp, NA),
    surf_diff_pres = ifelse(lev != first(lev), first(lev) - lev, NA)) %>%
  summarise(
    # Inversions have negative lapse rates
    lowest_inv_nosurf = (first(temp) - nth(temp,2)) / (first(lev) - nth(lev,2)),
    any_inv_nosurf = min(lapse / pd, na.rm=TRUE),
    surf_inv_nosurf = min(surf_diff_temp / surf_diff_pres, na.rm=TRUE))

inversions <- inner_join(
  inversions, inversions_nosurf
)

dat_levels <- bind_rows(dat_levels, inversions)

}
}

# Adjust time zone
dat <- dat_levels %>%
  mutate(local_time = format(time, tz="EST", usetz=TRUE),
         local_time = as.POSIXct(local_time)) %>%
  mutate(hour = hour(local_time),
         year = year(local_time),
         month = month(local_time),
         day = day(local_time)) %>%
  select(-time,
         -local_time)

write_csv(dat, "./data/inversion_data.csv")

