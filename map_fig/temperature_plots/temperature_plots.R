# building a sea and air temperature plots from NOAA weather stations

library(tidyverse)
library(lubridate)
library(tibbletime)


# load in different years with data. 

# Robbinston, closest is Station PSBM1 - 8410140 - Eastport, ME
# https://www.ndbc.noaa.gov/station_history.php?station=psbm1

RB_2020 = read.delim("../temp_data/Robbinston/psbm1h2020.txt.gz", na.strings='999', sep = "", header = 1)
RB_2021 = read.delim("../temp_data/Robbinston/psbm1h2021.txt.gz", na.strings='999', sep = "", header = 1)

RB_temp = bind_rows(
                    RB_2020) %>%
  dplyr::slice(-1) %>%
  mutate(Date = paste(paste(X.YY, MM, DD, sep = "-"), paste(hh, mm, sep = ":"), sep = "")) %>%
  mutate(Date = strptime(Date, format = "%Y-%m-%d %H:%M")) %>%
  select(Date, WTMP, ATMP) %>%
  mutate(WTMP = as.numeric(WTMP),
         ATMP = as.numeric(ATMP)) %>%
  as_tbl_time(Date) %>%
  mutate(Date = as.POSIXct(Date)) %>%
  collapse_by("hourly") %>%
  group_by(Date) %>%
  dplyr::summarise(Water_Temperature = mean(WTMP),
                   Air_Temperature = mean(ATMP)) %>%
  filter(Water_Temperature < 50, 
         Air_Temperature < 50)

# Kettle cove, closest station is CASM1 - 8418150 - Portland, ME
# https://www.ndbc.noaa.gov/station_page.php?station=casm1

KC_2020 = read.delim("../temp_data/Kettle_cove/casm1h2020.txt.gz", na.strings='999', sep = "", header = 1)
KC_2021 = read.delim("../temp_data/Kettle_cove/casm1h2021.txt.gz", na.strings='999', sep = "", header = 1)
KC_temp = bind_rows(
                    KC_2020) %>%
  dplyr::slice(-1) %>%
  mutate(Date = paste(paste(X.YY, MM, DD, sep = "-"), paste(hh, mm, sep = ":"), sep = "")) %>%
  mutate(Date = strptime(Date, format = "%Y-%m-%d %H:%M")) %>%
  select(Date, WTMP, ATMP) %>%
  mutate(WTMP = as.numeric(WTMP),
         ATMP = as.numeric(ATMP)) %>%
  as_tbl_time(Date) %>%
  mutate(Date = as.POSIXct(Date)) %>%
  collapse_by("hourly") %>%
  group_by(Date) %>%
  dplyr::summarise(Water_Temperature = mean(WTMP),
                   Air_Temperature = mean(ATMP)) %>%
  filter(Water_Temperature < 50,
         Air_Temperature < 50)

# Beverly, closest station is BHBM3 - 8443970 - Boston, MA
# https://www.ndbc.noaa.gov/station_page.php?station=bhbm3

BV_2020 = read.delim("../temp_data/Beverly/bhbm3h2020.txt.gz", na.strings='999', sep = "", header = 1)
BV_2021 = read.delim("../temp_data/Beverly/bhbm3h2021.txt.gz", na.strings='999', sep = "", header = 1)
BV_temp = bind_rows(
                    BV_2020) %>%
  dplyr::slice(-1) %>%
  mutate(Date = paste(paste(X.YY, MM, DD, sep = "-"), paste(hh, mm, sep = ":"), sep = "")) %>%
  mutate(Date = strptime(Date, format = "%Y-%m-%d %H:%M")) %>%
  select(Date, WTMP, ATMP) %>%
  mutate(WTMP = as.numeric(WTMP),
         ATMP = as.numeric(ATMP)) %>%
  as_tbl_time(Date) %>%
  mutate(Date = as.POSIXct(Date)) %>%
  collapse_by("hourly") %>%
  group_by(Date) %>%
  dplyr::summarise(Water_Temperature = mean(WTMP),
                   Air_Temperature = mean(ATMP)) %>%
  filter(Water_Temperature < 50, 
         Air_Temperature < 50)

# Newport, closest station is NWPR1 - 8452660 - Newport, RI
# https://www.ndbc.noaa.gov/station_page.php?station=nwpr1

NP_2020 = read.delim("../temp_data/Newport/nwpr1h2020.txt.gz", na.strings='999', sep = "", header = 1)
NP_2021 = read.delim("../temp_data/Newport/nwpr1h2021.txt.gz", na.strings='999', sep = "", header = 1)
NP_temp = bind_rows(
                    NP_2020) %>%
  dplyr::slice(-1) %>%
  mutate(Date = paste(paste(X.YY, MM, DD, sep = "-"), paste(hh, mm, sep = ":"), sep = "")) %>%
  mutate(Date = strptime(Date, format = "%Y-%m-%d %H:%M")) %>%
  select(Date, WTMP, ATMP) %>%
  mutate(WTMP = as.numeric(WTMP),
         ATMP = as.numeric(ATMP)) %>%
  as_tbl_time(Date) %>%
  mutate(Date = as.POSIXct(Date)) %>%
  collapse_by("hourly") %>%
  group_by(Date) %>%
  dplyr::summarise(Water_Temperature = mean(WTMP),
                   Air_Temperature = mean(ATMP)) %>%
  filter(Water_Temperature < 50, 
         Air_Temperature < 50)

# Cape May, closest station is CMAN4 - 8536110 - Cape May, NJ
# https://www.ndbc.noaa.gov/station_page.php?station=cman4

CM_2020 = read.delim("../temp_data/Cape_may/cman4h2020.txt.gz", na.strings='999', sep = "", header = 1)
CM_2021 = read.delim("../temp_data/Cape_may/cman4h2021.txt.gz", na.strings='999', sep = "", header = 1)
CM_temp = bind_rows(
                    CM_2020) %>%
  dplyr::slice(-1) %>%
  mutate(Date = paste(paste(X.YY, MM, DD, sep = "-"), paste(hh, mm, sep = ":"), sep = "")) %>%
  mutate(Date = strptime(Date, format = "%Y-%m-%d %H:%M")) %>%
  select(Date, WTMP, ATMP) %>%
  mutate(WTMP = as.numeric(WTMP),
         ATMP = as.numeric(ATMP)) %>%
  as_tbl_time(Date) %>%
  mutate(Date = as.POSIXct(Date)) %>%
  collapse_by("hourly") %>%
  group_by(Date) %>%
  dplyr::summarise(Water_Temperature = mean(WTMP),
                   Air_Temperature = mean(ATMP)) %>%
  filter(Water_Temperature < 50, 
         Air_Temperature < 50)

### Plots

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")


RB_plot = 
  ggplot(RB_temp, aes(x = Date, y = Air_Temperature, color = Air_Temperature)) +
  geom_line(size = 1.1) +
  scale_color_gradient2(low = "dodgerblue1", mid = "lightyellow2", high = "orangered1", midpoint = 15) +
  labs(x = "Year", y = "Temperature (°C)", color = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_line(data = RB_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature), color = "#264D59") 

ggsave("Robbinston_profile.pdf", plot = RB_plot, width = 4, height = 4, units = "in")

KC_plot = 
  ggplot(KC_temp, aes(x = Date, y = Air_Temperature, color = Air_Temperature)) +
  geom_line(size = 1.1) +
  scale_color_gradient2(low = "dodgerblue1", mid = "lightyellow2", high = "orangered1", midpoint = 15) +
  labs(x = "Year", y = "Temperature (°C)", color = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_line(data = KC_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature), color = "#43978D") 

ggsave("Kettle_Cove_profile.pdf", plot = KC_plot, width = 4, height = 4, units = "in")

BV_plot = 
  ggplot(BV_temp, aes(x = Date, y = Air_Temperature, color = Air_Temperature)) +
  geom_line(size = 1.1) +
  scale_color_gradient2(low = "dodgerblue1", mid = "lightyellow2", high = "orangered1", midpoint = 15) +
  labs(x = "Year", y = "Temperature (°C)", color = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_line(data = BV_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature), color = "#F9E07F") 

ggsave("Beverly_profile.pdf", plot = BV_plot, width = 4, height = 4, units = "in")

NP_plot = 
  ggplot(NP_temp, aes(x = Date, y = Air_Temperature, color = Air_Temperature)) +
  geom_line(size = 1.1) +
  scale_color_gradient2(low = "dodgerblue1", mid = "lightyellow2", high = "orangered1", midpoint = 15) +
  labs(x = "Year", y = "Temperature (°C)", color = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_line(data = NP_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature), color = "#F9AD6A") 

ggsave("Newport_profile.pdf", plot = NP_plot, width = 4, height = 4, units = "in")

CM_plot = 
  ggplot(CM_temp, aes(x = Date, y = Air_Temperature, color = Air_Temperature)) +
  geom_line(size = 1.1) +
  scale_color_gradient2(low = "dodgerblue1", mid = "lightyellow2", high = "orangered1", midpoint = 15) +
  labs(x = "Year", y = "Temperature (°C)", color = "Temperature (°C)") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_line(data =CM_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature), color = "#D46C4E") +
  geom_smooth(data = CM_temp, aes(x = Date, y = Water_Temperature, color = Water_Temperature))

ggsave("Cape_May_profile.pdf", plot = CM_plot, width = 4, height = 4, units = "in")

Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")


ggplot(RB_temp, aes(x = Date, y = Water_Temperature)) +
  geom_line(data = RB_temp, aes(x = Date, y = Water_Temperature), color = "#264D59", alpha = 0.3) +
  geom_line(data = KC_temp, aes(x = Date, y = Water_Temperature), color = "#43978D", alpha = 0.3) +
  geom_line(data = BV_temp, aes(x = Date, y = Water_Temperature), color = "#F9E07F", alpha = 0.3) +
  geom_line(data = NP_temp, aes(x = Date, y = Water_Temperature), color = "#F9AD6A", alpha = 0.3) +
  geom_line(data = CM_temp, aes(x = Date, y = Water_Temperature), color = "#D46C4E", alpha = 0.3) + 
  geom_smooth(data = RB_temp, aes(x = Date, y = Water_Temperature, color = "Robbinston")) +
  geom_smooth(data = KC_temp, aes(x = Date, y = Water_Temperature, color = "Kettle Cove")) +
  geom_smooth(data = BV_temp, aes(x = Date, y = Water_Temperature, color = "Beverly")) +
  geom_smooth(data = NP_temp, aes(x = Date, y = Water_Temperature, color = "Newport")) +
  geom_smooth(data = CM_temp, aes(x = Date, y = Water_Temperature, color = "Cape May")) + 
  scale_colour_manual(name="Site", values=c("Robbinston" =  "#264D59",
                                            "Kettle Cove" = "#43978D",
                                            "Beverly"=      "#F9E07F",
                                            "Newport" =     "#F9AD6A",
                                            "Cape May" =    "#D46C4E")) +
  theme_classic()

ggsave("temp_plot5.png", last_plot())

Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")




