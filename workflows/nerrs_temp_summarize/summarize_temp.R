pacman::p_load(dplyr,ggplot2)#,rnaturalearth)

data_nerrs_download_1 = read.csv("../../../data_nerrs/635425.csv")
data_nerrs_download_2 = read.csv("../../../data_nerrs/692129.csv")
data_nerrs_download_3 = read.csv("../../../data_nerrs/234241.csv")
data_nerrs_download_4 = read.csv("../../../data_nerrs/225083.csv")
data_nerrs_download_5 = read.csv("../../../data_nerrs/100351.csv")

data_nerrs_download_1_temp = read.csv("../../../data_nerrs/170944.csv")
data_nerrs_download_2_temp = read.csv("../../../data_nerrs/244450.csv")
data_nerrs_download_3_temp = read.csv("../../../data_nerrs/343369.csv")
data_nerrs_download_4_temp = read.csv("../../../data_nerrs/467564.csv")
data_nerrs_download_5_temp = read.csv("../../../data_nerrs/511994.csv")
data_nerrs_download_6_temp = read.csv("../../../data_nerrs/404467.csv")

metadata_nerrs_download_1 = read.csv("../../../data_nerrs/635425_sampling_stations_utf_8.csv")
metadata_nerrs_download_2 = read.csv("../../../data_nerrs/692129_sampling_stations_utf_8.csv")
metadata_nerrs_download_3 = read.csv("../../../data_nerrs/234241_sampling_stations_utf_8.csv")
metadata_nerrs_download_4 = read.csv("../../../data_nerrs/225083_sampling_stations_utf_8.csv")
metadata_nerrs_download_5 = read.csv("../../../data_nerrs/100351_sampling_stations_utf_8.csv")

metadata_nerrs_download_1_temp = read.csv("../../../data_nerrs/170944_sampling_stations_utf_8.csv")
metadata_nerrs_download_2_temp = read.csv("../../../data_nerrs/244450_sampling_stations_utf_8.csv")
metadata_nerrs_download_3_temp = read.csv("../../../data_nerrs/343369_sampling_stations_utf_8.csv")
metadata_nerrs_download_4_temp = read.csv("../../../data_nerrs/467564_sampling_stations_utf_8.csv")
metadata_nerrs_download_5_temp = read.csv("../../../data_nerrs/511994_sampling_stations_utf_8.csv")
metadata_nerrs_download_6_temp = read.csv("../../../data_nerrs/404467_sampling_stations_utf_8.csv")

merged_nerrs_1_temp = data_nerrs_download_1_temp %>% dplyr::left_join(metadata_nerrs_download_1_temp,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_2_temp = data_nerrs_download_2_temp %>% dplyr::left_join(metadata_nerrs_download_2_temp,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_3_temp = data_nerrs_download_3_temp %>% dplyr::left_join(metadata_nerrs_download_3_temp,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_4_temp = data_nerrs_download_4_temp %>% dplyr::left_join(metadata_nerrs_download_4_temp,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_5_temp = data_nerrs_download_5_temp %>% dplyr::left_join(metadata_nerrs_download_5_temp,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_6_temp = data_nerrs_download_6_temp %>% dplyr::left_join(metadata_nerrs_download_6_temp,
                                           by=c("StationCode"="Station.Code"))

merged_nerrs_1 = data_nerrs_download_1 %>% dplyr::left_join(metadata_nerrs_download_1,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_2 = data_nerrs_download_2 %>% dplyr::left_join(metadata_nerrs_download_2,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_3 = data_nerrs_download_3 %>% dplyr::left_join(metadata_nerrs_download_3,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_4 = data_nerrs_download_4 %>% dplyr::left_join(metadata_nerrs_download_4,
                                           by=c("StationCode"="Station.Code"))
merged_nerrs_5 = data_nerrs_download_5 %>% dplyr::left_join(metadata_nerrs_download_5,
                                           by=c("StationCode"="Station.Code"))

all_metadata = metadata_nerrs_download_1 %>%
    dplyr::bind_rows(metadata_nerrs_download_2)%>%
    dplyr::bind_rows(metadata_nerrs_download_3)%>%
    dplyr::bind_rows(metadata_nerrs_download_4)%>%
    dplyr::bind_rows(metadata_nerrs_download_5)%>%
    dplyr::mutate(NERR.Site.ID=trimws(NERR.Site.ID))

all_merged = merged_nerrs_1 %>%
    dplyr::bind_rows(merged_nerrs_2)%>%
    dplyr::bind_rows(merged_nerrs_3)%>%
    dplyr::bind_rows(merged_nerrs_4)%>%
    dplyr::bind_rows(merged_nerrs_5) %>%
    dplyr::mutate(NERR.Site.ID=trimws(NERR.Site.ID))

all_metadata_temp = metadata_nerrs_download_1_temp %>%
    dplyr::bind_rows(metadata_nerrs_download_2_temp)%>%
    dplyr::bind_rows(metadata_nerrs_download_3_temp)%>%
    dplyr::bind_rows(metadata_nerrs_download_4_temp)%>%
    dplyr::bind_rows(metadata_nerrs_download_5_temp)%>%
    dplyr::bind_rows(metadata_nerrs_download_6_temp)%>%
    dplyr::mutate(NERR.Site.ID=trimws(NERR.Site.ID))

all_merged_temp = merged_nerrs_1_temp %>%
    dplyr::bind_rows(merged_nerrs_2_temp)%>%
    dplyr::bind_rows(merged_nerrs_3_temp)%>%
    dplyr::bind_rows(merged_nerrs_4_temp)%>%
    dplyr::bind_rows(merged_nerrs_5_temp) %>%
    dplyr::bind_rows(merged_nerrs_6_temp) %>%
    dplyr::mutate(NERR.Site.ID=trimws(NERR.Site.ID))

temp_values = all_merged_temp %>%
    dplyr::filter(trimws(F_Temp)%in%c("<0>","<4>","<5>"))

averaged_temp = temp_values%>%
    dplyr::mutate(DateTime=as.Date(format(as.Date(DateTimeStamp,format="%m/%d/%Y"), 
                                                     "%Y-%m-%d"))) %>%
    dplyr::group_by(DateTime,NERR.Site.ID)%>%
    dplyr::summarize(Temp=mean(Temp))

write.csv(averaged_temp,
    "../../data/nerrs_data_cleaned_temp.csv")

