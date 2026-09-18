pacman::p_load(dplyr,ggplot2,tidyr,stringr,lubridate)

args <- commandArgs(trailingOnly = TRUE)

taxonomy_abundances=args[1]
output_file=args[2]
temperature_data=args[3]

## read in taxonomy linked to abundances
merged_w_tax = read.csv(taxonomy_abundances)

temps_ptown = read.csv(temperature_data) %>%
    dplyr::mutate(Date=as.POSIXct(tolower(STAT_ARRIV),format="%d-%b-%Y %H:%M:%S")) %>%
    dplyr::group_by(Date,STAT_ID) %>%
    dplyr::summarize(temperature=mean(FFU_RESULT)) %>%
    dplyr::mutate(DateMonthYr = as.POSIXct(paste0(format(Date,"%Y-%m"),"-01"),
                                          format="%Y-%m-%d")) %>%
    dplyr::arrange(DateMonthYr) %>%
    dplyr::mutate(Year = year(DateMonthYr), Month=month(DateMonthYr))

## summarize by transporter and Order
summarized_transporter = merged_w_tax %>%
    dplyr::mutate(Order=trimws(Order)) %>% 
    dplyr::group_by(Order,Sample)%>%
    dplyr::summarize(TPM = sum(TPM)) %>%
    tidyr::separate(Sample,sep="_",into=c("SampleNum","Rep",
                                          "Station","Year","Month"),
                    remove=FALSE) %>%
    dplyr::mutate(Month = str_replace_all(Month, "SEPT", "SEP")) %>%
    dplyr::mutate(SampleDate=as.POSIXct(strptime(tolower(paste0(Year,
                                                                "-",Month,
                                                                "-01")),
                                                 format="%Y-%b-%d")))

for_calculation = summarized_transporter %>%
    dplyr::left_join(temps_ptown,by=c("Station"="STAT_ID",
                                      "SampleDate"="DateMonthYr")) %>%
    dplyr::filter(!is.na(temperature)) %>%
    dplyr::group_by(Order) %>%
    dplyr::slice_max(order_by = TPM, n = 5) %>%
    dplyr::group_by(Order) %>%
    dplyr::summarize(mean_temp = mean(temperature,na.rm=TRUE),
                     std_temp = sd(temperature,na.rm=TRUE))

write.csv(for_calculation,output_file)