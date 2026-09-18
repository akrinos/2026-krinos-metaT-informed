pacman::p_load(dplyr,ggplot2,tidyr,stringr)

args <- commandArgs(trailingOnly = TRUE)

taxonomy_abundances=args[1]
taxonomy_name_links=args[2]
output_file=args[3]
nutrient_data_1=args[4]
nutrient_data_2=args[5]
output_orders=args[6]

## read in taxonomy linked to abundances
merged_w_tax = read.csv(taxonomy_abundances) %>%
    dplyr::filter(E.value<0.00001)

## separate into taxonomic labels
merged_w_tax = merged_w_tax %>%
    tidyr::separate(full_classification,sep=";",into=c("Domain","Supergroup","Phylum",
                                                       "Class","Order","Family",
                                                       "Genus","Species"))

tax_abundances = read.csv(taxonomy_name_links)

## read in environmental phosphate abundances
phos_ptown = read.csv(nutrient_data_1) %>%
    dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%d-%b-%Y %H:%M:%S")) %>%
    dplyr::bind_rows(read.csv(nutrient_data_2) %>%
                     dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%Y-%m-%d")))

phos_ptown_dt = phos_ptown%>% dplyr::filter(PARAM_CODE=="PO4") %>% 
    dplyr::filter(DEPTH<5)
phos_ptown_dt = phos_ptown_dt %>% dplyr::group_by(Date,STAT_ID) %>%
    dplyr::summarize(phosphate = mean(FFU_RESULT))%>%
    dplyr::mutate(SampleDate = as.POSIXct(paste0(format(as.Date(Date, 
                                                                 format = "%Y-%m-%d"),"%Y-%m"),"-01"),
                                          format="%Y-%m-%d"))


## read in environmental nitrate/nitrite/ammonium abundances
nitr_ptown = read.csv(nutrient_data_1) %>%
    dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%d-%b-%Y %H:%M:%S")) %>%
    dplyr::bind_rows(read.csv(nutrient_data_2) %>%
                     dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%Y-%m-%d")))
nitr_ptown_dt = nitr_ptown%>% dplyr::filter(DEPTH<5)

nitr_ptown_dt = nitr_ptown_dt %>% dplyr::filter(PARAM_CODE %in% 
                                             c("NO2","NO3","NO3+NO2")) %>%
    dplyr::group_by(PARAM_CODE,STAT_ARRIV,STAT_ID,Date) %>%
    dplyr::summarize(FFU_RESULT = mean(FFU_RESULT,na.rm=TRUE)) %>%
    pivot_wider(id_cols = c("STAT_ARRIV","STAT_ID","Date"),names_from="PARAM_CODE",
                values_from="FFU_RESULT")

nitr_ptown_dt = nitr_ptown_dt %>% 
    dplyr::rename(c("NO3_NO2"="NO3+NO2")) %>%
    dplyr::mutate(nitrate = case_when(is.na(NO3_NO2) ~ NO3+NO2,
                                      TRUE~NO3_NO2))%>%
    dplyr::mutate(SampleDate = as.POSIXct(paste0(format(as.Date(Date, 
                                                                format = "%Y-%m-%d"),"%Y-%m"),"-01"),
                                          format="%Y-%m-%d"))

ammon_ptown = read.csv(nutrient_data_1) %>%
    dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%d-%b-%Y %H:%M:%S")) %>%
    dplyr::bind_rows(read.csv(nutrient_data_2) %>%
                     dplyr::mutate(Date=as.Date(STAT_ARRIV, format="%Y-%m-%d")))
ammon_ptown_dt = ammon_ptown %>% dplyr::filter(DEPTH<5)

ammon_ptown_dt = ammon_ptown_dt %>% dplyr::filter(PARAM_CODE %in% 
                                             c("NH4")) %>%
    dplyr::group_by(PARAM_CODE,STAT_ARRIV,STAT_ID,Date) %>%
    dplyr::summarize(FFU_RESULT = mean(FFU_RESULT,na.rm=TRUE)) %>%
    pivot_wider(id_cols = c("STAT_ARRIV","STAT_ID","Date"),names_from="PARAM_CODE",
                values_from="FFU_RESULT")

ammon_ptown_dt = ammon_ptown_dt %>% 
    dplyr::rename(c("ammonium"="NH4")) %>%
    dplyr::mutate(SampleDate = as.POSIXct(paste0(format(as.Date(Date, 
                                                                 format = "%Y-%m-%d"),"%Y-%m"),"-01"),
                                          format="%Y-%m-%d"))


## filter taxa that pass threshold of presence of all 4 KOs
pass_threshold = merged_w_tax %>% dplyr::mutate(Order=trimws(Order)) %>%
    dplyr::filter(trimws(Phylum)%in%c("Ochrophyta","Dinoflagellata",
                                      "Haptophyta","Chlorophyta",
                                      "Rhodophyta","Cryptophyta"))%>%
    dplyr::group_by(Order) %>%
    dplyr::summarize(KO_merged=paste(sort(unique(KO)),collapse="|"),
                     num_KOs=length(unique(KO))) %>%
    dplyr::arrange(desc(num_KOs)) %>%
    dplyr::filter(num_KOs!=1) %>%
    dplyr::filter((grepl("K14430",KO_merged)&grepl("K14683",KO_merged))&
                  (grepl("K10534",KO_merged)&grepl("K02575",KO_merged))) %>%
    dplyr::filter(!is.na(Order)&!
                  (Order%in%c("Cyclotrichiida","Acanthoecida","Pedinellales",
                              "Flavobacteriales","Stramenopiles_XXX",
                              "Pinguiochrysidaceae")))

## filter based on whether at least 30 samples have reported abundances

passing_orders = merged_w_tax %>%
    dplyr::mutate(Order=trimws(Order)) %>% dplyr::group_by(Order,KO,Sample)%>%
    dplyr::summarize(TPM = sum(V4)) %>%
    tidyr::separate(Sample,sep="_",into=c("SampleNum","Rep","Station","Year","Month"),remove=FALSE) %>%
    dplyr::mutate(transporter_annot=KO)%>%
    dplyr::mutate(Month = str_replace_all(Month, "SEPT", "SEP")) %>%
    dplyr::mutate(SampleDate=as.POSIXct(strptime(tolower(paste0(Year,"-",Month,"-01")),
                                                 format="%Y-%b-%d"))) %>%
    dplyr::left_join(phos_ptown_dt,by=c("Station"="STAT_ID","SampleDate")) %>%
    dplyr::filter((transporter_annot %in% c("K14430","K14683")))%>%
    dplyr::mutate(transporter_interp = case_when(transporter_annot=="K14430"~"low_affinity",
                                                 transporter_annot=="K14683"~"high_affinity")) %>%
    dplyr::mutate(propTPM=TPM) %>%
    tidyr::pivot_wider(id_cols=c("Station","Order","SampleDate","Sample","phosphate"),
                       names_from="transporter_interp",
                       values_from="propTPM") %>%
    dplyr::filter(!is.na(low_affinity)&!is.na(high_affinity)) %>%
    dplyr::group_by(Order) %>% dplyr::tally() %>% dplyr::filter(n>30)

## summarize by transporter and Order
summarized_transporter = merged_w_tax %>%
    dplyr::mutate(Order=trimws(Order),Phylum=trimws(Phylum)) %>% dplyr::group_by(Order,Phylum,KO,Sample)%>%
    dplyr::summarize(TPM = sum(V4)) %>%
    tidyr::separate(Sample,sep="_",into=c("SampleNum","Rep","Station","Year","Month"),remove=FALSE) %>%
    dplyr::mutate(transporter_annot=KO)%>%
    dplyr::mutate(Month = str_replace_all(Month, "SEPT", "SEP")) %>%
    dplyr::mutate(SampleDate=as.POSIXct(strptime(tolower(paste0(Year,"-",Month,"-01")),
                                                 format="%Y-%b-%d"))) %>%
    dplyr::filter((Order%in%pass_threshold$Order)&(Order%in%passing_orders$Order))

for_calculation = summarized_transporter %>%
    dplyr::left_join(phos_ptown_dt,by=c("Station"="STAT_ID","SampleDate")) %>%
    dplyr::filter((transporter_annot %in% c("K14430","K14683")))%>%
    dplyr::mutate(transporter_interp = case_when(transporter_annot=="K14430"~"low_affinity",
                                                 transporter_annot=="K14683"~"high_affinity")) %>%
    dplyr::mutate(propTPM=TPM) %>%
    tidyr::pivot_wider(id_cols=c("Station","Order","SampleDate","Sample","phosphate"),
                       names_from="transporter_interp",
                       values_from="propTPM") %>%
    dplyr::mutate(low_affinity=case_when(is.na(low_affinity)~0,
                                         TRUE~low_affinity),
                  high_affinity=case_when(is.na(high_affinity)~0,
                                         TRUE~high_affinity)) %>%
    dplyr::group_by(Order) %>%
    dplyr::mutate(low_to_high = (high_affinity-low_affinity)/(low_affinity+high_affinity)) %>%
    dplyr::group_by(Order,phosphate,Station,SampleDate,Sample) %>%
    dplyr::summarize(low_to_high=mean(low_to_high))

stats_frame = data.frame()
for_paper_stats_frame = data.frame()
for (order_curr in unique(for_calculation$Order)) {
    subset_fr = for_calculation %>% dplyr::filter(Order==order_curr) %>%
        dplyr::filter(!is.na(phosphate)&!is.na(low_to_high)) %>%
        dplyr::mutate(phosphate = phosphate + 0.001)
    print(order_curr)
    spearman_corr = cor.test(subset_fr$phosphate,subset_fr$low_to_high,method="spearman")
    print(spearman_corr)
    linear_res = lm(subset_fr$low_to_high~log10(subset_fr$phosphate))
    linear_res = lm(subset_fr$low_to_high~subset_fr$phosphate)
    print(summary(linear_res))
    
    print(as.numeric(predict(linear_res)))
    min_y_pred <- min(predict(linear_res))
    max_y_pred <- max(predict(linear_res))
    max_y = max(subset_fr$low_to_high)
    min_y = min(subset_fr$low_to_high)
    half_decay_y <- 0#0.5#0.5#(min_y+max_y)/2 #0.5#(min_y+max_y)/2 #(min_y_pred) / 2 # y value at 1/2 of the max.
    intercept <- as.numeric(coef(linear_res)[1])
    slope <- as.numeric(coef(linear_res)[2])
    print(half_decay_y)
    print(intercept)
    print(slope)
    log10_x_at_half_decay <- (half_decay_y - intercept) / slope
    print(log10_x_at_half_decay)
    
    x_at_half_decay <- log10_x_at_half_decay#10^(log10_x_at_half_decay)
    
    print(paste0("x intercept is ",-linear_res$coefficients[1]/linear_res$coefficients[2]))
    stats_frame = stats_frame %>%
        dplyr::bind_rows(data.frame("Taxon"=order_curr,
                                    "X_intercept"=10^(-linear_res$coefficients[1]/linear_res$coefficients[2]),
                                    "half_decay_y"=half_decay_y,
                                    "x_at_half_decay"=x_at_half_decay,
                                    "linear_R_square"=summary(linear_res)$r.squared,
                                    "linear_p_val"=summary(linear_res)$coefficients[,4][2],
                                    "linear_slope"=linear_res$coefficients[2],
                                    "linear_intercept"=linear_res$coefficients[1],
                                    "spearman_corr" = spearman_corr$estimate,
                                    "spearman_p_val" = spearman_corr$p.value))
    
    
    for_paper_stats_frame = for_paper_stats_frame %>%
        dplyr::bind_rows(data.frame("Taxon"=order_curr,
                                    "estimated_parameter"=x_at_half_decay,
                                    "exponential_slope"=linear_res$coefficients[2],
                                    "exponential_intercept"=linear_res$coefficients[1],
                                    "exponential_p_val"=summary(linear_res)$coefficients[,4][2],
                                    "exponential_R_square"=summary(linear_res)$r.squared,
                                    "spearman_corr" = spearman_corr$estimate,
                                    "spearman_p_val" = spearman_corr$p.value))
                                    
}

write.csv(for_paper_stats_frame  %>% dplyr::filter(exponential_p_val < 0.15) ,output_file)
order_list=unique(c(pass_threshold$Order,passing_orders$Order))
write.csv(data.frame("Order"=order_list),output_orders,row.names=FALSE,quote=FALSE)


