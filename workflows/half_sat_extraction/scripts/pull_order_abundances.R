pacman::p_load(ggplot2,dplyr,ape,kmer,dendextend,tibble,data.table,dynamicTreeCut)

args <- commandArgs(trailingOnly = TRUE)
output_file = args[1] # get output directory from run name
taxonomy_file = args[2]
order_file = args[3]

## filter by order inline
filtered_taxonomy = data.frame(fread(cmd = paste0("rg -f ",
                                                  order_file," ",
                                                  taxonomy_file))) %>%
    tidyr::separate(V4,sep=";",
                    into=c("Domain","Supergroup","Phylum","Class",
                           "Order","Family","Genus","Species")) %>%
    dplyr::rename(c("transcript_name"="V2"))

order_list=read.csv(order_file)

filt_transcript_file="order_filtered_transcripts.txt"
write.csv(filtered_taxonomy%>%dplyr::select(transcript_name),
          filt_transcript_file,
          row.names=FALSE,quote=FALSE)

## pull the abundances of the transcripts of interest
sample_dir="/oscar/home/akrinos/data/arianna/2026-krinos-diatom-model/14-MAD-mapping/salmon_sample"

all_salmon=data.frame()
for (sample in list.files(sample_dir)) {
    dt = fread(cmd = paste0("rg -f ",filt_transcript_file," ",file.path(sample_dir,sample,"quant.sf")))
    print(head(dt))
    flush.console()
    dt = data.frame(dt) %>%
        dplyr::right_join(filtered_taxonomy %>% dplyr::select(transcript_name,
                                                              Order),
                          by=c("V1"="transcript_name"))
    dt = dt %>% dplyr::mutate(Sample=unlist(strsplit(sample,"_quant"))[1]) %>%
        dplyr::filter(V4>0) %>%
        dplyr::mutate(Order=trimws(Order))%>%
        dplyr::group_by(Order,Sample) %>%
        dplyr::summarize(V4=sum(V4),V5=sum(V5))%>%
        dplyr::rename(c("TPM"="V4","NumReads"="V5"))
        #dplyr::filter(V4>1)
    print(nrow(dt))
    print(head(dt))
    if (nrow(dt) > 0) {
        all_salmon=all_salmon%>%
            dplyr::bind_rows(dt)
    }
}
                          
## write to file
write.csv(all_salmon,output_file)
