library(devtools)
library(PheWAS)
library(dplyr)
library(ggrepel)


group.colors <- c(
  "Infectious Diseases" = "#9FC2DA",
  "Hematopoietic" = "#C22322",
  "Endocrine/metabolic" ="#489030",
  "Dermatologic" = "#3062A0",
  "Digestive" = "#AFD881",
  "Sense organs" = "#9b870c",
  "Symptoms" = "#934A28",
  "Musculoskeletal" = "#E77027",
  "Genitourinary" = "#E78B89",
  "Neoplasms" = "#B9A2C9",
  "Infectious diseases" = "#EFB568",
  "Neurological" = "#4F2D83",
  "Respiratory" = "#A8A2B2",
  "Mental disorders" = "#71455D",
  "Injuries & poisonings" = "#EAC9C1",
  "Congenital anomalies "  = "#C0C0C0",
  "Pregnancy complications" = "#C0C0C0"
    
)

#READING IN DATA
dat<-as.data.frame(read.table("all_results_file.txt"))
colnames(dat)<-c("Phecode","CHROM","POS",	"ID",	"REF",	"ALT",	"A1",	"FIRTH?",	"TEST",	"OBS_CT",	"OR",	"LOG(OR)_SE",	"Z_STAT",	"P",	"ERRCODE")

bonf<- -log10(0.05/length(unique(dat$Phecode)))

#Simple function that capitalizes the first letter of a string
#x-charater, a string
capitalize <-function(x){
  paste(toupper(substring(x, 1, 1)),
        tolower(substring(x, 2, nchar(x))),
        sep = "")
}

#This function takes a dataframe and produces a manhattan plot for each SNP inside it
#df-dataframe, a phewas output, must have the following columns: 
#"Phecode","CHROM","POS",	"ID",	"REF",	"ALT",	"A1",	"FIRTH?",	"TEST",	
#"OBS_CT",	"OR",	"LOG(OR)_SE",	"Z_STAT",	"P",	"ERRCODE"
#bonf-numeric, the threshold for the statistical test to determine significance 
graphsnpsindividually<-function(df, bonf){
  dat2<-addPhecodeInfo(df) #to add phecode description and categories. REQUIRE PheWAS package . Also make sure the phecode column name is "Phecode"
  dat2<-dat2[order(dat2$group, dat2$Phecode),] #ordering dataset to group all phecodes together despite gaps
  dat2$consec <- 1:nrow(test) #adding a consec # to graph on so that all the clusters are together
  snps<-unique(df$ID)
  yhat<-ifelse(max(dat2$P)>bonf,max(dat2$P), bonf+0.10)
  #generate graph for each individual unique SNP
  for(i in 1:length(snps)){
    temp<-dat2[dat2$ID == snps[i],]
    temp %>%
      # Add a column for the direction of effect (risk or protective)
      mutate(.,shape_p = ifelse(OR>1,'risk','protective')) %>% #if you don't have OR then compute it from beta OR= exp(beta)
      # Capitalize phenotype group names
      mutate(.,group = capitalize(group)) %>%
      ggplot(aes(x=as.factor(consec), y=-log10(P), color=group, shape=shape_p)) +  geom_point(size=4) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "#FAFAFA"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line = element_line(colour = "black")
      )+
      scale_color_manual(values = group.colors) +
      scale_shape_manual(values = c("\u25BC","\u25B2")) +
      #scale_x_continuous(expand = c(0, 0), limits = c(0, xhat)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, yhat)) +
      # Add labels for significant associations (above the Bonferroni threshold)
      theme(legend.position="bottom", axis.text.y = element_text(face="bold",size=15)) +
      labs(x = "Phenotype", y = "-log10(p)", title=paste(snps[i], "PheWAS", sep = " "), shape = "Direction of Effect:") +
      # Add a horizontal line at the Bonferroni threshold
      geom_hline(yintercept=bonf, color = "red", size=0.6) + 
      annotate("text", x=175, y=bonf + 0.05, label=paste("Bonferroni Threshold: ", bonf), color = "red")
      filename<-gsub(":", "_", snps[i])
    #save the output
    ggsave(paste0(filename, ".jpg"), width = 15,  height = 10, dpi = 450)
  }
}

#This function takes a dataframe and produces a manhattan plot for each SNP inside it
#df-dataframe, a phewas output, must have the following columns: 
#"Phecode","CHROM","POS",	"ID",	"REF",	"ALT",	"A1",	"FIRTH?",	"TEST",	
#"OBS_CT",	"OR",	"LOG(OR)_SE",	"Z_STAT",	"P",	"ERRCODE"
#bonf-numeric, the threshold for the statistical test to determine significance 
#filename-character, name of graph file for output
graphsnpscluster<-function(df,bonf,filename){
  df<-addPhecodeInfo(df) #to add phecode description and categories. REQUIRE PheWAS package . Also make sure the phecode column name is "Phecode"
  df<-df[order(dat2$group, df$Phecode),] #ordering dataset to group all phecodes together despite gaps
  df$consec <- 1:nrow(test) #adding a consec # to graph on so that all the clusters are together
  df %>%
    # Add a column for the direction of effect (risk or protective)
    mutate(.,shape_p = ifelse(OR>1,'risk','protective')) %>% #if you don't have OR then compute it from beta OR= exp(beta)
    # Capitalize phenotype group names
    mutate(.,group = capitalize(group)) %>%
    ggplot(aes(x=consec, y=-log10(P), color=group, shape=shape_p)) +  geom_point(size=4) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "#FAFAFA"),
      plot.title = element_text(hjust = 0.5)
    )+
    scale_color_manual(values = group.colors) +
    scale_shape_manual(values = c("\u25BC","\u25B2")) +
    # Add labels for significant associations (above the Bonferroni threshold)
    theme(legend.position="bottom", axis.text.y = element_text(face="bold",size=15)) +
    labs(x = "Phenotype", y = "-log10(p)", title=paste("test", "PheWAS", sep = " "), shape = "Direction of Effect:") +
    # Add a horizontal line at the Bonferroni threshold
    geom_hline(yintercept=bonf, color = "red", size=0.6) + 
    annotate("text", x=175, y=bonf + 0.05, label=paste("Bonferroni Threshold: ", bonf), color = "red")
  #save the output
  ggsave(paste0(filename, ".jpg"), width = 15,  height = 10, dpi = 450)
}

graphsnpsindividually(dat, bonf)
#graphsnpscluster(dat,bonf,"testsingle")
#write.csv(dat2, file="all_results.csv")















# Calculate the Bonferroni threshold for significance
# bonf_threshold<-dat %>%
#   mutate(., P_bonf = p.adjust(P, method = 'bonferroni', n = length(P) )) %>%
#   filter(P_bonf < 0.05) %>%
#   arrange(P_bonf) %>%
#   tail(.,n =1) %>%
#   select(P)


#  dat %>%
# #   # Add Phecode descriptions and categories to the data
#    addPhecodeInfo(.) %>% #to add phecode description and categories. REQUIRE PheWAS package . Also make sure the phecode column name is "Phecode"
# #   # Add a column for the direction of effect (risk or protective)
#    mutate(.,shape_p = ifelse(OR>1,'risk','protective')) %>% #if you don't have OR then compute it from beta OR= exp(beta)
# #   # Capitalize phenotype group names
#   mutate(.,group = capitalize(group)) %>%
#   ggplot(aes(x=as.factor(Phecode), y=-log10(P), color=group, shape=shape_p)) +  geom_point(size=4) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     plot.background = element_rect(fill = "#FAFAFA"))+
#   scale_color_manual(values = group.colors) +
#   scale_shape_manual(values = c("\u25BC","\u25B2")) +
#   # Add labels for significant associations (above the Bonferroni threshold)
#   theme(legend.position="bottom", axis.text.y = element_text(face="bold",
#                                                              size=11)) +
#  # geom_label_repel(aes(label=ifelse(-log10(P) > -log10(bonf_threshold$P), paste0(description),'')), size=3.2,
# #                   min.segment.length = 0, box.padding = 1, seed = 10,  force = 0.5, max.overlaps = 250, color='black'
# #  ) +
# #  labs(x = "Phenotype", y = "-log10(p)", title=paste0("r")
# #       , shape = "Direction of Effect") +
#   # Add a horizontal line at the Bonferroni threshold
#   geom_hline(yintercept=bonf, color = "red", size=0.6)
# 
# #save the output
# ggsave("test.jpg", width = 15,  height = 10, dpi = 450)
