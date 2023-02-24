
library(rcartocolor)
palette.colors()
library("colorspace")
install.packages("factoextra",)("rcartocolor")
#Attempt witj Batch 1#######
setwd("D:/Wolf project/Data processing")

dir <-"./RLC_rep2_byr2_Fitness/Batch_1"
Batch1_consol <-NA
Batch1_consol <- left_join( lookup,growth_batch1_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 
Batch1_consol <- left_join( Batch1_consol,ecological_sel_bottom_batch1_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,competition.x,comparison.x,type.x,G,G_log,relative_W_evo,ESB,ESB_log)%>%
  rename(s_evol_difflog_bottom = s_evol_difflog )%>%
  rename(relative_W_bottom=relative_W_evo)
Batch1_consol <- left_join( Batch1_consol,mating_unedit_batch1_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,s_evol_difflog_bottom,relative_v,competition.x,
         comparison.x,type.x,G,G_log,relative_W_bottom,ESB,ESB_log,viability,viability_log)

Batch1_consol <- left_join( Batch1_consol,ecological_sel_top_batch1_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,s_evol_difflog_bottom,relative_v,
         competition.x,comparison.x,type.x,G,G_log,relative_W_evo,relative_W_bottom,ESB,ESB_log,viability,viability_log,EST,EST_log)%>%
  rename(s_evol_difflog_top = s_evol_difflog)%>%
  rename(relative_W_top=relative_W_evo)%>%
  rename(competition = competition.x)%>%
  rename(comparison = comparison.x)%>%
  rename(type = type.x)


Batch1_final<- Batch1_consol %>%
  filter( competition == "Test" )
Batch1_Dataset <- Batch1_final %>%
  add_column(batch = "batch1")

write_csv(Batch1_Dataset, file = "./Databatch1/221202_data_batch1_final datasetconsolidated_final.csv")







#Attempt witj Batch 2#########
setwd("D:/Wolf project/Data processing")

dir <-"./RLC_rep2_byr2_Fitness/Batch_3"
batch2_consol <-NA
batch2_consol <- left_join( lookup,growth_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 
batch2_consol <- left_join( batch2_consol,ecological_sel_bottom_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,competition.x,comparison.x,type.x,G,G_log,relative_W_evo,ESB,ESB_log)%>%
  rename(s_evol_difflog_bottom = s_evol_difflog )%>%
  rename(relative_W_bottom=relative_W_evo)
batch2_consol <- left_join( batch2_consol,mating_unedit_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,s_evol_difflog_bottom,relative_v,competition.x,
         comparison.x,type.x,G,G_log,relative_W_bottom,ESB,ESB_log,viability,viability_log)

batch2_consol <- left_join( batch2_consol,ecological_sel_top_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,s_evol_difflog_bottom,relative_v,
         competition.x,comparison.x,type.x,G,G_log,relative_W_evo,relative_W_bottom,ESB,ESB_log,viability,viability_log,EST,EST_log)%>%
  rename(s_evol_difflog_top = s_evol_difflog)%>%
  rename(relative_W_top=relative_W_evo)%>%
  rename(competition = competition.x)%>%
  rename(comparison = comparison.x)%>%
  rename(type = type.x)


batch2_final<- batch2_consol %>%
  filter( competition == "Test" )
batch2_Dataset <- batch2_final %>%
  add_column(batch = "batch2")

write_csv(batch2_Dataset, file = "./Databatch2/221202_data_batch2_final datasetconsolidated_final.csv")
#Attempt witj Batch 2########
setwd("D:/Wolf project/Data processing")

dir <-"./RLC_rep2_byr2_Fitness/Batch_2"
batch2_consol <-NA
batch2_consol <- left_join( lookup,growth_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 
batch2_consol <- left_join( batch2_consol,ecological_sel_bottom_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,competition.x,comparison.x,type.x,G,G_log,relative_W_evo,ESB,ESB_log)%>%
  rename(s_evol_difflog_bottom = s_evol_difflog )%>%
  rename(relative_W_bottom=relative_W_evo)
batch2_consol <- left_join( batch2_consol,mating_unedit_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,s_evol_difflog_bottom,relative_v,competition.x,
         comparison.x,type.x,G,G_log,relative_W_bottom,ESB,ESB_log,viability,viability_log)

batch2_consol <- left_join( batch2_consol,ecological_sel_top_batch2_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,s_evol_difflog_bottom,relative_v,
         competition.x,comparison.x,type.x,G,G_log,relative_W_evo,relative_W_bottom,ESB,ESB_log,viability,viability_log,EST,EST_log)%>%
  rename(s_evol_difflog_top = s_evol_difflog)%>%
  rename(relative_W_top=relative_W_evo)%>%
  rename(competition = competition.x)%>%
  rename(comparison = comparison.x)%>%
  rename(type = type.x)


batch2_final<- batch2_consol %>%
  filter( competition == "Test" )
batch2_Dataset <- batch2_final %>%
  add_column(batch = "batch2")

write_csv(batch2_Dataset, file = "./Databatch2/221202_data_batch2_final datasetconsolidated_final.csv")





#Attempt witj Batch 3########
setwd("D:/Wolf project/Data processing")

dir <-"./RLC_rep2_byr2_Fitness/Batch_3"
batch3_consol <-NA
batch3_consol <- left_join( lookup,growth_batch3_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 



batch3_consol <- left_join( batch3_consol,  ecological_sel_bottom_batch3_normalised_std,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,competition.x,comparison.x,type.x,G,G_log,relative_W_evo,ESB,ESB_log)%>%
  rename(s_evol_difflog_bottom = s_evol_difflog )%>%
  rename(relative_W_bottom=relative_W_evo)




batch3_consol <- left_join( batch3_consol,mating_unedit_batch3_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,G,G_log,s_evol_difflog_bottom,relative_v,competition.x,
         comparison.x,type.x,G,G_log,relative_W_bottom,ESB,ESB_log,viability,viability_log)

batch3_consol <- left_join( batch3_consol,ecological_sel_top_batch3_normalised,by ="well")%>% 
  #filter( competition == "Test" )%>%
  select(genotype.x, well,m.evol,s_evol_difflog,s_evol_difflog_bottom,relative_v,
         competition.x,comparison.x,type.x,G,G_log,relative_W_evo,relative_W_bottom,ESB,ESB_log,viability,viability_log,EST,EST_log)%>%
  rename(s_evol_difflog_top = s_evol_difflog)%>%
  rename(relative_W_top=relative_W_evo)%>%
  rename(competition = competition.x)%>%
  rename(comparison = comparison.x)%>%
  rename(type = type.x)


batch3_final<- batch3_consol %>%
  filter( competition == "Test" )
batch3_Dataset <- batch3_final %>%
  add_column(batch = "batch3")

write_csv(batch3_Dataset, file = "./Databatch3/221202_data_batch3_final datasetconsolidated_final.csv")


######Creating a single Dataset######

dataset_1_2 <- rbind(Batch1_Dataset,batch2_Dataset)
Dataset_all_batches<-  rbind(dataset_1_2,batch3_Dataset)

Dataset_all_batches <- Dataset_all_batches %>%
  mutate(Gene_of_interest = case_when(
  endsWith(genotype.x, "RwBf") ~ "rep2 wild type",
  endsWith(genotype.x, "RfBw") ~ "byr2 wild type",
  endsWith(genotype.x, "RfBm") ~ "byr2 mutant",
  endsWith(genotype.x, "RmBm") ~ "Double Mutant",
  endsWith(genotype.x, "A") ~ "Wild type",
  endsWith(genotype.x, "RmBf") ~ "rep2 mutant"
))





####Growth#####
A<-Dataset_all_batches %>%
  filter( competition == 'Test') %>%
 #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = G_log,fill= Gene_of_interest), ) +
  geom_boxplot( ) +
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+
  stat_summary(fun.data = give.n, geom = "text")+
  #scale_fill_brewer()+
  #facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("A. Growth (All Batches Combined)")
ggsave("Growth_all_datasets.pdf")
Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = G_log,fill= Gene_of_interest)) +
            theme(legend.position="none") +
            #scale_color_manual(values= colorBlindGrey8 )+
            scale_fill_brewer(palette="BuPu")+
  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
  facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("B. Growth ")
ggsave("Growth_all_split_datasets.pdf")
A+B
ggsave("Growth_all_datasets_consol.pdf")
growth_total_anova<- 
  Dataset_all_batches %>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(G_log~Gene_of_interest+batch,.)

anova(growth_total_anova)
growth_total_tukey<-TukeyHSD(growth_total_anova)
growth_tukey_df<-as.data.frame(growth_total_tukey[1:1])
write_csv(growth_tukey_df, file="./RLC_rep2_byr2_Fitness/growth_total_tukey.csv")


growth_total_anova<- 
  Dataset_all_batches %>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(G_log~Gene_of_interest +batch,.)



#####Bottom Selection#######

Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = ESB_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() +   
  
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+
  
  stat_summary(fun.data = give.n, geom = "text")+
  #facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("A.Ecological Selection Bottom (All Batches combined) ")
ggsave("Ecological_Selection_Bottom_all_datasets.pdf")



Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = ESB_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() +   
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+ stat_summary(fun.data = give.n, geom = "text")+
  facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Ecological Selection Bottom ")
ggsave("Ecological_Selection_Bottom_all_split_datasets.pdf")



bottom_total_anova<- 
  Dataset_all_batches %>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(ESB_log~Gene_of_interest+batch,.)
bottom_total_tukey<-TukeyHSD(bottom_total_anova)
bottom_tukey_df<-as.data.frame(bottom_total_tukey[1:1])
write_csv(bottom_tukey_df, file="./RLC_rep2_byr2_Fitness/bottom_selection_total_tukey.csv")



#####Top Selectiom#######
Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  #na.omit( EST_log)%>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = EST_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() + 
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+   stat_summary(fun.data = give.n, geom = "text")+
  #facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("A.Ecological Selection Top(Combined All batches) ")
ggsave("Ecological_Selection_Top_all_datasets.pdf")

Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  #na.omit( EST_log)%>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = EST_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() + 
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+   stat_summary(fun.data = give.n, geom = "text")+
  facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Ecological Selection Top")
ggsave("Ecological_Selection_Top_all_split_datasets.pdf")


top_total_anova<- 
  Dataset_all_batches %>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(EST_log~Gene_of_interest+batch,.)
top_total_tukey<-TukeyHSD(top_total_anova)
top_tukey_df<-as.data.frame(top_total_tukey[1:1])
write_csv(top_tukey_df, file="./RLC_rep2_byr2_Fitness/top_total_tukey.csv")


####mating#######

Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = viability_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() +
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+   stat_summary(fun.data = give.n, geom = "text")+
  #facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("A.Mating Viability( All Batches Combined) ")
ggsave("Mating_viability_all_datasets.pdf")



Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = Gene_of_interest, y = viability_log,fill= Gene_of_interest), palette = "jco") +
  geom_boxplot() +
  theme(legend.position="none") +
  #scale_color_manual(values= colorBlindGrey8 )+
  scale_fill_brewer(palette="BuPu")+    stat_summary(fun.data = give.n, geom = "text")+
  facet_grid( .~ batch )+  
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Mating Viability ")
ggsave("Mating_viability_split_all_datasets.pdf")






mating_total_anova<- 
  Dataset_all_batches %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(viability_log~Gene_of_interest+batch,.)
mating_total_tukey<-TukeyHSD(mating_total_anova)
mating_tukey_df<-as.data.frame(mating_total_tukey[1:1])
write_csv(mating_tukey_df, file="./RLC_rep2_byr2_Fitness/mating_total_tukey.csv")




##########PCA######
Dataset_all_batches_NA<- Dataset_all_batches%>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
  select(genotype.x, well,batch,competition,type,Gene_of_interest,G,relative_W_bottom,
         ESB,viability,EST,G_log,ESB_log,EST_log,viability_log)%>%
  na.omit(.)

Dataset_all_batches_NA$G_log_z<-(Dataset_all_batches_NA$G_log - mean(Dataset_all_batches_NA$G_log))/sd(Dataset_all_batches_NA$G_log)

Dataset_all_batches_NA$viability_log_z<-(Dataset_all_batches_NA$viability_log - mean(Dataset_all_batches_NA$viability_log))/sd(Dataset_all_batches_NA$viability_log)

Dataset_all_batches_NA$ESB_log_z<-(Dataset_all_batches_NA$ESB_log - mean(Dataset_all_batches_NA$ESB_log))/sd(Dataset_all_batches_NA$ESB_log)

Dataset_all_batches_NA$EST_log_z<-(Dataset_all_batches_NA$EST_log - mean(Dataset_all_batches_NA$EST_log))/sd(Dataset_all_batches_NA$EST_log)





Allbatch.pca <- prcomp(Dataset_all_batches_NA[,c(16:19)], center = TRUE,scale. = TRUE)

allPCA<-autoplot(Allbatch.pca, data = Dataset_all_batches_NA, colour = 'Gene_of_interest',shape="batch",
                loadings=TRUE,loadings.label = TRUE,loadings.label.size = 4,x=1,y=2)+
  ggtitle("PCA ALL samples")
ggsave("Alldata_PCA_easy_byr2_rep2.pdf")
ggplotly(allPCA)

#####################


Dataset_all_batches_NA_2<- Dataset_all_batches%>%
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
  select(genotype.x, well,batch,competition,type,Gene_of_interest,G,relative_W_bottom,
         ESB,viability,EST,G_log,ESB_log,EST_log,viability_log)%>%
  na.omit(.)
Allbatch_2.pca <- prcomp(Dataset_all_batches_NA_2[,c(12:15)], center = TRUE,scale. = TRUE)

allPCA_2<-autoplot(Allbatch_2.pca, data = Dataset_all_batches_NA_2, colour = 'Gene_of_interest',shape="batch",
                 loadings=TRUE,loadings.label = TRUE,loadings.label.size = 4,x=1,y=2)+
  ggtitle("PCA ")
ggsave("Alldata_2_PCA_easy_byr2_rep2.pdf")
ggplotly(allPCA_2)
#########################
Dataset_all_batches_NA_3<- Dataset_all_batches%>%
  filter ( batch ==('batch3')) %>%
  filter( competition == 'Test') %>%
  filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
  select(genotype.x, well,batch,competition,type,Gene_of_interest,G,relative_W_bottom,
         ESB,viability,EST,G_log,ESB_log,EST_log,viability_log)%>%
  na.omit(.)
Allbatch_3.pca <- prcomp(Dataset_all_batches_NA_3[,c(12:15)], center = TRUE,scale. = TRUE)

allPCA_3<-autoplot(Allbatch_3.pca, data = Dataset_all_batches_NA_3, colour = 'Gene_of_interest',shape="batch",
                   loadings=TRUE,loadings.label = TRUE,loadings.label.size = 4,x=1,y=2)+
  ggtitle("PCA ")
ggsave("Alldata_3_PCA_easy_byr2_rep2.pdf")
ggplotly(allPCA_3)
####correlation analysis########
Dataset_all_batches_NA_2<- Dataset_all_batches%>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  filter( competition == 'Test') %>%
  #filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
  select(genotype.x, well,batch,competition,type,Gene_of_interest,G,relative_W_bottom,
         ESB,viability,EST,G_log,ESB_log,EST_log,viability_log,m.evol,s_evol_difflog_top,
         s_evol_difflog_bottom,relative_v, G,G_log,relative_W_top,relative_W_bottom)%>%
  na.omit(.)

 
 F<- Dataset_all_batches %>% 
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  select(competition,Gene_of_interest,viability_log,viability,G,G_log)%>%
  filter( competition == 'Test') %>%
  #group_by(Gene_of_interest) %>%
  na.omit( .)
 
 
 colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                 ############
 
 

   o<-Dataset_all_batches%>%
     filter ( batch %in% c('batch2' ,'batch3')) %>%
     filter( competition == 'Test') %>%
     filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
   select(competition,Gene_of_interest,viability_log,viability,G,G_log,ESB,ESB_log,EST,EST_log) 
   #%>%
   #na.omit( .)
   #(competition,Gene_of_interest,viability_log,viability,G,G_log,ESB,ESB_log,EST,EST_log)
  Correlation<- rcorr( as.matrix(o[,c(3:10)]))
n<-as.data.frame(Correlation$n)
  write_csv(n, file="./RLC_rep2_byr2_Fitness/Full_correlation_Matrix_sample_size.csv")
  
  
  
  library("Hmisc")
  ######################
  Larger <- Dataset_all_batches%>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter( competition == 'Test') %>%
    #filter(genotype.x %in% c("RfBm","RmBf","RmBm")  )%>%
    select(genotype.x, well,batch,competition,type,Gene_of_interest,G,
           ESB,viability,EST,G_log,ESB_log,EST_log,viability_log,m.evol,s_evol_difflog_top,
           s_evol_difflog_bottom,relative_v,relative_W_top,relative_W_bottom)
 
  #%>%
  #na.omit( .)
  #(competition,Gene_of_interest,viability_log,viability,G,G_log,ESB,ESB_log,EST,EST_log)
  
  Correlation_L<- rcorr( as.matrix(Larger[,c(15:20)]))
  
  N_L<-as.data.frame(Correlation_L$P)
  write_csv(N_L, file="./RLC_rep2_byr2_Fitness/Full_correlation_Matrix_significance_Larger.csv")
  
  
  
  
####consolidating medians######


growth.all.median.values<- 
  Dataset_all_batches %>% 
  #filter(type=="Wild type") %>% 
  filter( competition == 'Test') %>%
  group_by(Gene_of_interest) %>%
  dplyr::summarize( 
    G_log.median = median( G_log),
   standard.deviation.G_log = sd( G_log ),
    G.median= median(G),
   standard.deviation.G = sd( G )
  )
viability.all.median.values <- 
  Dataset_all_batches %>% 
  #filter ( batch %in% c('batch2' ,'batch3')) %>%
  select(competition,Gene_of_interest,viability_log,viability)%>%
  filter( competition == 'Test') %>%
  group_by(Gene_of_interest) %>%
  na.omit( .)%>%
  dplyr::summarize( 
    viability_log.median = median( viability_log),
   standard.deviation.viability_log = sd( viability_log ),
    viability.median = median( viability),
   standard.deviation.viability = sd( viability )
     )
ESB.all.median.values <- 
  Dataset_all_batches %>% 
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  select(competition,Gene_of_interest,ESB_log,ESB)%>%
  filter( competition == 'Test') %>%
  group_by(Gene_of_interest) %>%
  na.omit( .)%>%
  dplyr::summarize( 
    ESB_log.median = median( ESB_log),
   standard.deviation.ESB_log = sd( ESB_log ),
    ESB.median = median( ESB),
   standard.deviation.ESB = sd( ESB )
  )
EST.all.median.values <- 
  Dataset_all_batches %>% 
  filter ( batch %in% c('batch2' ,'batch3')) %>%
  select(competition,Gene_of_interest,EST_log,EST)%>%
  filter( competition == 'Test') %>%
  group_by(Gene_of_interest) %>%
  na.omit( .)%>%
  dplyr::summarize( 
    EST_log.median = median( EST_log),
   standard.deviation.EST_log = sd( EST_log ),
    EST.median = median( EST),
   standard.deviation.EST = sd( EST )
  )
  All.median.values <- growth.all.median.values%>%
  
  left_join( . , 
             EST.all.median.values, by=c( "Gene_of_interest" ) ) %>% 
    left_join( . , 
               ESB.all.median.values, by=c( "Gene_of_interest" ) ) %>%
    left_join( . , 
               viability.all.median.values , by=c( "Gene_of_interest" ) )
  
  write_csv(All.median.values, file="./RLC_rep2_byr2_Fitness/All_medians.csv")
  #EST.all.median.values <- 
  #  Dataset_all_batches %>% 
   # select(competition,Gene_of_interest,ESB_log)%>%
    #filter( competition == 'Test') %>%
    #group_by(Gene_of_interest) %>%
    #na.omit(.)%>%
    #count(.)
 #a<- Dataset_all_batches %>%
  #  filter( competition == 'Test') %>%
   # na.omit( EST_log)%>%
  ######tukey consolidated#####
  tukey_all=growth_total_tukey+bottom_total_tukey+top_total_tukey +mating_total_tukey
  summary(growth_total_tukey)
  as.data.frame(growth_total_tukey[1:1])
  ############
  top_total_anova_2<- 
    Dataset_all_batches %>%
   # filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter( competition == 'Test') %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    aov(EST_log~Gene_of_interest,.)
  TukeyHSD(top_total_anova_2)
  Dataset_all_batches %>%
    filter( competition == 'Test') %>%
    #na.omit( EST_log)%>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = Gene_of_interest, y = EST_log), palette = "jco") +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
    #facet_grid( .~ batch )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
    ggtitle("A.Ecological Selection Top(Combined All batches) ")
  
  
##### All batch figures ######
  Dataset_all_batches %>%
    filter( competition == 'Test') %>%
    #filter ( batch == 'batch1') %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    #filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = type , y = log(m.evol))) +
    #theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    #scale_fill_brewer(palette="Set1")+
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( batch~comparison)+     
    #facet_grid(  . ~ comparison )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
    ggtitle("Growth ")
  ggsave("Growth_all_batches.pdf")
  
 
  
  
  Dataset_all_batches %>%
    filter( competition == 'Test') %>%
    #filter ( batch == 'batch1') %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    #filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = type , y = log(relative_v))) +
    #theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    #scale_fill_brewer(palette="Set1")+
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( batch~comparison)+     
    #facet_grid(  . ~ comparison )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")
  #+
    #ggtitle("Growth ")
  ggsave("Mating_all_batches.pdf")
  
  
  
  
  Dataset_all_batches %>%
    filter( competition == 'Test') %>%
    #filter ( batch == 'batch1') %>%
    filter ( batch %in% c('batch2' ,'batch3')) %>%
    #filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = type ,y =log(relative_W_top))) +
    #theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    #scale_fill_brewer(palette="Set1")+
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( batch~comparison)+     
    #facet_grid(  . ~ comparison )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")
  #+
  #ggtitle("Growth ")
  ggsave("Top_selection_all_batches.pdf")
  
  Dataset_all_batches %>%
    filter( competition == 'Test') %>%
    #filter ( batch == 'batch1') %>%
    filter ( batch %in% c('batch2' ,'batch3')) %>%
    #filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = type ,y =log(relative_W_bottom))) +
    #theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    #scale_fill_brewer(palette="Set1")+
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( batch~comparison)+     
    #facet_grid(  . ~ comparison )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")
  #+
  #ggtitle("Growth ")
  ggsave("Bottom_selection_all_batches.pdf")
  ################# tukwysss
  growth_batch_2_anova<- 
    Dataset_all_batches %>%
    filter ( batch ==('batch3')) %>%
    filter( competition == 'Test') %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    aov(G_log~Gene_of_interest,.)
  growth2_tukey<-TukeyHSD(growth_batch_2_anova)
  growth2_df<-as.data.frame(growth2_tukey[1:1])
  write_csv(growth2_df, file="./RLC_rep2_byr2_Fitness/growth_batch3_tukey.csv")
  #######fullcycle#############
  ####batch2####
  Fullcycle_2_consol <- left_join( lookup,growth_batch2_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 
  Fullcycle_2_consol <- left_join( Fullcycle_2_consol,Fullcycle_top_batch2_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,m.evol,competition.x,comparison.x,type.x,s,S_T_log, S_T) %>%
  rename(Top_selection_differential=s)
  Fullcycle_2_consol <- left_join(  Fullcycle_2_consol,Fullcycle_bottom_batch2_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,competition.x,comparison.x,type.x,Top_selection_differential,S_T_log, S_T,s,S_B,S_B_log)%>%
    rename(Bottom_selection_differential=s )%>%
    rename(competition = competition.x)%>%
    rename(comparison = comparison.x)%>%
    rename(type = type.x)
  
  
  Fullcycle_2_final<- Fullcycle_2_consol %>%
    filter( competition == "Test" )
  Fullcycle_2_Dataset <- Fullcycle_2_final %>%
    add_column(batch = "batch2")
  #####Batch3####
  Fullcycle_3_consol <- left_join( lookup,growth_batch3_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,m.evol,G,G_log,competition,comparison,type) 
  Fullcycle_3_consol <- left_join( Fullcycle_3_consol,Fullcycle_top_batch3_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,m.evol,competition.x,comparison.x,type.x,s,S_T_log, S_T) %>%
    rename(Top_selection_differential=s)
  Fullcycle_3_consol <- left_join(  Fullcycle_3_consol,Fullcycle_bottom_batch3_normalised,by ="well")%>% 
    #filter( competition == "Test" )%>%
    select(genotype.x, well,competition.x,comparison.x,type.x,Top_selection_differential,S_T_log, S_T,s,S_B,S_B_log)%>%
    rename(Bottom_selection_differential=s )%>%
    rename(competition = competition.x)%>%
    rename(comparison = comparison.x)%>%
    rename(type = type.x)
  
  
  Fullcycle_3_final<- Fullcycle_3_consol %>%
    filter( competition == "Test" )
  Fullcycle_3_Dataset <- Fullcycle_3_final %>%
    add_column(batch = "batch3")
  
  Full_cycle_all_batches<-  rbind(Fullcycle_3_Dataset ,Fullcycle_2_Dataset)
  
  Full_cycle_all_batches <- Full_cycle_all_batches %>%
    mutate(Gene_of_interest = case_when(
      endsWith(genotype.x, "RwBf") ~ "rep2 wild type",
      endsWith(genotype.x, "RfBw") ~ "byr2 wild type",
      endsWith(genotype.x, "RfBm") ~ "byr2 mutant",
      endsWith(genotype.x, "RmBm") ~ "Double Mutant",
      endsWith(genotype.x, "A") ~ "Wild type",
      endsWith(genotype.x, "RmBf") ~ "rep2 mutant"
    ))
  
  
  
 ####plots#####
  Full_cycle_all_batches  %>%
    filter( competition == 'Test') %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = Gene_of_interest, y = S_B_log,fill= Gene_of_interest), palette = "jco") +
    geom_boxplot() +   
    theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    scale_fill_brewer(palette="BuPu")+ stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( .~ batch )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
    ggtitle("Full Cycle Bottom ")
  
  
  Full_cycle_all_batches  %>%
    filter( competition == 'Test') %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    ggplot( aes(x = Gene_of_interest, y = S_T_log,fill= Gene_of_interest), palette = "jco") +
    geom_boxplot() +   
    theme(legend.position="none") +
    #scale_color_manual(values= colorBlindGrey8 )+
    scale_fill_brewer(palette="BuPu")+ stat_summary(fun.data = give.n, geom = "text")+
    facet_grid( .~ batch )+  
    #facet_grid(  . ~ type )+ 
    stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
    ggtitle("Full Cycle Top ")
  
  
  full_cycle_top_total_anova<- 
    Full_cycle_all_batches  %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter( competition == 'Test') %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    aov(S_T_log~Gene_of_interest+batch,.)
  full_cycle_top_total_tukey<-TukeyHSD( full_cycle_top_total_anova)
  full_cycle_top_tukey_df<-as.data.frame( full_cycle_top_total_tukey[1:1])
  write_csv( full_cycle_top_tukey_df, file="./RLC_rep2_byr2_Fitness/ full_cycle_top_total_tukey.csv")
  
  full_cycle_bottom_total_anova<- 
    Full_cycle_all_batches  %>%
    #filter ( batch %in% c('batch2' ,'batch3')) %>%
    filter( competition == 'Test') %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    aov(S_B_log~Gene_of_interest+batch,.)
  full_cycle_bottom_total_tukey<-TukeyHSD( full_cycle_bottom_total_anova)
  full_cycle_bottom_tukey_df<-as.data.frame( full_cycle_bottom_total_tukey[1:1])
  write_csv( full_cycle_bottom_tukey_df, file="./RLC_rep2_byr2_Fitness/ full_cycle_bottom_total_tukey.csv")
  
  
  full_cycle_top_batch2_anova<- 
    Full_cycle_all_batches  %>%
    filter ( batch =='batch2' ) %>%
    filter( competition == 'Test') %>%
    filter(genotype.x %in% c("A","RfBm","RmBf","RmBm")  )%>%
    aov(S_T_log~Gene_of_interest,.)
  full_cycle_top_batch2_tukey<-TukeyHSD( full_cycle_top_batch2_anova)
  full_cycle_top_tukey_df_2<-as.data.frame( full_cycle_top_batch2_tukey[1:1])
  write_csv( full_cycle_top_tukey_df_2, file="./RLC_rep2_byr2_Fitness/ full_cycle_top_batch2_tukey.csv")
  