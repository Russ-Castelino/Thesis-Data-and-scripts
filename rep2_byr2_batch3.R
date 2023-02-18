####install packages######
#install.packages("ggpubr","ggbiplot")
#install_github("vqv/ggbiplot")


library( tidyverse )
library( gridExtra )
library( patchwork )
library( ggpubr )
library(devtools)
library(cluster)
library(ggfortify)
library(plotly)
#install.packages( c('tidyverse','gridExtra', 'patchwork' ) )

# library(flowCore)
library( sf )

library(tidyverse)
library(gridExtra)
library(ggpmisc)
library(ggplot2)
library(BBmisc)
#####Attempt witj Batch 2######
setwd("D:/Wolf project/Data processing")

dir <-"./RLC_rep2_byr2_Fitness/Batch_3"
give.n <- function(x){
  return(c(y = median(x), label = length(x)))}
# ######rename files so they have the same pattern
filesafter <- list.files( path = dir, pattern = "bottom_AG_230120_Batch3_undiluted_*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("bottom_","B",f) )
}
filesafter <- list.files( path = dir, pattern = "TopAG_230120_Batch3_undiluted_*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("TopAG","TAG",f) )
}
filesafter <- list.files( path = dir, pattern = "230119_Before_growth*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("230119_Before_growth","BG_230119_batch3",f) )
}
filesafter <- list.files( path = dir, pattern = "MatingAG_230126_batch3_undiluted_*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("MatingAG_230126_batch3_undiluted","MAG_230126_batch3_undiluted",f) )
}
filesafter <- list.files( path = dir, pattern = "Full_Cycle_Bottom_230126_batch3_undiluted_*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("Full_Cycle_Bottom_230126_batch3_undiluted_","FCB_230126_batch3_undiluted",f) )
}
filesafter <- list.files( path = dir, pattern = "Full _Cycle_Top  230126_batch3_diluted_*.*.fcs$", full.names = T, recursive = T)
for(f in filesafter){
  file.rename(f, gsub("Full _Cycle_Top  230126_batch3_diluted","FCT_230126_batch3_diluted",f) )
}
# this should be a file that contains the information about which file is what
#cor_dfbatch3 <- read_csv("./input/all_files_before.csv", col_types = "cddd")


source("./scripts/functions_fcs.R")

# get fcs file names
fcsFiles3 <- list.files( path = dir, pattern = "*.*.fcs$", full.names = T, recursive = T)
fcsFiles3

# We also need to initiate a dataframe to store all data in 
allMeasures_batch3 <- tibble(sample = NA, plate = NA,
                             totalCounts = NA,
                             blue = NA, wt = NA, all = NA)[0,] # the [0,] is used to remove all data


draw =    TRUE #    FALSE  #

f = fcsFiles3[1]

#### Start loop ####
for(f in fcsFiles3 ){ # start the loop #csvFiles[ !str_starts(basename(csvFiles),"_" ) ]
  
  #### load data ####
  dfbatch3 <- getFCS( f ) # read the data
  
  print(basename(f))  
  allCounts <- nrow(dfbatch3)
  
  dfbatch3 <- dfbatch3 %>%
    mutate(GFP =    log(FITC_A + abs(min(FITC_A)) + 1, 10),
           mCherry = log(PE_Texas_Red_A + abs(min(PE_Texas_Red_A)) + 1, 10),
           bfp = log( BV421_A + abs(min( BV421_A )) + 1, 10)) %>%
    filter(FSC_A < 240000, SSC_A < 240000)
  
  
  #### Gating for cells #### 
  mF = 100000
  cs = 50000 # cells_size
  
  # make the gate
  FSC_SCC <-  matrix( 
    c( 50000  , 100000, 150000, 200000,200000,100000,50000, 50000   ,# FSC_A
       150000 , 150000, 120000, 90000, 30000, 30000, 50000, 150000   # SSC_A
    ), ncol = 2  )
  
  # make a plot
  p_all <- plot_oneGate( dfbatch3, "FSC_A", "SSC_A", FSC_SCC )
  p_all
  
  # add the filter variable 'cells'
  dfbatch3 <- dfbatch3 %>% mutate( cells = gateCells_2d(., xCh = "FSC_A", yCh = "SSC_A", FSC_SCC) )
  
  #### Gating for singleton ####
  # make the gate
  FSCH_FSCW <-  matrix( 
    c( 70000, 90000, 135000, 135000,   90000,   75000,  60000, 70000            ,# FSC_H
       60000, 35000, 25000,  45000,    80000, 185000, 185000, 60000             # FSC_W
    ), ncol = 2  )
  
  # make a plot
  p_cells <- dfbatch3 %>%
    filter(cells) %>% 
    plot_oneGate( "FSC_H", "FSC_W", FSCH_FSCW, 
                  setX = c( 60000, 140000), 
                  setY = c( 20000, 185000) )
  p_cells
  
  # add the filter variable 'single_cells'
  dfbatch3 <- dfbatch3 %>% mutate( single_cells = gateCells_2d(., xCh = "FSC_H", yCh = "FSC_W", FSCH_FSCW) )
  
  
  #### BFP ####
  cval = 0
  
  bfp_FSC_neg<-  matrix( 
    c(
      c( mF - cs, mF + cs*1.2, mF + cs*1.2, mF - cs, mF - cs ),   # FSC_A
      c(  0    , 0     , 3.8     , 3.8    , 0     )  + cval   # bfp
    ), ncol = 2  )
  
  bfp_FSC_pos<-  matrix( 
    c( 
      c( mF - cs, mF + cs*1.2, mF + cs*1.2, mF - cs, mF - cs ),   # FSC_A
      c( 5.5    , 5.5     , 3.85     , 3.85     , 5.5      )  + cval   # bfp
    ), ncol = 2  )
  
  dfbatch3 <- dfbatch3 %>% mutate( bfpPOS = gateCells_2d(., xCh = "FSC_A", yCh = "bfp", bfp_FSC_pos) )
  dfbatch3 <- dfbatch3 %>% mutate( bfpNEG = gateCells_2d(., xCh = "FSC_A", yCh = "bfp", bfp_FSC_neg) )
  
  
  # plot_noGate(dfbatch3, y = "bfp", x = "FSC_A",
  #             setX = c( mF - 1.2 * cs, mF + cs * 1.2),
  #             setY = c( 1 + cval, 5.5 + cval) )
  
  p_bfp <- dfbatch3 %>% 
    filter(single_cells & cells ) %>%
    plot_twoGates(., "FSC_A", "bfp", bfp_FSC_neg, bfp_FSC_pos,
                  setX = c( mF - 1.2 * cs, mF + cs * 1.2),
                  setY = c( 0 + cval, 5.5 + cval) )
  
  p_bfp
  
  plots <- p_all + p_cells +  p_bfp 
  
  #plot(dens)
  if(draw) { 
    dir.create(file.path(dirname(f),"/figures_batch3/"), showWarnings = FALSE)
    ggsave(filename = paste0(dirname(f),"/figures_batch3/",basename(f),".jpg"),
           plot = plots,
           width = 12, height = 6, units = "cm"
    )
    
  }
  
  blue <-  sum( dfbatch3$bfpPOS & dfbatch3$cells & dfbatch3$single_cells )
  wt <-  sum(   dfbatch3$bfpNEG & dfbatch3$cells & dfbatch3$single_cells )
  all <- blue + wt
  
  ## add the measurements to the dataframe allMeasurements
  allMeasures_batch3 <- rbind(allMeasures_batch3, # that's what we already had
                              tibble( # here we make a 'table of one row to add
                                sample = basename( f ),
                                plate = 1,
                                totalCounts = allCounts,
                                blue = blue,
                                wt = wt,
                                all = all
                              ))
  
}
#### end loop ####





###Data wrangling####

dfbatch3 <- allMeasures_batch3 %>%
  separate( sample, into = c("timepoint","date","batch","well1","well","rest"), sep = "_") %>%
  select(-rest, -well1) %>%
  mutate( row = str_sub( well, 1,1),
          col = as.numeric( str_sub( well, 2) )  )
dfbatch3

dir.create(file.path( "Databatch3" ), showWarnings = FALSE)
write_csv(dfbatch3, file = "./Databatch3/230119_data_batch3.csv")
#dfbatch3<- read.csv(file = 'D:/Wolf project/Data processing/Databatch3/221202_data_batch3.csv')

#upload the order of the plates 

#df0 <- read.csv(file = 'D:/Wolf project/Data processing/plate_order.csv')

#rep_byr_batch3 <- merge(dfbatch3, df0, by = 'well')
#write_csv(rep_byr_batch3, file = "./Databatch3/221202_data_batch3_order.csv")
#rep_byr_batch3 <- read.csv(file = "D:/Wolf project/Data processing/Databatch3/221202_data_batch3_order.csv")
#rep_byr_Batch1_BG <-filter(rep_byr_Batch1, timepoint == "BG")
#rep_byr_Batch1_AG <-filter(rep_byr_Batch1, timepoint == "AG")

#write_csv(dfbatch3, file = "./Data/2212_data_batch3.csv")

#dfbatch3 <- read_csv("./Data/221202_data_batch1.csv")

lookup <- read.csv(file = 'D:/Wolf project/Data processing/plate_order.csv')






######changing somethings#####
dfbatch3$well <- str_replace(dfbatch3$well, "A10", "AX")
dfbatch3$well <- str_replace(dfbatch3$well, "B10","BX")
dfbatch3$well <- str_replace(dfbatch3$well, "C10","CX")
dfbatch3$well <- str_replace(dfbatch3$well, "D10","DX")
dfbatch3$well <- str_replace(dfbatch3$well, "E10","EX")
dfbatch3$well <- str_replace(dfbatch3$well, "F10","FX")
dfbatch3$well <- str_replace(dfbatch3$well, "G10","GX")
dfbatch3$well <- str_replace(dfbatch3$well, "H10","HX")
####

dfbatch3$well <- str_replace(dfbatch3$well, "A1","A01")
dfbatch3$well <- str_replace(dfbatch3$well, "A2","A02")
dfbatch3$well <- str_replace(dfbatch3$well, "A3","A03")
dfbatch3$well <- str_replace(dfbatch3$well, "A4","A04")
dfbatch3$well <- str_replace(dfbatch3$well, "A5","A05")
dfbatch3$well <- str_replace(dfbatch3$well, "A6","A06")
dfbatch3$well <- str_replace(dfbatch3$well, "A7","A07")
dfbatch3$well <- str_replace(dfbatch3$well, "A8","A08")
dfbatch3$well <- str_replace(dfbatch3$well, "A9","A09")
dfbatch3$well <- str_replace(dfbatch3$well, "B1","B01")
dfbatch3$well <- str_replace(dfbatch3$well, "B2","B02")
dfbatch3$well <- str_replace(dfbatch3$well, "B3","B03")
dfbatch3$well <- str_replace(dfbatch3$well, "B4","B04")
dfbatch3$well <- str_replace(dfbatch3$well, "B5","B05")
dfbatch3$well <- str_replace(dfbatch3$well, "B6","B06")
dfbatch3$well <- str_replace(dfbatch3$well, "B7","B07")
dfbatch3$well <- str_replace(dfbatch3$well, "B8","B08")
dfbatch3$well <- str_replace(dfbatch3$well, "B9","B09")
dfbatch3$well <- str_replace(dfbatch3$well, "C1","C01")
dfbatch3$well <- str_replace(dfbatch3$well, "C2","C02")
dfbatch3$well <- str_replace(dfbatch3$well, "C3","C03")
dfbatch3$well <- str_replace(dfbatch3$well, "C4","C04")
dfbatch3$well <- str_replace(dfbatch3$well, "C5","C05")
dfbatch3$well <- str_replace(dfbatch3$well, "C6","C06")
dfbatch3$well <- str_replace(dfbatch3$well, "C7","C07")
dfbatch3$well <- str_replace(dfbatch3$well, "C8","C08")
dfbatch3$well <- str_replace(dfbatch3$well, "C9","C09")
dfbatch3$well <- str_replace(dfbatch3$well, "D1","D01")
dfbatch3$well <- str_replace(dfbatch3$well, "D2","D02")
dfbatch3$well <- str_replace(dfbatch3$well, "D3","D03")
dfbatch3$well <- str_replace(dfbatch3$well, "D4","D04")
dfbatch3$well <- str_replace(dfbatch3$well, "D5","D05")
dfbatch3$well <- str_replace(dfbatch3$well, "D6","D06")
dfbatch3$well <- str_replace(dfbatch3$well, "D7","D07")
dfbatch3$well <- str_replace(dfbatch3$well, "D8","D08")
dfbatch3$well <- str_replace(dfbatch3$well, "D9","D09")
dfbatch3$well <- str_replace(dfbatch3$well, "E1","E01")
dfbatch3$well <- str_replace(dfbatch3$well, "E2","E02")
dfbatch3$well <- str_replace(dfbatch3$well, "E3","E03")
dfbatch3$well <- str_replace(dfbatch3$well, "E4","E04")
dfbatch3$well <- str_replace(dfbatch3$well, "E5","E05")
dfbatch3$well <- str_replace(dfbatch3$well, "E6","E06")
dfbatch3$well <- str_replace(dfbatch3$well, "E7","E07")
dfbatch3$well <- str_replace(dfbatch3$well, "E8","E08")
dfbatch3$well <- str_replace(dfbatch3$well, "E9","E09")
dfbatch3$well <- str_replace(dfbatch3$well, "F1","F01")
dfbatch3$well <- str_replace(dfbatch3$well, "F2","F02")
dfbatch3$well <- str_replace(dfbatch3$well, "F3","F03")
dfbatch3$well <- str_replace(dfbatch3$well, "F4","F04")
dfbatch3$well <- str_replace(dfbatch3$well, "F5","F05")
dfbatch3$well <- str_replace(dfbatch3$well, "F6","F06")
dfbatch3$well <- str_replace(dfbatch3$well, "F7","F07")
dfbatch3$well <- str_replace(dfbatch3$well, "F8","F08")
dfbatch3$well <- str_replace(dfbatch3$well, "F9","F09")
dfbatch3$well <- str_replace(dfbatch3$well, "G1","G01")
dfbatch3$well <- str_replace(dfbatch3$well, "G2","G02")
dfbatch3$well <- str_replace(dfbatch3$well, "G3","G03")
dfbatch3$well <- str_replace(dfbatch3$well, "G4","G04")
dfbatch3$well <- str_replace(dfbatch3$well, "G5","G05")
dfbatch3$well <- str_replace(dfbatch3$well, "G6","G06")
dfbatch3$well <- str_replace(dfbatch3$well, "G7","G07")
dfbatch3$well <- str_replace(dfbatch3$well, "G8","G08")
dfbatch3$well <- str_replace(dfbatch3$well, "G9","G09")
dfbatch3$well <- str_replace(dfbatch3$well, "H1","H01")
dfbatch3$well <- str_replace(dfbatch3$well, "H2","H02")
dfbatch3$well <- str_replace(dfbatch3$well, "H3","H03")
dfbatch3$well <- str_replace(dfbatch3$well, "H4","H04")
dfbatch3$well <- str_replace(dfbatch3$well, "H5","H05")
dfbatch3$well <- str_replace(dfbatch3$well, "H6","H06")
dfbatch3$well <- str_replace(dfbatch3$well, "H7","H07")
dfbatch3$well <- str_replace(dfbatch3$well, "H8","H08")
dfbatch3$well <- str_replace(dfbatch3$well, "H9","H09")
###
dfbatch3$well <- str_replace(dfbatch3$well, "AX", "A10")
dfbatch3$well <- str_replace(dfbatch3$well, "BX","B10")
dfbatch3$well <- str_replace(dfbatch3$well, "CX","C10")
dfbatch3$well <- str_replace(dfbatch3$well, "DX","D10")
dfbatch3$well <- str_replace(dfbatch3$well, "EX","E10")
dfbatch3$well <- str_replace(dfbatch3$well, "FX","F10")
dfbatch3$well <- str_replace(dfbatch3$well, "GX","G10")
dfbatch3$well <- str_replace(dfbatch3$well, "HX","H10")




######make a joined matrix that combines the before to each correct after ####
dfbatch3j <- dfbatch3  %>%  filter( timepoint != "BG" ) %>% select( well, plate, blue, wt, all, timepoint ) %>%
  left_join(., 
            dfbatch3 %>% filter( timepoint == "BG" ) %>% select(-timepoint, -totalCounts), 
            by = c( "well" , "plate" ), 
            suffix = c(".after",".before") ) 






##### add the information what is in which well####
dfbatch3j <- left_join(dfbatch3j, lookup, by ="well") %>% 
  select(-rep2_allele, -byr2_allele, starts_with("all")) %>%
  select(-date, -batch, -row, -col)
dfbatch3j$competition <- NA
dfbatch3j<-dfbatch3j %>%
  mutate(competition = case_when(
    endsWith(genotype, "(control)") ~ "Control",
    endsWith(genotype, "(Control)") ~ "Control",
    endsWith(genotype, "") ~ "Test"
  ))
dfbatch3j <- dfbatch3j %>% mutate(type = case_when(
  endsWith(genotype, "Bm") ~ "Mutant",
  endsWith(genotype, "RmBf") ~ "Mutant",
  endsWith(genotype, "Bw") ~ "Wild type",
  endsWith(genotype, "RwBf") ~ "Wild type",
  endsWith(genotype, "A") ~ "Wild type"
))
dfbatch3j <- dfbatch3j %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))




#####Competitions#####
## this will make a very simple plot of just the fitness of each, without any compensations for growth
summarise(dfbatch3js)
dfbatch3js<-
  dfbatch3j %>% 
  mutate( s = log(wt.after/blue.after) - log(wt.before/blue.before)) %>%
  filter( competition == "Test" )%>%
  filter( competition == "Test" )%>%
  filter( all.before >500 )%>%
  filter( all.after >150 ) 
write_csv(dfbatch3j, file = "./Databatch3/rep2_byr2_main_batch3.csv")
write_csv(dfbatch3js, file = "./Databatch3/rep2_byr2_analysis_batch3.csv")
dfbatch3js%>%
  ggplot(aes( x = type, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ comparison, scales = "free_y")+  
  stat_compare_means(method = "wilcox.test",label.x=2.2, label.y=0,label="p.signif")+
  ggtitle("Competitions (Batch 2)")
ggsave("batch3_competitions_byr2_rep2.pdf")
#ggsave("batch3_competitions_byr2_rep2.pdf")
#dfbatch3j<-read.csv(file = 'D:/Wolf project/Data processing/Databatch3/rep2_byr2_main_batch3.csv')
#dfbatch3js<-read.csv(file = 'D:/Wolf project/Data processing/Databatch3/rep2_byr2_analysis_batch3.csv')





######## Growth - fitness calculations ##########
### the main idea of the next part is that we want to calculate the malthusian parameter
### for each evolved line, using the absolute values (number of cells) along the way.
### Using the known values for the reference we can do these calculations, and then
### we can use these later to compensate for growth in the selection and mating measurements.
growth_batch3 <- 
  dfbatch3js %>% 
  filter(timepoint == "AG") %>% # only look at growth
  mutate(# define total cells when saturated
    Sat_ref_total=30000000, 
    # calculate the number we started growth with for blue after all the dilutiones
    N_ref_BG=((((Sat_ref_total/300)*100)/350)*5)*30/505 ,
    # calculate the number of cells for not-blue
    
    N_evo_BG=( (wt.before/all.before) * N_ref_BG)/( blue.before / all.before ) , 
    # define the number of cell cycles that occurred during growth
    tAG = 6 ,
    ###m.ref is assumed to be 1 for N_AG = N_BG*exp(tAG*m)
    # calculate how many cells there were after growth for the ref
    N_ref_AG = N_ref_BG * exp( tAG ) ,
    # calculate the fraction of evolved cells after growth
    f.evol.AG = wt.after / all.after ,
    # and from that we can calculate the number of evolved cells after growth
    N_evol_AG =( f.evol.AG * N_ref_AG ) / ( 1 - f.evol.AG ), 
    # and now that we know the before and the after and the number of generations,
    # we can use these to calculate the Malthusian parameter from the 
    # equations N_t = N_0 * e^(m * t)
    m.evol = log( N_evol_AG / N_evo_BG ) / tAG ,
    m.evol_log = log(m.evol))


growth_batch3 %>%
  select( - starts_with( "wt." ),
          - starts_with( "blue."),
          - starts_with( "all." ) )

#G_batch3 <- 
growth_batch3 %>%
  filter( competition == 'Test') %>%
  ggplot( aes(x = type, y = m.evol_log), palette = "jco") +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+ 
  facet_grid(  . ~ comparison )+ 
  stat_compare_means(method = "t.test",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Growth (batch3)")
ggsave("batch3_growth_byr2_rep2.pdf")
ggpubr::compare_means(m.evol_log~genotype,  growth_batch3 , method = "t.test", paired = FALSE,
                      group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
growth_batch3_anova<-aov(m.evol_log~genotype,data=growth_batch3 )
TukeyHSD(growth_batch3_anova)
#G_batch3 + stat_compare_means(method = "wilcox.test")
ggsave("batch3_growth_byr2_rep2.pdf")

growth_batch3 %>%
  filter( competition == 'Test') %>%
  #filter( comparison %in% c( 'rep2',"byr2")) %>%
  ggplot( aes(x = comparison, y = m.evol), palette = "jco") +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+ 
  facet_grid( type ~ . )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=1,label="p.signif")+
  ggtitle("Growth (batch3)")
ggsave("batch3_growth_mutant_byr2_rep2.pdf")
growth_batch3 %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  filter( competition == 'Test') %>%
  #filter( comparison %in% c( 'rep2',"byr2")) %>%
  ggplot( aes(x = genotype, y = m.evol), palette = "jco") +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+ 
  #facet_grid( type ~ . )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=1,label="p.signif")+
  ggtitle("Growth (batch3)")
####significance growth###


growth_batch3_pairwise<-ggpubr::compare_means(m.evol_log~genotype,  growth_batch3 , method = "t.test", paired = FALSE,
                      group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(growth_batch3_pairwise, file = "./Databatch3/growth_pairwise-t-test_batch3.csv")
#####eliminate the effect of fluorophores###
growth.median.values.batch3<- 
  growth_batch3%>% 
  filter(type=="Wild type") %>% 
  filter( competition == 'Test') %>%
  group_by(genotype) %>%
  dplyr::summarize( 
    m.evol.median.wild = median( m.evol))
#add the comparison
growth.median.values.batch3 <- growth.median.values.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


growth.median.values.batch3<-select(growth.median.values.batch3, -genotype)
growth_batch3_normalised<- left_join(growth_batch3,growth.median.values.batch3,by="comparison")
growth_batch3_normalised$G <- growth_batch3_normalised$m.evol/growth_batch3_normalised$m.evol.median.wild
growth_batch3_normalised$G_log<-log(growth_batch3_normalised$G)
###box plot###




growth_batch3_normalised %>%
  filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y = G_log), palette = "jco") +
  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +    
  #facet_grid(  . ~ type )+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Growth (batch3)")
ggsave("batch3_growth_byr2_rep2_adjusted.pdf")


##signficance growth##
growth_batch3_pairwise <-growth_batch3_normalised %>%
filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
 ggpubr::compare_means(G_log~genotype,  . , method = "t.test", paired = FALSE,
                                              group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(growth_batch3_pairwise, file = "./Databatch3/growth_adjusted_pairwise-t-test_batch3.csv")


## calculate the median, mean and sd values for each genotype
growth.summary.values_batch3 <- 
  growth_batch3 %>%
  group_by(genotype) %>%
  dplyr::summarize( 
    median.m = median( m.evol ),
    mean.m = mean( m.evol ),
    sd.m = sd( m.evol ))
write_csv(growth_batch3, file = "./Databatch3/Growth_adjusted_batch3.csv")
write_csv(growth.summary.values_batch3, file = "./Databatch3/Growth_adjusted_summary_batch3.csv")
#summary(growth_batch3)
#growth normalising
#get median values of wild type
#growth.median.values.batch3 <-  growth_batch3 %>% 
 # filter(type=="Wild type") %>% 
  #group_by(genotype) %>%
  #dplyr::summarize( 
   # median.m.wild = median( m.evol ))
#add the comparison
#growth.median.values.batch3 <- growth.median.values.batch3 %>% mutate(comparison = case_when(
 # endsWith(genotype, "Bf") ~ "rep2",
#  endsWith(genotype, "RfBw") ~ "byr2",
 # endsWith(genotype, "RfBm") ~ "byr2",
#  endsWith(genotype, "RmBm") ~ "Double Mutant",
 # endsWith(genotype, "A") ~ "Double Mutant"
#))

#growth.median.values.batch3<-select(growth.median.values.batch3, -genotype)
#growth_batch3_normalised<- left_join(growth_batch3,growth.median.values.batch3,by="comparison")%>% 
 # filter(type=="Mutant") 
#growth_batch3_normalised$m.evol_divided <- growth_batch3_normalised$m.evol/growth_batch3_normalised$median.m.wild
#growth_batch3_normalised$m.evol_divided_log<-log(growth_batch3_normalised$m.evol_divided)

#growth_batch3_normalised %>%
 # filter( competition == 'Test') %>%
  #filter( comparison %in% c( 'rep2',"byr2")) %>%
  #ggplot( aes(x = comparison, y = m.evol_divided_log), palette = "jco") +
  #  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text")+ 
  #facet_grid( type ~ . )+ 
  #stat_compare_means(method = "anova",label.x=0.5, label.y=-0.03,label="p.signif")+
  #ggtitle("Growth (batch3)")










#### Selection - fitness calculations ####
### similar to how we calculated growth, we will again use absolute numbers to 
### back calculate the values directly after selection. We use the Malthusian parameter
### calculated above to compensate for growth.
ecological_sel_batch3 <-
  dfbatch3js  %>% 
  # keep only the data from top and bottom selection
  filter( timepoint %in% c("BAG", "TAG") ) %>% 
  # and only the competitions
  filter( competition == "Test"  ) %>%
  # now join for each of these the median value for the Malthusian parameter we calculated above
  left_join( . , 
             growth.summary.values_batch3, by=c( "genotype" ) ) %>% 
  # and now again, we will do the calculations for each step
  mutate(
    # Set absolute numbers before
    Sat_ref_total=30000000, 
    # set growth cycles
    tAG=6 ,
    # set expected value for the reference 
    Expected_ref_BS = Sat_ref_total / 6, 
    # fraction of cells that should be selected. This number is assumed,
    # as we measured this for the wildtype before. Doesn't matter too much if not 100% correct.
    selection_co_0 = 0.01, 
    # absolute number of reference cells that are selected (i.e. 1%). 
    Expected_ref_After_Sel = Sat_ref_total * selection_co_0, 
    # incorporate al the dilutions we did to get to absolute number of cells
    N_ref_BS = (Sat_ref_total*100/300)*50/350, 
    # calculate the number of cells for wt strains before selection
    N_evo_BS = ( ( blue.before / all.before)  * N_ref_BS )/ ( wt.before / all.before ), 
    # calculate the number of reference cells after selection
    N_ref_AS = N_ref_BS * Expected_ref_After_Sel / Expected_ref_BS,
    # Now we can calculate the number of cells that were selected after growth,  
    # because we know the start number, and the fraction that were selection and from 
    # that can now calculate the number that we should have after 6 asexual growth cycles. 
    N_ref_S=N_ref_AS*exp(tAG) , 
    # and now calculate the fraction of wild type strains
    f.evol.S= wt.after / all.after , 
    # the fraction can then be used to calculate the absolute number after selection and growth
    N_evo_S = ( f.evol.S * N_ref_S ) / ( 1 - f.evol.S ), 
    # and then we can calculate back the number after selection, before growth
    N_evo_AS = N_evo_S / (exp( tAG * median.m ) ), 
    # calculate the ratio in frequencies over the selection step
    v_ref = N_ref_AS / N_ref_BS, 
    v_evo = N_evo_AS / N_evo_BS, 
    # and from there for the evolved strains the relative fitness
    relative_W_evo = v_evo / v_ref, 
    # from which you can calculate the selection coefficients
    s_evo = 1 - relative_W_evo,
    # As this is not the most direct way to calculate fitness, we also calculate using
    # the frequencies and taking the log of their ratios:
    # calculate the frequencies at the different time points
    fq_BS_evo = N_evo_BS / ( N_evo_BS + N_ref_BS ),
    fq_AS_evo = N_evo_AS / ( N_evo_AS + N_ref_AS ), 
    fq_BS_ref = N_ref_BS / ( N_evo_BS + N_ref_BS ), 
    fq_AS_ref = N_ref_AS / ( N_evo_AS + N_ref_AS ),
    # calculate the selection coefficient estimates
    s_evol_difflog = log( ( fq_AS_evo / ( 1 - fq_AS_evo ) ) ) - log( ( fq_BS_evo / ( 1 - fq_BS_evo) ) )
  )

summary(ecological_sel_batch3$genotype)
ecological_sel_batch3 %>% 
  filter(competition=="Test")%>% 
  filter(timepoint=="TAG") %>% 
  group_by(genotype) %>% 
  count()



ecological_sel_batch3 %>%
  filter( competition == 'Test') %>%
  filter(timepoint=="TAG")%>%
  ggplot( aes(x = type, y =  log(relative_W_evo) ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ comparison ) +  
  ggpubr::stat_compare_means(method = "t.test",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Ecological Selection (batch3)")
ggsave("batch3_eco_selection_top_byr2_rep2.pdf")
write_csv(ecological_sel_batch3, file = "./Databatch3/eco_selection_adjusted_batch3.csv")


ecological_sel_bottom_batch3<- ecological_sel_batch3 %>%
  filter( competition == 'Test') %>%
  filter( timepoint == 'BAG')
ecological_sel_top_batch3<- ecological_sel_batch3 %>%
  filter( competition == 'Test') %>%
  filter( timepoint == 'TAG')



####t.test##
ecological_sel_bottom_batch3_pairwise<-ggpubr::compare_means(relative_W_evo~genotype,  ecological_sel_bottom_batch3 , method = "t.test", paired = FALSE,
                                              group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(ecological_sel_bottom_batch3_pairwise, file = "./Databatch3/ecological_sel_bottom_pairwise-t-test_batch3.csv")

ecological_sel_top_batch3_pairwise<-ggpubr::compare_means(relative_W_evo~genotype,  ecological_sel_top_batch3 , method = "t.test", paired = FALSE,
                                                             group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(ecological_sel_top_batch3_pairwise, file = "./Databatch3/ecological_sel_top_pairwise-t-test_batch3.csv")
##normalising selection
#eco selectiom
#bottom
#get median values of wild type

eco.median.values.bottom.batch3 <- 
  ecological_sel_bottom_batch3 %>% 
  filter(type=="Wild type") %>% 
  group_by(genotype) %>%
  dplyr::summarize( 
    WB.median.wild = median( relative_W_evo ))
#add the comparison
eco.median.values.bottom.batch3 <- eco.median.values.bottom.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


eco.median.values.bottom.batch3<-select(eco.median.values.bottom.batch3, -genotype)
ecological_sel_bottom_batch3_normalised<- left_join(ecological_sel_bottom_batch3,eco.median.values.bottom.batch3,by="comparison")
ecological_sel_bottom_batch3_normalised$ESB <- ecological_sel_bottom_batch3_normalised$relative_W_evo/ecological_sel_bottom_batch3_normalised$WB.median.wild
ecological_sel_bottom_batch3_normalised$ESB_log<-log(ecological_sel_bottom_batch3_normalised$ESB)
write_csv(ecological_sel_bottom_batch3_normalised, file = "./Databatch3/ecological_sel_bottom_batch3_normalised.csv")
###box plot####
ecological_sel_bottom_batch3_normalised%>%
  #filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y =  ESB_log ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ timepoint ) +  
  ggpubr::stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Ecological Selection (batch3)")

ecological_sel_bottom_batch3_adjusted_pairwise<-ecological_sel_bottom_batch3_normalised%>%
filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
ggpubr::compare_means(ESB_log~genotype,  ., method = "t.test", paired = FALSE,
                      group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(ecological_sel_bottom_batch3_adjusted_pairwise, file = "./Databatch3/ecological_sel_bottom_adjusted_pairwise-t-test_batch3.csv")


#Top

eco.median.values.top.batch3 <- 
  ecological_sel_top_batch3 %>% 
  filter(type=="Wild type") %>% 
  group_by(genotype) %>%
  dplyr::summarize( 
    WT.median.wild = median( relative_W_evo))
#add the comparison
eco.median.values.top.batch3 <- eco.median.values.top.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


eco.median.values.top.batch3<-select(eco.median.values.top.batch3, -genotype)
ecological_sel_top_batch3_normalised<- left_join(ecological_sel_top_batch3,eco.median.values.top.batch3,
                                                 by="comparison")
ecological_sel_top_batch3_normalised$EST<- ecological_sel_top_batch3_normalised$relative_W_evo/ecological_sel_top_batch3_normalised$WT.median.wild
ecological_sel_top_batch3_normalised$EST_log<-log(ecological_sel_top_batch3_normalised$EST)
write_csv(ecological_sel_top_batch3_normalised, file = "./Databatch3/ecological_sel_top_batch3_normalised.csv")
###box plot####
ecological_sel_top_batch3_normalised%>%
  #filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y =  EST_log ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ timepoint ) +  
  ggpubr::stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Ecological Selection (batch3)")
ggsave("batch3_eco_selection_top_byr2_rep2_adjusted.pdf")
ecological_sel_top_batch3_adjusted_pairwise<- ecological_sel_top_batch3_normalised%>%
filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
ggpubr::compare_means(EST_log~genotype,  ., method = "t.test", paired = FALSE,
                      group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(ecological_sel_top_batch3_adjusted_pairwise, file = "./Databatch3/ecological_sel_top_adjusted_pairwise-t-test_batch3.csv")


#ecological_sel_batch3 %>%
#filter( competition == 'Test') %>%
#filter( comparison %in% c( 'rep2',"byr2")) %>%
#ggplot( aes(x = comparison, y = s_evol_difflog ) ) +
#geom_boxplot() +
#facet_grid( type ~ timepoint ) +  
#stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")+
#ggtitle("Ecological Selection (batch3)")
#ggsave("batch3_eco_selection_mutation_byr2_rep2.pdf") 

#comparing mutants
#ecological_sel_batch3 %>%
# filter( competition == 'Test') %>%
#filter( type == "Mutant") %>%
#filter(comparison==c("rep2","byr2")) %>%
#ggplot( aes(x = genotype, y = s_evol_difflog ) ) +
#  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
#facet_grid(  ~ timepoint ) +  
#stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
###### mating######
# mating_data_ed<-rbind(mating_data_evolved_pop, mating_data_parental_pop) %>% 
mating_unedit_batch3<-  dfbatch3js %>% 
  filter( timepoint %in% c("MAG") ) %>% 
  # now join for each of these the median value for the Malthusian parameter we calculated above
  left_join( . , 
             growth.summary.values_batch3, by=c( "genotype" ) ) %>% 
  mutate(saturate_ref=30000000, 
         N_ref_BM=(((saturate_ref/300)*20)*10/100)*25/1000, 
         N_evol_BM=((blue.before/all.before)*N_ref_BM)/( wt.before / all.before ), 
         N_cell_AM=100000, #it does not matter which number is used here
         N_evo_AMG=( blue.after / all.after)*N_cell_AM, 
         N_ref_AMG=N_cell_AM-N_evo_AMG, 
         growth_time_generations=(21/3.5)-1, 
         N_evol_AM=N_evo_AMG/(exp(median.m*growth_time_generations)), 
         N_ref_AM=N_ref_AMG/(exp(growth_time_generations)), 
         r.AM=N_evol_AM/(N_evol_AM+N_ref_AM), 
         relative_v=(r.AM*N_ref_BM)/(N_evol_BM*(1-r.AM)))
#mating_unedit_batch3%>%
#filter( competition == 'Test') %>%
#filter(comparison=="rep2")

mating_unedit_batch3 %>%
  filter( competition == 'Test') %>%
  ggplot( aes(x = type, y = log(relative_v))) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( comparison ~ . )+
  stat_compare_means(method = "t.test",label.x=0.5, label.y=0,label="p.signif",)+
  ggtitle("Mating (Batch 3)")
ggsave("batch3_mating_byr2_rep2.pdf")
write_csv(mating_unedit_batch3, file = "./Databatch3/mating_batch3.csv")


mating_unedit_batch3_pairwise<- mating_unedit_batch3%>%
  #filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggpubr::compare_means(relative_v~genotype,  ., method = "t.test", paired = FALSE,
                        group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(mating_unedit_batch3_pairwise, file = "./Databatch3/mating_batch3_pairwise-t-test.csv")
#mating_unedit_batch3 %>%
 # filter( competition == 'Test') %>%
  #filter( comparison %in% c( 'rep2',"byr2")) %>%
  #ggplot( aes(x = comparison, y = relative_v)) +
  #  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ . )+
  #stat_compare_means(method = "anova",label.x=0.5, label.y=1,label="p.signif")+
  #ggtitle("Mating (Batch 2)")
#ggsave("batch3_mating_mutant_byr2_rep2.pdf")


#Eliminating effect of fluorophores
mating.median.values.batch3 <- 
  mating_unedit_batch3 %>% 
  filter(type=="Wild type") %>% 
  group_by(genotype) %>%
  dplyr::summarize( 
    RV.median.wild = median( relative_v ))
#add the comparison
mating.median.values.batch3 <- mating.median.values.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


mating.median.values.batch3<-select(mating.median.values.batch3, -genotype)
mating_unedit_batch3_normalised<- left_join(mating_unedit_batch3,mating.median.values.batch3,
                                            by="comparison")
mating_unedit_batch3_normalised$viability <- mating_unedit_batch3_normalised$relative_v/mating_unedit_batch3_normalised$RV.median.wild
mating_unedit_batch3_normalised$viability_log<-log(mating_unedit_batch3_normalised$viability)
write_csv(mating_unedit_batch3_normalised, file = "./Databatch3/mating_batch3_normalised.csv")
###boxplot#
mating_unedit_batch3_normalised%>%
  #filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y =  viability_log) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ timepoint ) +  
  ggpubr::stat_compare_means(method = "anova",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Mating viability (batch3)")
ggsave("batch3_mating_byr2_rep2_adjusted.pdf")
mating_unedit_batch3_adjusted_pairwise<- mating_unedit_batch3_normalised%>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggpubr::compare_means(viability_log~genotype,  ., method = "t.test", paired = FALSE,
                        group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(mating_unedit_batch3_adjusted_pairwise, file = "./Databatch3/mating_pairwise-t-test_batch3.csv")



#full cycle
#dfbatch3j<- read.csv(file = 'D:/Wolf project/Data processing/Databatch3/rep2_byr2_main_batch3.csv')
#dfbatch3js<- read.csv(file = 'D:/Wolf project/Data processing/Databatch3/rep2_byr2_analysis_batch3.csv')

#####Batch 2 Full cycle seperating the data points#######
dfbatch3jfull<- dfbatch3j %>% 
  mutate( s = log(wt.after/blue.after) - log(wt.before/blue.before)) %>%
  filter( competition == "Test" )%>%
  filter( timepoint %in% c("FCT","FCB") )%>%
  filter( all.before >500 ) %>%
  filter( all.after >200 )
write.csv(dfbatch3jfull,file = "./Databatch3/Full_cycle_batch3.csv")
dfbatch3jfull%>% 
  filter(competition=="Test")%>% 
  filter(timepoint=="FCB") %>% 
  group_by(genotype) %>% 
  count()
dfbatch3jfull%>%
  filter( competition == 'Test') %>%
  filter( timepoint =="FCB" )%>%
  ggplot(aes( x = type, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid(  timepoint~comparison , scales = "free_y")+ 
  stat_compare_means(method = "t.test",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Full cycle")
ggsave("batch3_competitions_full_cycle_bottom_byr2_rep2.pdf")
dfbatch3jfull%>%
  filter( competition == 'Test') %>%
  filter( timepoint =="FCT" )%>%
  ggplot(aes( x = type, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid(  timepoint~comparison , scales = "free_y")+ 
  stat_compare_means(method = "t.test",label.x=0.5, label.y=0,label="p.signif")+
  ggtitle("Full cycle")
ggsave("batch3_competitions_full_cycle_top_byr2_rep2.pdf")
##significance full cycle####
Fullcycle_top_pairwise_batch3<-dfbatch3jfull %>%
  filter(timepoint=="FCT")%>%
  ggpubr::compare_means(s~genotype,  ., method = "t.test", paired = FALSE,
   group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(Fullcycle_top_pairwise_batch3, file = "./Databatch3/Fullcycle_top_pairwise_batch3.csv")

Fullcycle_bottom_pairwise_batch3<-dfbatch3jfull %>%
  filter(timepoint=="FCB")%>%
  ggpubr::compare_means(s~genotype,  ., method = "t.test", paired = FALSE,
                        group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(Fullcycle_bottom_pairwise_batch3, file = "./Databatch3/Fullcycle_bottom_pairwise_batch3.csv")
####eliminate fluorophore#####
Fullcycle_top_batch3<-dfbatch3jfull %>%
  filter(timepoint=="FCT")
Fullcycle_top_batch3$s_adj=  Fullcycle_top_batch3$s -  min(Fullcycle_top_batch3$s)+1
Fullcycle_bottom_batch3<-dfbatch3jfull %>%
  filter(timepoint=="FCB")
Fullcycle_bottom_batch3$s_adj=  Fullcycle_bottom_batch3$s -  min(Fullcycle_bottom_batch3$s)+1
####################
##normalising selection
#full cycle
#bottom
#get median values of wild type

Full.cycle.median.values.bottom.batch3 <- 
  Fullcycle_bottom_batch3%>% 
  filter(type=="Wild type") %>% 
  group_by(genotype) %>%
  dplyr::summarize( 
    s_B_median_wild = median( s_adj))
#add the comparison
Full.cycle.median.values.bottom.batch3 <- Full.cycle.median.values.bottom.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


Full.cycle.median.values.bottom.batch3<-select(Full.cycle.median.values.bottom.batch3, -genotype)
Fullcycle_bottom_batch3_normalised<- left_join( Fullcycle_bottom_batch3,Full.cycle.median.values.bottom.batch3,by="comparison")
Fullcycle_bottom_batch3_normalised$S_B <- Fullcycle_bottom_batch3_normalised$s_adj/ Fullcycle_bottom_batch3_normalised$s_B_median_wild
Fullcycle_bottom_batch3_normalised$S_B_log<-log(Fullcycle_bottom_batch3_normalised$S_B)
write_csv(Fullcycle_bottom_batch3_normalised, file = "./Databatch3/Fullcycle_bottom_batch3_normalised.csv")
###box plot####
Fullcycle_bottom_batch3_normalised%>%
  #filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y =  S_B_log ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ timepoint ) +  
  ggpubr::stat_compare_means(method = "anova",label.x=0.5, label.y=0.75,label="p.signif")+
  ggtitle("Full cycle Bottom (batch3)")
ggsave("batch3_competitions_full_cycle_bottom_byr2_rep2.pdf")



Fullcycle_bottom_batch3_normalised_pairwise<-Fullcycle_bottom_batch3_normalised%>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggpubr::compare_means(S_B_log~genotype,  ., method = "t.test", paired = FALSE,
                        group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(Fullcycle_bottom_batch3_normalised_pairwise, file = "./Databatch3/Fullcycle_bottom_batch3_normalised_pairwise-t-test_batch3.csv")
####top

Full.cycle.median.values.top.batch3 <- 
  Fullcycle_top_batch3%>% 
  filter(type=="Wild type") %>% 
  group_by(genotype) %>%
  dplyr::summarize( 
    s_T_median_wild = median( s_adj))
#add the comparison
Full.cycle.median.values.top.batch3 <- Full.cycle.median.values.top.batch3 %>% mutate(comparison = case_when(
  endsWith(genotype, "Bf") ~ "rep2",
  endsWith(genotype, "RfBw") ~ "byr2",
  endsWith(genotype, "RfBm") ~ "byr2",
  endsWith(genotype, "RmBm") ~ "Double Mutant",
  endsWith(genotype, "A") ~ "Double Mutant"
))


Full.cycle.median.values.top.batch3<-select(Full.cycle.median.values.top.batch3, -genotype)
Fullcycle_top_batch3_normalised<- left_join( Fullcycle_top_batch3,Full.cycle.median.values.top.batch3,by="comparison")
Fullcycle_top_batch3_normalised$S_T <- Fullcycle_top_batch3_normalised$s_adj/ Fullcycle_top_batch3_normalised$s_T_median_wild
Fullcycle_top_batch3_normalised$S_T_log<-log(Fullcycle_top_batch3_normalised$S_T)
write_csv(Fullcycle_top_batch3_normalised, file = "./Databatch3/Fullcycle_top_batch3_normalised.csv")
###box plot####
Fullcycle_top_batch3_normalised%>%
  #filter( competition == 'Test') %>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggplot( aes(x = genotype, y =  S_T_log ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  #facet_grid( type ~ timepoint ) +  
  ggpubr::stat_compare_means(method = "anova",label.x=0.5, label.y=0.75,label="p.signif")+
  ggtitle("Full cycle top (batch3)")
ggsave("batch3_competitions_full_cycle_top_byr2_rep2_normalised.pdf")



Fullcycle_top_batch3_normalised_pairwise<-Fullcycle_top_batch3_normalised%>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  ggpubr::compare_means(S_T~genotype,  ., method = "t.test", paired = FALSE,
                        group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")
write_csv(Fullcycle_top_batch3_normalised_pairwise, file = "./Databatch3/ecological_sel_top_adjusted_pairwise-t-test_batch3.csv")

B<-Fullcycle_top_batch3_normalised%>%
  filter(genotype %in% c("A","RfBm","RmBf","RmBm")  )%>%
  aov(S_T~genotype, .)
Tukey
#dfbatch3jfull %>%
#filter( competition == 'Test') %>%
#filter( type == "Mutant") %>%
#filter(comparison==c("rep2","byr2")) %>%
#ggplot( aes(x = genotype, y = s ) ) +
#  geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
#facet_grid(  ~ timepoint ) +  
#stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")

#####Significance#######
t.test.growth_batch3 <-compare_means(m.evol_log~genotype,  growth_batch3, method = "t.test", paired = FALSE,
                                     group.by ="type" , ref.group = NULL, p.adjust.method = "bonferroni")

growth_batch3_normalised %>%
  filter(type=="Mutant")%>%
  compare_means(m.evol_divided_log~genotype, ., method = "t.test", paired = FALSE,
                group.by =NULL , ref.group = NULL, p.adjust.method = "bonferroni")

ecological_sel_bottom_batch3 <-ecological_sel_batch3 %>%
  filter( timepoint == 'BAG') 
wilcox.eco.select_bottom_batch3<-compare_means(s_evol_difflog~genotype, ecological_sel_bottom_batch3, method = "wilcox.test", paired = FALSE,
                                               group.by = "comparison"  , ref.group = NULL)

ecological_sel_top_batch3<-ecological_sel_batch3 %>%
  filter( timepoint == 'TAG') 
wilcox.eco.select_top_batch3<-
  compare_means(WB_divided_log~genotype, ecological_sel_bottom_batch3_normalised, method = "anova", 
                paired = FALSE,group.by = NULL, ref.group = NULL,p.adjust.method = "bonferroni")
summary (aov(WB_divided_log~genotype,ecological_sel_bottom_batch3_normalised))

pairwise.t.test(ecological_sel_bottom_batch3_normalised$WB_divided_log,
                ecological_sel_bottom_batch3_normalised$genotype,p.adjust.method="bonferroni")




wilcox.mating_batch3<-compare_means(relative_v~genotype,mating_unedit_batch3, method = "wilcox.test", paired = FALSE,
                                    group.by = "comparison"  , ref.group = NULL)
full_cycle_bottom_batch3 <-dfbatch3jfull %>%
  filter( timepoint == 'FCB') 
write_csv(full_cycle_bottom_batch3, file = "./Databatch3/full_cycle_bottom_batch3.csv")

wilcox.full_cycle_bottom_batch3<-compare_means(s~genotype,full_cycle_bottom_batch3, method = "t.test", paired = FALSE,
                                               group.by = "comparison" , ref.group = NULL,p.adjust.method = "bonferroni")
full_cycle_top_batch3<-dfbatch3jfull %>%
  filter( timepoint == 'FCT') 
write_csv(full_cycle_top_batch3, file = "./Databatch3/full_cycle_top_batch3.csv")
wilcox.full_cycle_top_batch3<-compare_means(s~genotype,full_cycle_top_batch3, method = "wilcox.test", paired = FALSE,
                                            group.by = "comparison"  , ref.group = NULL)
write_csv(wilcox.mating_batch3, file = "./Databatch3/wilcox_mating_batch3.csv")
write_csv(wilcox.growth_batch3, file = "./Databatch1/wilcox_growth_batch3.csv")
write_csv(wilcox.eco.select_bottom_batch3, file = "./Databatch3/wilcox_eco.select_bottom_batch3.csv")
write_csv(wilcox.eco.select_top_batch3, file = "./Databatch3/wilcox_eco.select_top_batch3.csv")
write_csv(wilcox.full_cycle_bottom_batch3, file = "./Databatch3/wilcox_full_cycle_bottom_batch3.csv")
write_csv(wilcox.full_cycle_top_batch3, file = "./Databatch3/wilcox_full_cycle_top_batch3.csv")

#summarise(aov.growth_batch3)
#aov.growth_batch3 <-kruskal.test(m.evol~genotype,growth,comparison$growth)
#wilcox.test(m.evol~genotype,growth)

#summarise(aov.growth_batch3)
#shapiro.test(m.evol$growth)
#attach(growth)
#CConsolidate data ?
#dfbatch3_consolidate <- dfbatch3
growth_batch3%>%
  filter( competition == "Test" )%>%
  filter( genotype %in% c("RmBf","RfBm") )%>%
  #filter( type == "Mutant" )%>%
  ggplot(aes( x = genotype, y = m.evol ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
ggsave("batch3_growth_mutant_byr2_rep2.pdf")
ecological_sel_batch3%>%
  filter( competition == "Test" )%>%
  #filter( genotype == c("RmBf","RfBm") )%>%
  filter( type == "Mutant" )%>%
  ggplot(aes( x = genotype, y = s_evol_difflog ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
mating_unedit_batch3%>%
  filter( competition == "Test" )%>%
  #filter( genotype == c("RmBf","RfBm") )%>%
  filter( type == "Mutant" )%>%
  ggplot(aes( x = genotype, y = relative_v ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
dfbatch3jfull%>%
  filter( competition == "Test" )%>%
  #filter( genotype == c("RmBf","RfBm") )%>%
  filter( type == "Mutant" )%>%
  ggplot(aes( x = genotype, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=0.5,label="p.signif")

growth_batch3%>%
  filter( competition == "Test" )%>%
  filter( type == "Wild type" )%>%
  ggplot(aes( x = genotype, y = m.evol ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
ecological_sel_batch3%>%
  filter( competition == "Test" )%>%
  filter( type == "Wild type" )%>%
  ggplot(aes( x = genotype, y = s_evol_difflog ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")

mating_unedit_batch3%>%
  filter( competition == "Test" )%>%
  filter( type == "Wild type" )%>%
  ggplot(aes( x = genotype, y = relative_v ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=2,label="p.signif")
dfbatch3jfull%>%
  filter( competition == "Test" )%>%
  filter( timepoint == "FCB" )%>%
  #filter( type == "Mutant" )%>%
  ggplot(aes( x = type, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid(  timepoint~ comparison , scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=0.5,label="p.signif")+
  ggtitle("Full cycle bottom")
dfbatch3jfull%>%
  filter( competition == "Test" )%>%
  filter( timepoint == "FCT" )%>%
  #filter( type == "Mutant" )%>%
  ggplot(aes( x = comparison, y = s ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid(  timepoint~ type , scales = "free_y")+ 
  stat_compare_means(method = "wilcox.test",label.x=0.5, label.y=0.5,label="p.signif")+
  ggtitle("Full cycle top")
mating_unedit_batch3%>%
  filter( competition == "Test" )%>%
  #filter( genotype %in% c("RmBf","RfBm") )%>%
  filter( type == "Mutant" )%>%
  ggplot(aes( x = genotype, y = relative_v ) ) +
    geom_boxplot() +    stat_summary(fun.data = give.n, geom = "text") +
  facet_grid( timepoint ~ type, scales = "free_y")+ 
  stat_compare_means(method = "anova",label.x=0.5, label.y=2,label="p.signif")\