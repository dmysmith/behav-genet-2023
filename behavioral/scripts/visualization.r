# bar graphs of FEMA random effects (sig2mat)
# Diana Smith
# May 2022

library(ggplot2)
library(viridis)
library(hrbrthemes)

library(dplyr)
library(tidyr)
library(tidyverse)
library(weights)

# define paths
inpath = '/Users/dsmith/OneDrive\ -\ UC\ San\ Diego/PhD/ABCD\ Projects/random_effects/results/behavioral'
outpath = '/Users/dsmith/OneDrive\ -\ UC\ San\ Diego/PhD/ABCD\ Projects/random_effects/plots/behavioral'

designmat_list = c('designMat1_allcovs',
                   'designMat2_agesex',
                   'designMat3_agesexPCs',
                   'designMat4_agesexhispSES')

random_effects = c('FAE',
                   'FSE',
                   'FASE',
                   'FADSE')

phenotypes = c('nihtbx picvocab','nihtbx flanker','nihtbx pattern','nihtbx picture','nihtbx reading','nihtbx cryst','lmt scr perc correct','anthroheightcalc')

# generate key-value pairs for color coding
keys = c('F','A','D','S','E')
values = c('#003f5c',  '#58508d',  '#bc5090',  '#ff6361',  '#ffa600')
colors = list()
for (i in 1:length(keys)){colors[[keys[i]]] <- values[i]}

for (i in 1:length(designmat_list)){
  for (j in 1:length(random_effects)){
    sig2mat_file = paste('sig2mat_',random_effects[j], '.csv',sep='')
    infile = paste(inpath,designmat_list[i],sig2mat_file,sep='/')
    df <- read.csv(infile, header = F)
    colnames(df) = phenotypes
    rownames(df) = as.list(strsplit(random_effects[j], "")[[1]])
    
    p <- df %>%
      rownames_to_column(var = 'x') %>%
      pivot_longer(-x) %>%
      filter(value > 0) %>% 
      mutate(x = factor(x, levels = c('F', 'A', 'D', 'S', 'E'))) %>% 
      ggplot(aes(fill = x, y = value, x = forcats::fct_rev(name), label = x)) +
      geom_bar(position="stack", stat="identity") +
      geom_text(aes(label=rd(value,digits=2)), position = position_stack(vjust = 0.5), color="black", size=3.5) +
      theme_ipsum() +
      xlab("") +
      ylab("") +
      theme(legend.title = element_blank()) +
      ggtitle(paste('Model including: ',substr(designmat_list[i], 12,nchar(designmat_list[i])), ', ', random_effects[j], sep = ''))+
      scale_fill_manual(values = colors[strsplit(random_effects[j], "")[[1]]])
    p
    
    plotname = paste('sig2mat_',substr(designmat_list[i], 12,nchar(designmat_list[i])),'_',random_effects[j], '.jpg',sep='')
    outfile = paste(outpath, plotname, sep='/')
    ggsave(p, filename = outfile, width = 12)
  }
}

