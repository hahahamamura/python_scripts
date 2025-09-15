library(tidyverse)

dados <- "/home/lab/Downloads/piloto_hla/analise_estatistica_microhaplotipos.csv"


dados <- read_csv(dados)  


dados %>% select(Posicao,N_Alelos) %>% mutate(Start = str_split_i(Posicao,"-",1)) %>%
    ggplot(aes(x=as.numeric(Start), y = N_Alelos, group = 1)) +
    geom_point() +
    geom_line()


dados %>% select(Heterozigosidade_Esperada:Poder_Discriminacao) %>%
   pivot_longer(cols = everything()) %>%
   ggplot(aes(x = value, fill = name)) +
   geom_histogram() +
   facet_wrap(~name,scale = "free_x")

dados %>% select(Posicao,Heterozigosidade_Esperada:Poder_Discriminacao) %>% 
    mutate(Start = str_split_i(Posicao,"-",1),
           Start = as.numeric(Start)) %>%
    pivot_longer(cols = Heterozigosidade_Esperada:Poder_Discriminacao) %>%
    ggplot(aes(x=Start, y = value, color = name)) +
    geom_point() +
    geom_line() +
    facet_wrap(~name, scale = "free_y",nrow = 5)
