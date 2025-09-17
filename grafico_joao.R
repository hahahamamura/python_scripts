library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

# Caminho do arquivo
dados <- "/home/lab/Downloads/piloto_hla/summaries/all_windows_summary.tsv"

# Leitura (TSV, não CSV)
dados <- read_tsv(dados)

# Gráfico 1: número de alelos ao longo da posição inicial
dados %>%
  mutate(
    allele_group = case_when(
      num_alleles <= 5  ~ "3-5",
      num_alleles <= 10 ~ "6-10",
      num_alleles <= 20 ~ "11-20",
      TRUE              ~ "21+"
    ),
    allele_group = factor(allele_group, levels = c("21+", "11-20", "6-10", "3-5")) # ordem da maior para menor # nolint
  ) %>%
  ggplot(aes(x = start, y = num_alleles, color = allele_group)) +
  geom_point(alpha = 0.7) +            
  geom_line(aes(group = 1), color = "grey70", alpha = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "red", linetype = "dashed") + # nolint
  labs(x = "Posição inicial", y = "Número de alelos",
       title = "Número de alelos por microhaplótipo",
       color = "Número de alelos") +
  theme_minimal() +
  theme(legend.position = "right")

#salva
ggsave(filename = "plot_num_alleles.png", bg = "white", width = 20, limitsize = FALSE) # nolint


# Gráfico 2: Ae
dados %>%
  mutate(
    Ae_group = case_when(
      Ae <= 2  ~ "0-2",
      Ae <= 4  ~ "2-4",
      Ae <= 6  ~ "4-6",
      TRUE     ~ "6+"
    ),
    Ae_group = factor(Ae_group, levels = c("6+", "4-6", "2-4", "0-2")) # ordem da maior para menor
  ) %>%
  ggplot(aes(x = start, y = Ae, color = Ae_group)) +
  geom_point(alpha = 0.7) +            
  geom_line(aes(group = 1), color = "grey70", alpha = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "red", linetype = "dashed") +
  labs(x = "Posição inicial", y = "Ae",
       title = "Ae por microhaplótipo",
       color = "Ae (intervalos)") +
  theme_minimal() +
  theme(legend.position = "right")
#salva
ggsave(filename = "plot_ae.png", bg = "white", width = 20, limitsize = FALSE) # nolint



# Gráfico 3: Fis e PIC
dados %>%
  select(start, Fis, PIC) %>%
  pivot_longer(cols = -start, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = start, y = value, color = metric)) +
  geom_line(size = 1) +      # linhas
  labs(x = "Posição inicial", y = "Valor",
       title = "Ae, Fis e PIC por microhaplótipo",
       color = "Métrica") +  # legenda única
  theme_minimal() +
  theme(legend.position = "right")
#salva
ggsave(filename = "plot_fis_pic.png", bg = "white", width = 20, limitsize = FALSE) # nolint



# Gráfico He, Ho e HW 
dados %>%
  select(start, He, Ho) %>%
  pivot_longer(cols = -start, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = start, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_point(alpha = 0.6) +
  labs(x = "Posição inicial", y = "Valor",
       title = "He, Ho e HWE por microhaplótipo",
       color = "Métrica") +
  theme_minimal()
#salva
ggsave(filename = "plot_he_ho.png", bg = "white", width = 20, limitsize = FALSE) # nolint

# Gráfico MP e PD 
dados %>%
  select(start, MP, PD) %>%
  pivot_longer(cols = -start, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = start, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_point(alpha = 0.6) +
  labs(x = "Posição inicial", y = "%",
       title = "MP e PD por microhaplótipo",
       color = "Métrica") +
  theme_minimal()
#salva
ggsave(filename = "plot_mp_pd.png", bg = "white", width = 20, limitsize = FALSE) # nolint




# Gráfico PROB com linhas coloridas e pontos tipo heatmap
dados %>%
  select(start, prob_2_diff, prob_3_diff, prob_4_diff) %>%
  pivot_longer(cols = -start, names_to = "metric", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%  # garante que value é numérico
  ggplot(aes(x = start, y = value)) +
  geom_line(aes(color = metric), size = 1) +              # linhas coloridas por métrica
  geom_point(aes(fill = value), shape = 21, color = "black", size = 2) +  # pontos heatmap
  scale_fill_viridis_c(option = "C") +                   # escala contínua para os pontos
  labs(x = "Posição inicial", y = "%",
       title = "Probabilidades de alelos diferentes (2, 3 e 4) por microhaplótipo",
       color = "Métrica", fill = "Valor (heatmap)") +
  theme_minimal()
#salva
ggsave(filename = "plot_probs.png", bg = "white", width = 20, limitsize = FALSE) # nolint
