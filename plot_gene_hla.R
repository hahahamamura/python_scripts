library(ggplot2)
library(cowplot)
library(tidyverse)

# ----------------------------
# Dados
# ----------------------------
dados <- read_tsv("/home/lab/Downloads/piloto_hla/summaries/all_windows_summary.tsv")

gene_start <- 29942254
gene_end   <- 29945755

exons <- tibble(
  start = c(29942554, 29942757, 29943268, 29944122,
            29944500, 29945059, 29945234, 29945451),
  end   = c(29942626, 29943026, 29943543, 29944397,
            29944616, 29945091, 29945281, 29945455),
  exon  = paste0("Exon", 1:8)
)

utr <- tibble(
  start = c(exons$start[1] - 22, exons$end[nrow(exons)]),
  end   = c(exons$start[1], exons$end[nrow(exons)] + 300)
)

# ----------------------------
# Gráfico 1 - Número de alelos
# ----------------------------
p1 <- dados %>%
  mutate(
    allele_group = case_when(
      num_alleles <= 5  ~ "3-5",
      num_alleles <= 10 ~ "6-10",
      num_alleles <= 20 ~ "11-20",
      TRUE              ~ "21+"
    ),
    allele_group = factor(allele_group, levels = c("21+", "11-20", "6-10", "3-5"))
  ) %>%
  ggplot(aes(x = start, y = num_alleles, color = allele_group)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = 1), color = "grey70", alpha = 0.5) +
  labs(
    x = NULL, y = "Número de alelos",
    title = "Número de alelos por microhaplótipo",
    color = "Número de alelos"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# ----------------------------
# Gráfico 2 - Estrutura do gene HLA-A
# ----------------------------
p2 <- ggplot() +
  geom_segment(aes(x = gene_start, xend = gene_end, y = 0, yend = 0), size = 0.4) +
  geom_rect(data = exons,
            aes(xmin = start, xmax = end, ymin = -0.08, ymax = 0.08),
            fill = "#0D00C0") +
  geom_rect(data = utr,
            aes(xmin = start, xmax = end, ymin = -0.04, ymax = 0.04),
            fill = "#0D00C0") +
  scale_x_continuous(
    breaks = (exons$start + exons$end) / 2,
    labels = exons$exon,
    limits = c(gene_start, gene_end),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  labs(x = "Estrutura do gene HLA-A", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, 5, 0, 5)
  )

# ----------------------------
# Combinar gráficos
# ----------------------------
final_plot <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(3, 0.6))
print(final_plot)

ggsave("plot_num_alleles_gene_hla_a.png",
       final_plot, bg = "white", width = 20, limitsize = FALSE)

