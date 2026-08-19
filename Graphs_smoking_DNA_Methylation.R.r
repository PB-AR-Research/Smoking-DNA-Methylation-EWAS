library(readr)
library(ggplot2)
library(ggrepel)

data <- read_csv("volcano_data.csv")
names(data)
data$Category <- "Not Significant"

data$Category[data$`Mean delta beta` >= 0.02 & data$FDR < 0.05] <- "Hypermethylated"

data$Category[data$`Mean delta beta` <= -0.02 & data$FDR < 0.05] <- "Hypomethylated"
#plots
ggplot(data, aes(x = `Mean delta beta`,
                 y = -log10(FDR),
                 color = Category)) +

  geom_point(size = 2, alpha = 0.8) +

  scale_color_manual(values = c(
    "Hypermethylated" = "red",
    "Hypomethylated" = "blue",
    "Not Significant" = "gray"
  )) +

  geom_vline(xintercept = c(-0.02, 0.02), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +

  theme_minimal(base_size = 14)
  library(ggplot2)
library(ggrepel)
top_hits <- subset(data, FDR < 0.001)

ggplot(data, aes(x = `Mean delta beta`,
                 y = -log10(FDR),
                 color = Category)) +

  geom_point(size = 2, alpha = 0.8) +

  geom_text_repel(data = top_hits,
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = 15) +

  scale_color_manual(values = c(
    "Hypermethylated" = "#D73027",
    "Hypomethylated" = "#4575B4",
    "Not Significant" = "grey70"
  )) +

  geom_vline(xintercept = c(-0.02, 0.02),
             linetype = "dashed",
             linewidth = 0.5) +

  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.5) +

  theme_minimal(base_size = 14) +

  labs(
    title = "Volcano Plot of Smoking-Associated DNA Methylation",
    x = expression(Delta*beta),
    y = expression(-log[10](FDR)),
    color = "Methylation Status"
  ) +

  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
  ggsave("volcano_thesis.pdf", width = 8, height = 6)
  ggsave("volcano_thesis_300dpi.png", dpi = 300, width = 8, height = 6)
  ggsave("volcano_thesis_600dpi.png", dpi = 600, width = 8, height = 6)
  #manhattan plot
  data$logFDR <- -log10(data$FDR)
data$index <- 1:nrow(data)

ggplot(data, aes(x = index, y = logFDR)) +

  geom_point(aes(color = logFDR), alpha = 0.8, size = 1.8) +

  scale_color_gradient(low = "blue", high = "red") +

  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +

  theme_classic(base_size = 14) +

  labs(
    title = "CpG Manhattan-Style Plot",
    x = "CpG Index",
    y = expression(-log[10](FDR))
  ) +

  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
  p <- ggplot(data, aes(x = index, y = logFDR)) +

  geom_point(aes(color = logFDR), alpha = 0.8, size = 1.8) +

  scale_color_gradient(low = "blue", high = "red") +

  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +

  theme_classic(base_size = 14) +

  labs(
    title = "CpG Manhattan-Style Plot",
    x = "CpG Index",
    y = expression(-log[10](FDR))
  ) +

  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
  ggsave("manhattan_plot_thesis.pdf",
       plot = p,
       width = 10,
       height = 5)
       ggsave("manhattan_plot_300dpi.png",
       plot = p,
       width = 10,
       height = 5,
       dpi = 300)
       ggsave("manhattan_plot_600dpi.png",
       plot = p,
       width = 10,
       height = 5,
       dpi = 600)
       library(ggplot2)

top <- data[order(data$FDR), ][1:30, ]

p <- ggplot(top, aes(x = Gene, y = 1, fill = `Mean delta beta`)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) +
  labs(title = "Top Differentially Methylated CpGs")
  ggsave("heatmap_thesis.pdf",
       plot = p,
       width = 10,
       height = 5)
       ggsave("heatmap_300dpi.png",
       plot = p,
       width = 10,
       height = 5,
       dpi = 300)
       ggsave("heatmap_600dpi.png",
       plot = p,
       width = 10,
       height = 5,
       dpi = 600)
       #effect size distribution plot
library(ggplot2)

p1 <- ggplot(data, aes(x = `Mean delta beta`)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of DNA Methylation Changes",
    x = "Mean Delta Beta",
    y = "Count"
  )
  p1
  ggsave("histogram_thesis.pdf",
       plot = p1,
       width = 8,
       height = 5)
       ggsave("histogram_300dpi.png",
       plot = p1,
       width = 8,
       height = 5,
       dpi = 300)
       #boxplot: hyper vs hypo methylation
  library(ggplot2)

top_gene <- data[order(abs(data$`Mean delta beta`), decreasing = TRUE), ][1:20, ]

p3 <- ggplot(top_gene, aes(x = reorder(Gene, `Mean delta beta`),
                           y = `Mean delta beta`,
                           fill = `Mean delta beta`)) +

  geom_bar(stat = "identity") +
  coord_flip() +

  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +

  theme_classic(base_size = 13) +

  labs(
    title = "Top Differentially Methylated Genes",
    x = "Gene",
    y = "Mean Delta Beta"
  )
  p3
  ggsave("gene_barplot_thesis.pdf",
       plot = p3,
       width = 8,
       height = 6)
       ggsave("gene_barplot_300dpi.png",
       plot = p3,
       width = 8,
       height = 6,
       dpi = 300)
#FDR vs Effect Size Scatter Plot
       library(ggplot2)

p5 <- ggplot(data, aes(x = `Mean delta beta`,
                       y = -log10(FDR),
                       color = Category)) +

  geom_point(alpha = 0.7, size = 2) +

  scale_color_manual(values = c(
    "Hypermethylated" = "red",
    "Hypomethylated" = "blue",
    "Not Significant" = "gray"
  )) +

  theme_classic(base_size = 14) +

  labs(
    title = "Effect Size vs Statistical Significance in CpG Methylation",
    x = "Mean Delta Beta (Effect Size)",
    y = expression(-log[10](FDR)),
    color = "Methylation Status"
  )