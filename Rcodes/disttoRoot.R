#!/usr/bin/Rscript
# Marco G
# Proyecto: aeneus
# Eval. de posibles candidatos de LBA (Long Branch Attraction)
# 30/sep/2026

# ============================================================
# 0. LIBRERÍAS
# ============================================================
library(ape)          # árboles filogenéticos
library(adephylo)     # distancias raíz-punta
library(pheatmap)     # heatmap de distancias
library(ggplot2)      # gráficas
library(ggrepel)      # etiquetas sin solapamiento
library(dplyr)        # manipulación de datos

# ============================================================
# 1. DIRECTORIO DE TRABAJO
# ============================================================
setwd("/run/media/pinkpunk/Myriad/aeneus/results/phylo")

# ============================================================
# 2. ÁRBOL ASC — DISTANCIAS, HEATMAP Y BOXPLOT
# ============================================================
aeneusASC <- read.tree("asc/amex-RAD-asc.treefile")

# --- 2.1 Distancias raíz-punta ---
dist_root_tip <- distRoot(aeneusASC, method = "patristic")
hist(dist_root_tip, breaks = 20,
     main = "Distancias de raíz a punta (ASC)",
     xlab = "Distancia")

umbral   <- mean(dist_root_tip) + 2 * sd(dist_root_tip)
outliers <- names(dist_root_tip[dist_root_tip > umbral])
outliers

# --- 2.2 Matriz de distancias patrísticas y anotación ---
D <- as.matrix(distTips(aeneusASC, method = "patristic"))

ann <- data.frame(
  root_to_tip = dist_root_tip,
  outlier     = ifelse(dist_root_tip > umbral, "Outlier", "Normal")
)
rownames(ann) <- names(dist_root_tip)
ann$samples   <- names(dist_root_tip)

# --- 2.3 Heatmap ---
pheatmap(
  mat = D,
  clustering_distance_rows = as.dist(D),
  clustering_distance_cols = as.dist(D),
  annotation_row           = ann[, c("root_to_tip", "outlier")],
  annotation_col           = ann[, c("root_to_tip", "outlier")],
  show_rownames            = FALSE,
  show_colnames            = FALSE,
  fontsize_row             = 4,
  main                     = "Distancias patrísticas + distancia raíz-punta (ASC)"
)

# --- 2.4 Boxplot de distancias raíz-punta (outliers etiquetados) ---
outliers_df <- subset(ann, root_to_tip > umbral)

ggplot(ann, aes(x = factor(0), y = root_to_tip)) +
  geom_boxplot(width = 0.25, fill = NA, outlier.colour = NA) +
  geom_jitter(width = 0.05, color = "grey60", size = 1, alpha = 0.6) +
  geom_point(data = outliers_df, aes(x = factor(0)),
             color = "red", size = 2.5) +
  geom_label_repel(
    data = outliers_df, aes(label = samples),
    size = 3, color = "red",
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3,
    min.segment.length = 0, direction = "y", nudge_x = 0.6,
    segment.color = "red", segment.size = 0.3, force = 3
  ) +
  labs(x = NULL, y = "Distancia raíz-punta (patristic)",
       title = "ASC") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ============================================================
# 3. ÁRBOL NOASC — DISTANCIAS, HEATMAP Y BOXPLOT
# ============================================================
aeneusNOASC <- read.tree("noasc/amex-RAD-NOASC.treefile")

# --- 3.1 Distancias raíz-punta ---
dist_root_tip <- distRoot(aeneusNOASC, method = "patristic")
hist(dist_root_tip, breaks = 20,
     main = "Distancias de raíz a punta (NOASC)",
     xlab = "Distancia")

umbral   <- mean(dist_root_tip) + 2 * sd(dist_root_tip)
outliers <- names(dist_root_tip[dist_root_tip > umbral])
outliers

# --- 3.2 Matriz y anotación ---
D <- as.matrix(distTips(aeneusNOASC, method = "patristic"))

ann <- data.frame(
  root_to_tip = dist_root_tip,
  outlier     = ifelse(dist_root_tip > umbral, "Outlier", "Normal")
)
rownames(ann) <- names(dist_root_tip)
ann$samples   <- names(dist_root_tip)

# --- 3.3 Heatmap ---
pheatmap(
  mat = D,
  clustering_distance_rows = as.dist(D),
  clustering_distance_cols = as.dist(D),
  annotation_row           = ann[, c("root_to_tip", "outlier")],
  annotation_col           = ann[, c("root_to_tip", "outlier")],
  show_rownames            = FALSE,
  show_colnames            = FALSE,
  fontsize_row             = 4,
  main                     = "Distancias patrísticas + distancia raíz-punta (NOASC)"
)

# --- 3.4 Boxplot ---
outliers_df <- subset(ann, root_to_tip > umbral)

ggplot(ann, aes(x = factor(0), y = root_to_tip)) +
  geom_boxplot(width = 0.25, fill = NA, outlier.colour = NA) +
  geom_point(data = outliers_df, aes(x = factor(0)),
             color = "red", size = 2.5) +
  geom_label_repel(
    data = outliers_df, aes(label = samples),
    size = 3, color = "red",
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3,
    min.segment.length = 0, direction = "y", nudge_x = 0.6,
    segment.color = "red", segment.size = 0.3, force = 3
  ) +
  labs(x = NULL, y = "Distancia raíz-punta (patristic)",
       title = "NOASC") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ============================================================
# 4. CLADOS Y FIRMAS DE TREESHRINK
# ============================================================
# --- 4.1 Asignación de clados ---
clados <- read.delim(file = "./cladosInfo.tsv", header = FALSE, sep = "\t")
colnames(clados) <- c("sample", "clado")

# --- 4.2 Firmas de TreeShrink (solo taxa que consideró problemáticos) ---
# El archivo output_summary.txt tiene columnas: Gene, Species, Taxon, Signature
firmas <- read.table("treeshrink_noasc/output_summary.txt",
                     header = TRUE, stringsAsFactors = FALSE)
firmas <- firmas[, c("Taxon", "Signature")]
colnames(firmas) <- c("sample", "signature")

# ============================================================
# 5. DATA FRAME MAESTRO (distancias + clados + firmas)
# ============================================================
df <- data.frame(
  sample = names(dist_root_tip),
  dist   = as.numeric(dist_root_tip)
) %>%
  left_join(clados, by = "sample") %>%
  left_join(firmas, by = "sample") %>%
  mutate(
    clado     = ifelse(is.na(clado), "Sin asignar", clado),
    outlier   = ifelse(dist > umbral, "Outlier", "Normal"),
    has_firma = !is.na(signature)
  ) %>%
  arrange(desc(dist)) %>%
  mutate(sample = factor(sample, levels = sample))

# Verificación
table(df$clado, useNA = "ifany")
cat("Tips totales:", nrow(df), "| Con firma:", sum(df$has_firma), "\n")

# ============================================================
# 6. PALETA DE COLORES
# ============================================================
paleta_final <- c(
  "#0077BB",  # Azul        -> Clado 1
  "#33BBEE",  # Celeste     -> Clado 2
  "#009988",  # Teal        -> Clado 3
  "#44BB99",  # Verde azul  -> Clado 4
  "#EE7733",  # Naranja     -> Clado 5
  "#CC3311",  # Rojo anaranjado -> outgroup
  "#EE3377",  # Magenta
  "#BBBBBB",  # Gris
  "#000000",  # Negro
  "#DDDDDD",  # Gris claro
  "#882255"   # Vino
)

# ============================================================
# 7. FIGURA: BARRAS POR CLADO + FIRMAS
# ============================================================
# Etiquetas de firmas: se dibujan encima de la barra, rotadas 90°.
# Solo aparecen para taxa que TreeShrink incluyó en output_summary.txt.
offset <- max(df$dist) * 0.015

p_barras <- ggplot(df, aes(x = sample, y = dist, fill = clado)) +
  geom_col(width = 1) +
  geom_hline(yintercept = umbral, linetype = "dashed",
             linewidth = 0.7, color = "black") +
  # Etiquetas de firmas (solo taxa con firma)
  geom_text(
    data = subset(df, has_firma),
    aes(label = round(signature, 4), y = dist + offset),
    angle = 90, hjust = 0, vjust = 0.5,
    size = 1.8, color = "black"
  ) +
  scale_fill_manual(values = paleta_final) +
  labs(
    x = NULL,
    y = "Distancia raíz-punta (patristic)",
    title = "Distancias raíz-punta por clado (firmas TreeShrink)",
    fill  = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 90, size = 6),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_barras)

# ============================================================
# 8. FIGURA: BOXPLOT POR CLADO
# ============================================================
p_box <- ggplot(df, aes(x = reorder(clado, dist, FUN = median),
                        y = dist, fill = clado)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  geom_hline(yintercept = umbral, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = paleta_final) +
  labs(x = NULL, y = "Distancia raíz-punta (patristic)",
       title = "Distribución de distancias por clado") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_box)

  # ============================================================
# 9. TABLA RESUMEN POR CLADO
# ============================================================
resumen <- df %>%
  group_by(clado) %>%
  summarise(
    n          = n(),
    mediana    = round(median(dist), 4),
    media      = round(mean(dist), 4),
    sd         = round(sd(dist), 4),
    n_outliers = sum(outlier == "Outlier"),
    prop_out   = paste0(round(100 * mean(outlier == "Outlier"), 1), "%")
  ) %>%
  arrange(desc(mediana))

print(resumen)

# ============================================================
# 10. GUARDAR FIGURAS
# ============================================================
ggsave("fig_barras_clados_firmas.png", p_barras,
       width = 18, height = 7, dpi = 300)
ggsave("fig_boxplot_clados.png", p_box,
       width = 8, height = 6, dpi = 300)