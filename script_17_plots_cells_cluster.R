library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

md <- readRDS("output/seurat_metadata_subset.rds")

md %>% filter(celltype %in% c("micro.0", "micro.1", "micro.3", "micro.4", "micro.5", "micro.6",
                              "micro.7", "micro.9", "micro.12", "micro.13", "micro.17", "micro.23",
                              "micro.25", "micro.36")) %>% nrow()

md %>% filter(celltype %in% c("micro.36")) %>% nrow()

#md %>% filter(celltype == "micro.4") %>% rownames() %>% write.table("output/files/micro_4_root_cells.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plt_dt <- md %>% 
  select(orig.ident, Treatment, Region, celltype) %>% 
  count(orig.ident, celltype, Treatment, Region, .drop=FALSE) %>% 
  mutate(percent = 0)
rm(md)

for (i in unique(plt_dt$orig.ident)) {
  print(i)
  total <- sum(plt_dt %>% filter(orig.ident==i) %>% pull(n))
  plt_dt[plt_dt$orig.ident==i,]$percent <- plt_dt[plt_dt$orig.ident==i,]$n/total
}

plt_dt[is.na(plt_dt)]

#Add meta data to missing values
tmp <- readxl::read_excel("library_meta_data.xlsx") %>%
  rename(orig.ident = Library)
plt_dt <- plt_dt %>% left_join(tmp, by = "orig.ident")
plt_dt <- plt_dt %>% select(orig.ident, celltype, Region.y, Treatment.y, n, percent)
plt_dt <- plt_dt %>% rename(Region = Region.y, Treatment = Treatment.y)
rm(tmp)
plt_dt[is.na(plt_dt)]

plt_dt <- plt_dt %>%
  mutate(Treatment.Region = paste(Treatment, Region, sep = "_")) %>%
  mutate(celltype = factor(celltype, levels = c("astro.10", "astro.20", "endo.2", "endo.11",
                                                "peri.14", "smc.21", "endo.27", "vlmc.30",
                                                "endo.38", "opc.19", "opc.37", "pvm.22",
                                                "pvm.29", "epen.26", "lymph.28", "micro.3",
                                                "micro.4", "micro.25", "micro.0", "micro.1",
                                                "micro.5", "micro.6", "micro.7", "micro.9",
                                                "micro.12", "micro.13", "micro.17", "micro.23",
                                                "micro.36", "g.neu.16", "g.neu.31", "g.neu.33",
                                                "g.neu.34", "im.neu.18", "im.neu.24", "im.neu.32")))

plt_dt <- plt_dt %>%
  mutate(Treatment.Region = factor(Treatment.Region, levels = c("Mock_cx", "RML_cx","Mock_hp", "RML_hp")))

aggregate(percent~celltype+Treatment.Region, plt_dt, mean)

ggplot(plt_dt, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(height=0, width=0.1, alpha = 0.5) +
  #geom_point(data = aggregate(percent~celltype+Treatment.Region, plt_dt, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
  scale_color_manual(values = c("cornflowerblue", "firebrick1", "navy", "firebrick4")) +
  facet_wrap(~celltype, scales="free_y") +
  theme_classic()
ggsave("output/plots/population plots/all.pdf", width = 10, height = 7.5, units = "in")

ggplot(plt_dt, aes(x=orig.ident, y= percent, fill=celltype)) +
  geom_bar(stat = "identity", position = "stack")

out <- list()
for (i in unique(plt_dt$celltype)) {
  for (j in unique(plt_dt$Region)) {
    mock <- plt_dt %>% filter(celltype == i, Region == j, Treatment == "Mock") %>% pull(percent)
    rml <- plt_dt %>% filter(celltype == i, Region == j, Treatment == "RML") %>% pull(percent)
    
    tmp <- wilcox.test(rml, mock)
    
    out[[paste0(i, j)]] <- data.frame(celltype = i,
                                     Region = j,
                                     pval = tmp$p.value)
  }
}
out <- do.call(rbind, out)
write.csv(out, "output/files/statistics/full_atlas_population_stats.csv")

out %>% filter(pval < 0.12)

micro_clusters <- c("micro.0", "micro.1", "micro.3", "micro.4", "micro.5", "micro.6",
                    "micro.7", "micro.9", "micro.12", "micro.13", "micro.17", "micro.23",
                    "micro.25", "micro.36")
tmp <- plt_dt %>% filter(celltype %in% micro_clusters)

out <- list()
for (i in unique(tmp$celltype)) {
  print(i)
  out[[i]] <- data.frame(celltype = i,
                         enrichment_hp_cx = (tmp %>% filter(tmp$celltype == i, Region == "hp") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx") %>% pull(percent) %>% mean()),
                         enrichment_hp_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "Mock") %>% pull(percent) %>% mean()),
                         enrichment_cx_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "Mock") %>% pull(percent) %>% mean())
                         )
  }
out <- do.call(rbind, out)

out[is.infinite(out$enrichment_hp_RML_Mock),]$enrichment_hp_RML_Mock <- max(out[is.infinite(out$enrichment_hp_RML_Mock)==FALSE,]$enrichment_hp_RML_Mock)
out$enrichment_hp_cx <- log2(out$enrichment_hp_cx)
out$enrichment_hp_RML_Mock <- log2(out$enrichment_hp_RML_Mock)
out$enrichment_cx_RML_Mock <- log2(out$enrichment_cx_RML_Mock)

write.csv(out, "output/files/statistics/enrichment_microglial_clusters.csv")

tmp %>% filter(tmp$celltype == i)

(tmp %>% filter(tmp$celltype == i, Region == "hp") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx") %>% pull(percent) %>% mean())
(tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "Mock") %>% pull(percent) %>% mean())
(tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "Mock") %>% pull(percent) %>% mean())


# keep <- c("astro.10", "astro.20")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic() +
#   theme(legend.position = "none")
# ggsave("output/plots/population plots/astro.pdf", width = 4, height = 2.5, units = "in")
# 
# keep <- c("endo.2", "endo.11", "peri.14", "smc.21", "endo.27", "vlmc.30", "endo.38")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic() +
#   theme(legend.position = "none")
# ggsave("output/plots/population plots/endo.pdf", width = 6, height = 6, units = "in")
# 
# keep <- c("opc.19", "opc.37")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("epen.26", "lymph.28")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("pvm.22", "pvm.29")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("micro.3", "micro.4", "micro.25")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("micro.0", "micro.1", "micro.5", "micro.6", "micro.7", "micro.9",
#           "micro.12", "micro.13", "micro.23", "micro.36")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("g.neu.16", "im.neu.18", "im.neu.24", "g.neu.31", "im.neu.32", "g.neu.33", "g.neu.34")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()

#############################################################################
########                        Neuron subclusters                  #########
#############################################################################
md <- readRDS("output/seurat_metadata_neuron_subsets_2.rds")

plt_dt <- md %>% 
  select(orig.ident, Treatment, Region, celltype) %>% 
  count(orig.ident, celltype, Treatment, Region, .drop=FALSE) %>% 
  mutate(percent = 0)
rm(md)

for (i in unique(plt_dt$orig.ident)) {
  print(i)
  total <- sum(plt_dt %>% filter(orig.ident==i) %>% pull(n))
  plt_dt[plt_dt$orig.ident==i,]$percent <- plt_dt[plt_dt$orig.ident==i,]$n/total
}

plt_dt[is.na(plt_dt)]

#Add meta data to missing values
tmp <- readxl::read_excel("library_meta_data.xlsx") %>%
  rename(orig.ident = Library)
plt_dt <- plt_dt %>% left_join(tmp, by = "orig.ident")
plt_dt <- plt_dt %>% select(orig.ident, celltype, Region.y, Treatment.y, n, percent)
plt_dt <- plt_dt %>% rename(Region = Region.y, Treatment = Treatment.y)
rm(tmp)
plt_dt[is.na(plt_dt)]

plt_dt <- plt_dt %>%
  mutate(Treatment.Region = paste(Treatment, Region, sep = "_")) %>%
  mutate(celltype = factor(celltype, levels = c("npc.1", "diff.neu.2", "diff.neu.3", "diff.neu.5",
                                                "diff.neu.12", "diff.neu.14", "cr.neu.7", "m.neu.0",
                                                "m.neu.11", "ex.neu.4", "ex.neu.6", "ex.neu.8",
                                                "ex.neu.9", "ex.neu.15", "inh.neu.10", "inh.neu.13")))

plt_dt <- plt_dt %>%
  mutate(Treatment.Region = factor(Treatment.Region, levels = c("Mock_cx", "RML_cx","Mock_hp", "RML_hp")))

aggregate(percent~celltype+Treatment.Region, plt_dt, mean)

ggplot(plt_dt, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(height=0, width=0.1, alpha = 0.5) +
  #geom_point(data = aggregate(percent~celltype+Treatment.Region, plt_dt, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
  scale_color_manual(values = c("cornflowerblue", "firebrick1", "navy", "firebrick4")) +
  facet_wrap(~celltype, scales="free_y") +
  theme_classic()
ggsave("output/plots/population plots/all_neurons.pdf", width = 8, height = 6, units = "in")

out <- list()
for (i in unique(plt_dt$celltype)) {
  for (j in unique(plt_dt$Region)) {
    mock <- plt_dt %>% filter(celltype == i, Region == j, Treatment == "Mock") %>% pull(percent)
    rml <- plt_dt %>% filter(celltype == i, Region == j, Treatment == "RML") %>% pull(percent)
    
    tmp <- wilcox.test(rml, mock)
    
    out[[paste0(i, j)]] <- data.frame(celltype = i,
                                      Region = j,
                                      pval = tmp$p.value)
  }
}
out <- do.call(rbind, out)
write.csv(out, "output/files/statistics/neuron_subclusters_population_stats.csv")

out %>% filter(pval < 0.12)

tmp <- plt_dt

out <- list()
for (i in unique(tmp$celltype)) {
  print(i)
  out[[i]] <- data.frame(celltype = i,
                         enrichment_hp_cx = (tmp %>% filter(tmp$celltype == i, Region == "hp") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx") %>% pull(percent) %>% mean()),
                         enrichment_hp_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "Mock") %>% pull(percent) %>% mean()),
                         enrichment_cx_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "Mock") %>% pull(percent) %>% mean())
  )
}
out <- do.call(rbind, out)

out$enrichment_hp_cx <- log2(out$enrichment_hp_cx)
out$enrichment_hp_RML_Mock <- log2(out$enrichment_hp_RML_Mock)
out$enrichment_cx_RML_Mock <- log2(out$enrichment_cx_RML_Mock)

write.csv(out, "output/files/statistics/enrichment_neuron_clusters.csv")



# 
# keep <- c("1")
# tmp <- plt_dt %>% filter(celltype %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("5", "3", "2", "12", "14")
# tmp <- plt_dt %>% filter(seurat_clusters %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("7")
# tmp <- plt_dt %>% filter(seurat_clusters %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("0", "11")
# tmp <- plt_dt %>% filter(seurat_clusters %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# keep <- c("4", "8", "10", "13", "15", "9", "6")
# tmp <- plt_dt %>% filter(seurat_clusters %in% keep)
# ggplot(tmp, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
#   geom_jitter(height=0, width=0.1, alpha = 0.5) +
#   geom_point(data = aggregate(percent~celltype+Treatment.Region, tmp, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
#   scale_color_manual(values = c("cornflowerblue", "navy", "firebrick1", "firebrick4")) +
#   facet_wrap(~celltype, scales="free_y") +
#   theme_classic()
# 
# t.test(plt_dt %>% filter(seurat_clusters == "10", Region == "hp", Treatment == "Mock") %>% pull(percent),
#        plt_dt %>% filter(seurat_clusters == "10", Region == "hp", Treatment == "RML") %>% pull(percent))


###########################################################################################
########                             Astrocyte subclusters                        #########
###########################################################################################
md <- readRDS("output/seurat_metadata_astrocyte_subset.rds")

plt_dt <- md %>% 
  select(orig.ident, Treatment, Region, seurat_clusters) %>% 
  count(orig.ident, seurat_clusters, Treatment, Region, .drop=FALSE) %>% 
  mutate(percent = 0)
rm(md)

for (i in unique(plt_dt$orig.ident)) {
  print(i)
  total <- sum(plt_dt %>% filter(orig.ident==i) %>% pull(n))
  plt_dt[plt_dt$orig.ident==i,]$percent <- plt_dt[plt_dt$orig.ident==i,]$n/total
}

plt_dt[is.na(plt_dt)]

#Add meta data to missing values
tmp <- readxl::read_excel("library_meta_data.xlsx") %>%
  rename(orig.ident = Library)
plt_dt <- plt_dt %>% left_join(tmp, by = "orig.ident")
plt_dt <- plt_dt %>% select(orig.ident, seurat_clusters, Region.y, Treatment.y, n, percent)
plt_dt <- plt_dt %>% rename(Region = Region.y, Treatment = Treatment.y)
rm(tmp)
plt_dt[is.na(plt_dt)]

plt_dt <- plt_dt %>%
  mutate(Treatment.Region = paste(Treatment, Region, sep = "_"))
 
plt_dt <- plt_dt %>%
  mutate(Treatment.Region = factor(Treatment.Region, levels = c("Mock_cx", "RML_cx","Mock_hp", "RML_hp")))

aggregate(percent~seurat_clusters+Treatment.Region, plt_dt, mean)
plt_dt$seurat_clusters <- factor(plt_dt$seurat_clusters, levels = c("0", "2", "3", "4", "5", "7",
                                                                    "1", "9", "10",
                                                                    "6", "8"))

ggplot(plt_dt, aes(x=Treatment.Region, y=percent, shape=Treatment.Region, color=Treatment.Region)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(height=0, width=0.1, alpha = 0.5) +
  #geom_point(data = aggregate(percent~celltype+Treatment.Region, plt_dt, mean), aes(x=Treatment.Region, y=percent, color=Treatment.Region), shape=95, size=8) +
  scale_color_manual(values = c("cornflowerblue", "firebrick1", "navy", "firebrick4")) +
  facet_wrap(~seurat_clusters, scales="free_y", nrow=3) +
  theme_classic()
ggsave("output/plots/population plots/all_astro.pdf", width = 8, height = 4.5, units = "in")

out <- list()
for (i in unique(plt_dt$seurat_clusters)) {
  for (j in unique(plt_dt$Region)) {
    mock <- plt_dt %>% filter(seurat_clusters == i, Region == j, Treatment == "Mock") %>% pull(percent)
    rml <- plt_dt %>% filter(seurat_clusters == i, Region == j, Treatment == "RML") %>% pull(percent)
    
    tmp <- wilcox.test(rml, mock)
    
    out[[paste0(i, j)]] <- data.frame(celltype = i,
                                      Region = j,
                                      pval = tmp$p.value)
  }
}
out <- do.call(rbind, out)
write.csv(out, "output/files/statistics/astrocyte_subclusters_population_stats.csv")

out %>% filter(pval < 0.12)

tmp <- plt_dt
tmp$celltype <- tmp$seurat_clusters

out <- list()
for (i in unique(tmp$celltype)) {
  print(i)
  out[[i]] <- data.frame(celltype = i,
                         enrichment_hp_cx = (tmp %>% filter(tmp$celltype == i, Region == "hp") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx") %>% pull(percent) %>% mean()),
                         enrichment_hp_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "hp", Treatment == "Mock") %>% pull(percent) %>% mean()),
                         enrichment_cx_RML_Mock = (tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "RML") %>% pull(percent) %>% mean())/(tmp %>% filter(tmp$celltype == i, Region == "cx", Treatment == "Mock") %>% pull(percent) %>% mean())
  )
}
out <- do.call(rbind, out)

out$enrichment_hp_cx <- log2(out$enrichment_hp_cx)
out$enrichment_hp_RML_Mock <- log2(out$enrichment_hp_RML_Mock)
out$enrichment_cx_RML_Mock <- log2(out$enrichment_cx_RML_Mock)

write.csv(out, "output/files/statistics/enrichment_astrocyte_clusters.csv")
