df_res<-fread("H:/analysis/test/global_virus_IMGVR/IMGVR_all_Sequence_information-high_confidence.tsv",sep="\t",header=T)
df_res<-df_res %>% select(genome=1,votu=6,size=7,Ecosystem= 5, completeness=11,contamination=12,quality= 13,taxonomy=15,provirus=8 ,host=17,source=16);gc()
df_res<-df_res %>% filter(Ecosystem!=";;;")
df_res<-df_res %>%
  separate(col = "Ecosystem", sep=";",into = c("Ecosystem1","Ecosystem2","Ecosystem3","Ecosystem4"))

df2 <- df_res %>%
  mutate(
    Ecosystem2 = case_when(
      Ecosystem2 %in% c("Human", "Mammals: Human") ~ "Human",
      Ecosystem2 %in% c("Aquatic", "Fish", "Algae", "Plankton") ~ "Aquatic",
      str_detect(Ecosystem2,"Arthropoda") ~ "Arthropoda",
      Ecosystem2 %in% c(
        "Mammals", "Birds", "Reptilia", "Amphibia",
        "Mollusca", "Porifera", "Annelida", "Tunicates",
        "Cnidaria", "Cephalochordata", "Nematoda",
        "Invertebrates","Fish"
      ) ~ "Animals",
      Ecosystem2 %in% c("Protists", "Protozoa", "Amoebozoa") ~ "Protists",
      Ecosystem2 %in% c("Lab enrichment", "Lab synthesis",
                        "Laboratory developed", "Lab culture") ~ "Laboratory",
      Ecosystem2 %in% c("Wastewater", "WWTP",
                        "Sewage treatment plant") ~ "Wastewater treatment",
      Ecosystem2 %in% c("Food production", "Animal feed production",
                        "Feedstock") ~ "Food production",
      Ecosystem2 %in% c("", "Unclassified") | is.na(Ecosystem2) ~ "Unclassified",
      TRUE ~ Ecosystem2
    )
  ) %>% filter(!Ecosystem2%in% c("Plankton","Artificial ecosystem","Modeled","Biotransformation","Bioremediation","Microbial","Fungi","Unclassified","Laboratory","Industrial production","Bioreactor")) %>%
  mutate(Ecosystem2=ifelse(Ecosystem3=="Soil","Soil",Ecosystem2)) %>% filter(Ecosystem2!="Terrestrial")
df2<-df2 %>% unite("ecosystem", Ecosystem1:Ecosystem2, sep = ";", na.rm = TRUE)
df2 %>% count(ecosystem)



# summarise
votu_summary <- df2 %>%
  group_by(votu) %>%
  summarise(
    n_genomes = n(),
    n_ecosystems = n_distinct(ecosystem),
    ecosystems = paste(unique(ecosystem), collapse = " | "),
    size= max(size)
  )
summary(votu_summary $n_ecosystems)
table(votu_summary $n_ecosystems)

df <- votu_summary %>% filter(votu !="vOTU_00217951") %>% ##filter phix
  group_by(n_ecosystems) %>%
  summarise(count = n()) %>%
  mutate(count_log = count + 1, count_plus1=count + 1)   # 👈 关键：+1

fit_exp <- lm(log10(count_plus1) ~ n_ecosystems, data = df)
fit_pow <- lm(log10(count_plus1) ~ log10(n_ecosystems), data = df)

get_stats <- function(fit) {
  tibble(
    r2 = summary(fit)$r.squared,
    p  = summary(fit)$coefficients[2,4]
  )
}
stats_exp <- get_stats(fit_exp)
stats_pow <- get_stats(fit_pow)
fmt <- function(r2, p) {
  if (p < 0.001) {
    paste0(
      "atop(R^2 == ", round(r2, 2), ", italic(P) < 0.001)"
    )
  } else {
    paste0(
      "atop(R^2 == ", round(r2, 2), ", italic(P) == ", signif(p, 2), ")"
    )
  }
}
lab_exp <- paste0("atop('Exponential', ", fmt(stats_exp$r2, stats_exp$p), ")")
lab_pow <- paste0("atop('Power-law', ", fmt(stats_pow$r2, stats_pow$p), ")")

df_pred <- tibble(
  n_ecosystems = seq(min(df$n_ecosystems), max(df$n_ecosystems), length.out = 100)
) %>%
  mutate(
    exp = 10^(predict(fit_exp, newdata = .)),
    pow = 10^(predict(fit_pow, newdata = .))
  )

ggplot(df, aes(n_ecosystems, count_plus1)) +
  geom_col(
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_line(data = df_pred, aes(y = exp),
            color = "blue", linetype = "dotted", linewidth = 0.9) +
  geom_line(data = df_pred, aes(y = pow),
            color = "red", linetype = "dashed", linewidth = 0.9) +
  scale_y_log10() +
  labs(
    x = "Extent of viral cross-habitat distribution",
    y = "Count of viral populations (log10 scale)"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = lab_exp,
           hjust = 1.1, vjust = 1.5,
           color = "#4657b7", size = 5,
           parse = TRUE) +
  annotate("text",
           x = Inf, y = Inf,
           label = lab_pow,
           hjust = 1.1, vjust = 3.5,
           color = "#f20c00", size = 5,
           parse = TRUE)+
  scale_y_log10(
    breaks = c(1, 2, 11, 101, 1001,10001,100001,1000001,10000001,100000001),
    labels = c("0", "1",  "10", "100", "1000", "10000", "100000", "1000000", "10000000", "100000000"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_continuous(
    breaks = c(1:max(df$n_ecosystems))
  ) +
  labs(
    x = "Extent of viral cross-habitat distribution",
    y = "Count of viral populations (log10 scale)"
  ) +
  theme(aspect.ratio = 1)+
  coord_cartesian(ylim = c(1,1000001))

ggsave("p_vir_IMGVR_cross_ecosystem_distribution_histogram.png", width = 5, height = 6);rstudioapi::viewer("p_vir_IMGVR_cross_ecosystem_distribution_histogram.png")








##############

library(terra)
data_gis<- rast("H:/analysis/test/新建文件夹/Entropy_01_05_1km_uint16.tif")# as.data.frame(data_gis) %>% head()
df_res<-read.xlsx2("../SJYT_MF2G/mabincollected.xlsx",sheetIndex = 1) %>% mutate(latlon=paste0(lat,lon)) %>%
  filter(Source=="Public Data",!ecosystem%in%c("Artificial Surfaces","Water bodies","Tundra"))  %>%
  mutate(Habitat="Soil",paired_seqs=as.numeric(paired_seqs)) %>% column_to_rownames("sample_name") %>%
  mutate(lon=as.numeric(lon),lat=as.numeric(lat))
df_res <- fread("H:/analysis/test/global_soil/global_soil_viruses_mabin_data_test/sample_data.txt") %>%
  dplyr::select(sample_name= Sample.id,lat=Lat,lon=Lon,ecosystem=Biome)


points <- as.data.frame(cbind(df_res$lon, df_res$lat)) %>% type_convert()
df_res$raster_value <- extract(data_gis, points, method="bilinear")[,2]
df_res$MAT <- extract(clim,points, method="bilinear")[,1]
df_res$MAP <- extract(clim,points, method="bilinear")[,12]

otu<-fread("H:/analysis/test/global_soil/global_soil_viruses_mabin_data_test/OTU_GSV_abundance_table.csv",data.table = F) %>% column_to_rownames("V1")
t_otu<-t(otu)
tmp1<-vegan::diversity(t_otu, index = "shannon")
df_res<-df_res %>% mutate(adiv_shannon=tmp1[sample_name]) %>%
  filter(!ecosystem%in%c("Permanent snow and ice","Cultivated land","Shrubland"))

df_tmp<-df_res %>% filter(!is.na(raster_value)) %>% filter(!ecosystem%in%c("Wetland","Artificial Surfaces","Tundra","Bare Land")) %>% filter(!is.na(adiv_shannon));df_tmp
df_tmp %>% count(ecosystem)

res<-lmer(adiv_shannon ~ raster_value + (1|ecosystem) + MAT + MAP , data = df_tmp )
summary(res)
tmp<-Anova(res);tmp

library(ggeffects)
pred <- ggpredict(res, terms = "raster_value")
head(pred)
col=c("Agricultural Land"="#c49c94","Grassland"="#42b540","Forest"="#ecc37a")

df_tmp %>% #
  mutate(value = adiv_shannon, subtype = ecosystem) %>%
  ggplot(aes(x = raster_value, y = value, color = subtype)) +
  geom_point() +
  geom_ribbon(data = pred, inherit.aes = FALSE, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = pred, inherit.aes = FALSE, aes(x = x, y = predicted), color = "black", linewidth = 1)+
  annotate(geom = "text", label = paste("~italic(P)==",tmp$`Pr(>Chisq)`[1] %>% signif(3)),parse=T ,
           x = -Inf, y = Inf, color = "black",size = 5,hjust = -0.2, vjust = 2.6) +
  annotate(geom = "text", label = paste("n ==",nrow(df_tmp)),parse=T ,
           x = -Inf, y = -Inf, color = "black",size = 5,hjust = -0.2, vjust = -1.8) +
  scale_y_continuous(limits = c(0,NA))+
  scale_color_manual(values=col)+
  labs(x="Habitat diversity",y="Soil viral diversity (Shannon index)",color="Land type")+
  theme(legend.position = "bottom",aspect.ratio = 1)
ggsave("p_virus_shannon_habitat_heterogenity.pdf",width = 5,height = 7);rstudioapi::viewer("p_virus_shannon_habitat_heterogenity.pdf")
