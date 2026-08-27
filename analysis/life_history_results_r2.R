library(tidyverse)
library(vegan)
library(factoextra)
library(ggrepel)
library(ggridges)
library(patchwork)
library(ggforce)
##read in life history outputs
lifehistorypost<-read.csv("analysis/lifehistorypost_r2.csv")
## add posterior draw -- we used 500 samples for each species
n_post<-1000
lifehistorypost$draw<-rep(1:n_post,times=7)

## check that matrices were ergodic and irreducible
table(lifehistorypost$isergodic)
table(lifehistorypost$isirreducible)

# PCA ---------------------------------------------------------------------
## need to pivot trait data from wide to long
## select traits that will go into PCA
lifehistorypost %>% 
  select(draw,species,R0_em,R0_ep) %>% 
  pivot_longer(R0_em:R0_ep,names_to="endo",values_to="R0") %>% 
  mutate(endo=case_when(endo=="R0_em"~"S-",endo=="R0_ep"~"S+"))->R0
lifehistorypost %>% 
  select(draw,species,G_em,G_ep) %>% 
  pivot_longer(G_em:G_ep,names_to="endo",values_to="G") %>% 
  mutate(endo=case_when(endo=="G_em"~"S-",endo=="G_ep"~"S+"))->G
lifehistorypost %>% 
  select(draw,species,meanelexp_em,meanelexp_ep) %>% 
  pivot_longer(meanelexp_em:meanelexp_ep,names_to="endo",values_to="meanelexp") %>% 
  mutate(endo=case_when(endo=="meanelexp_em"~"S-",endo=="meanelexp_ep"~"S+"))->meanelexp
lifehistorypost %>% 
  select(draw,species,matlifexp_em,matlifexp_ep) %>% 
  pivot_longer(matlifexp_em:matlifexp_ep,names_to="endo",values_to="matlifexp") %>% 
  mutate(endo=case_when(endo=="matlifexp_em"~"S-",endo=="matlifexp_ep"~"S+"))->matlifexp
lifehistorypost %>% 
  select(draw,species,longevity_em,longevity_ep) %>% 
  pivot_longer(longevity_em:longevity_ep,names_to="endo",values_to="longevity") %>% 
  mutate(endo=case_when(endo=="longevity_em"~"S-",endo=="longevity_ep"~"S+"))->longevity
lifehistorypost %>% 
  select(draw,species,entropyd_em,entropyd_ep) %>% 
  pivot_longer(entropyd_em:entropyd_ep,names_to="endo",values_to="entropyd") %>% 
  mutate(endo=case_when(endo=="entropyd_em"~"S-",endo=="entropyd_ep"~"S+"))->entropyd
lifehistorypost %>% 
  select(draw,species,firstrepro_em,firstrepro_ep) %>% 
  pivot_longer(firstrepro_em:firstrepro_ep,names_to="endo",values_to="firstrepro") %>% 
  mutate(endo=case_when(endo=="firstrepro_em"~"S-",endo=="firstrepro_ep"~"S+"))->firstrepro
lifehistorypost %>% 
  select(draw,species,shape_surv_em,shape_surv_ep) %>% 
  pivot_longer(shape_surv_em:shape_surv_ep,names_to="endo",values_to="shape_surv") %>% 
  mutate(endo=case_when(endo=="shape_surv_em"~"S-",endo=="shape_surv_ep"~"S+"))->shape_surv
lifehistorypost %>% 
  select(draw,species,shape_rep_em,shape_rep_ep) %>% 
  pivot_longer(shape_rep_em:shape_rep_ep,names_to="endo",values_to="shape_rep") %>% 
  mutate(endo=case_when(endo=="shape_rep_em"~"S-",endo=="shape_rep_ep"~"S+"))->shape_rep

pca.dat<-bind_cols(R0,G$G,meanelexp$meanelexp,longevity$longevity,entropyd$entropyd, firstrepro$firstrepro,shape_surv$shape_surv,shape_rep$shape_rep)
names(pca.dat)<-c("Draw","Species","Endo","R0","GenTime","LifeExpect","Longevity","EntropyD","FirstRepro","ShapeSurv","ShapeRep")

pca.dat %>% select(-Draw) %>% 
  group_by(Species,Endo) %>% 
  summarise_all(mean) -> mean.pca.dat

## reference PCA and loadings based on posterior mean
pca_ref <- prcomp(mean.pca.dat[,-(1:2)], scale. = TRUE)
ref_loadings <- pca_ref$rotation  # this is your fixed basis
ref_center <- pca_ref$center
ref_scale <- pca_ref$scale
summary(pca_ref)

## align each posterior draw to the reference loadings
## calculate shift along PC1 and PC2
posterior_shifts <- list()
for (d in unique(pca.dat$Draw)) {
  draw_df <- pca.dat %>% filter(Draw == d)
  # scale using reference
  mat <- draw_df %>% select(-c(Draw,Species,Endo)) %>% 
    scale(center = ref_center, scale = ref_scale)
  # project onto fixed PCA basis
  scores <- mat %*% ref_loadings[, 1:2]  # project to PC1 and PC2
  scores_df <- bind_cols(draw_df %>% select(Species, Endo), as.data.frame(scores))
  # compute shift vector for each species
  shift_df <- scores_df %>%
    group_by(Species) %>%
    summarise(
      dPC1 = PC1[Endo == "S+"] - PC1[Endo == "S-"],
      dPC2 = PC2[Endo == "S+"] - PC2[Endo == "S-"]
    ) %>%
    mutate(Draw = d)
  posterior_shifts[[length(posterior_shifts)+1]] <- shift_df
}
posterior_shift_df <- bind_rows(posterior_shifts)

posterior_shift_means<-posterior_shift_df %>% 
  group_by(Species) %>% 
  summarise(mean1 = mean(dPC1),
            mean2 = mean(dPC2))%>%
    mutate(direction1 = ifelse(mean1 >= 0, "positive", "negative"),
           direction2 = ifelse(mean2 >= 0, "positive", "negative"))

##Figure of shift in PC space (manuscript Figure 3)
scores_df <- as_tibble(pca_ref$x[, 1:2], .name_repair = "unique") %>%
  mutate(Species = mean.pca.dat$Species,
         Endo = mean.pca.dat$Endo)
arrows_df <- scores_df %>%
  pivot_wider(names_from = Endo, values_from = c(PC1, PC2)) %>% 
  mutate(Species1 = recode(
    Species,
    POAL = "P.al.",
    FESU = "F.s.",
    AGPE = "A.p.",
    POSY = "P.s.",
    POAU = "P.au.",
    ELVI = "E.vir.",
    ELRI = "E.vil."),
    Species2 = recode(
      Species,
      POAL = "Poa alsodes",
      FESU = "Festuca subverticillata",
      AGPE = "Agrostis perennans",
      POSY = "Poa sylvestris",
      POAU = "Poa autumnalis",
      ELVI = "Elymus virginicus",
      ELRI = "Elymus villosus"))
loadings <- as_tibble(pca_ref$rotation[, 1:2], rownames = "Trait")

pca_xlim <- range(scores_df$PC1, na.rm = TRUE)
pca_ylim <- range(scores_df$PC2, na.rm = TRUE)
pca_xlim <- pca_xlim + c(-1, 1) * diff(pca_xlim) * 0.05
pca_ylim <- pca_ylim + c(-1, 1) * diff(pca_ylim) * 0.08

##plot of PC vectors (panel A)
## claude helped with spacing the labels
ggplot() +
  geom_blank(data = scores_df, aes(x = PC1, y = PC2)) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2,
                                    color = Trait), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8) +
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = Trait, color = Trait),
                  size = 4,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,        # always draw a leader line back to the arrow tip
                  segment.size = 0.3,
                  segment.color = "grey50",
                  max.overlaps = Inf,            # never silently drop a label
                  force = 3,                     # push harder against neighbors
                  seed = 42                      # reproducible layout across re-runs
  ) +
  scale_color_brewer(palette = "Dark2") +
  guides(color = "none") +
  labs(title = "A)", x = "PC1 (80.1%)", y = "PC2 (7.8%)") +
  coord_fixed(ratio = 1, xlim = pca_xlim, ylim = pca_ylim, clip = "off") +
  theme_classic() +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5)) -> PCA_vectors

##plot of S+ and S- by species (panel B)
label_adj <- tribble(
  ~Species1, ~nudge_x, ~nudge_y, ~hjust,
  "P.au.",   -0.25,     0.15,     1,
  "P.al.",   -0.25,    -0.05,     1,
  "P.s.",    -0.20,     0.25,     1,
  "A.p.",     0.25,    0.10,     0,
  "F.s.",     0.25,     0.10,     0,
  "E.vir.",   0.00,    0.3,     0.5,
  "E.vil.",   0.25,     0.10,     0
)

arrows_df <- arrows_df %>%
  left_join(label_adj, by = "Species1") %>%
  mutate(PC1_label = `PC1_S-` + nudge_x,
         PC2_label = `PC2_S-` + nudge_y)

ggplot() +
  geom_point(data = scores_df, aes(x = PC1, y = PC2, color = Endo), size = 3) +
  geom_segment(data = arrows_df, aes(x = `PC1_S-`, y = `PC2_S-`, xend = `PC1_S+`, yend = `PC2_S+`),
               linewidth = 0.75, color = "black") +
  geom_text(data = arrows_df, aes(x = PC1_label, y = PC2_label, label = Species1, hjust = hjust),
            size = 4, color = "black", fontface = "italic") +
  scale_color_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue")) +
  guides(color = guide_legend(title = NULL)) +
  labs(title = "B)", x = "PC1 (80.1%)", y = "PC2 (7.8%)") +
  coord_fixed(ratio = 1, xlim = pca_xlim, ylim = pca_ylim, clip = "off") +
  theme_classic() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)) -> PCA_points

## Figure of posterior shifts
lims1 <- quantile(posterior_shift_df$dPC1, probs = c(0.01, 0.99))
posterior_shift_df$Species2<-rep(arrows_df$Species2,n_post)
posterior_shift_df$Species2 <-
  factor(posterior_shift_df$Species2,
         levels = sort(unique(posterior_shift_df$Species2)))
posterior_shift_means$Species2<-arrows_df$Species2
posterior_shift_means$Species2 <-
  factor(posterior_shift_means$Species2,
         levels = levels(posterior_shift_df$Species2))

ggplot(posterior_shift_df, aes(x = dPC1, y = Species2, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black") +
  scale_fill_gradient2(low = "tomato",mid = "white",high = "cornflowerblue",
    midpoint = 0,limits = c(-5, 5),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift")+
  coord_cartesian(xlim = lims1) + 
  theme_minimal(base_size = 13) + 
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),panel.grid = element_blank()) +
  labs(title = "C)",x=expression(Delta~PC1),y=NULL)+ 
  geom_vline(xintercept = 0)+
  geom_point(data = posterior_shift_means,
    aes(x = mean1, y = Species2, color = direction1),  # now fill is mapped!
    inherit.aes = FALSE,shape = 16,size = 2) + 
  scale_color_manual(values = c("positive" = "cornflowerblue", "negative" = "tomato"))-> PC1_plot

lims2 <- quantile(posterior_shift_df$dPC2, probs = c(0.01, 0.99))
ggplot(posterior_shift_df, aes(x = dPC2, y = Species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black") +
  scale_fill_gradient2(low = "tomato",mid = "white",high = "cornflowerblue",
    midpoint = 0,limits = c(-5, 5),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift")+
  coord_cartesian(xlim = lims2) + 
  theme_minimal(base_size = 13) + theme(legend.position = "none") +
  labs(title = "D)",x = expression(Delta~PC2))+ geom_vline(xintercept = 0)+
  geom_point(data = posterior_shift_means,
    aes(x = mean2, y = Species, color = direction2),  # now fill is mapped!
    inherit.aes = FALSE,shape = 16,size = 2)+
  theme(panel.grid = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()) + 
  scale_color_manual(values=c("positive"="cornflowerblue","negative"="tomato"))-> PC2_plot

top_row <- PCA_vectors | PCA_points
bottom_row <- PC1_plot | PC2_plot
PCA_combo <-
  patchwork::wrap_elements(full = top_row,clip = FALSE) /
  patchwork::wrap_elements(full = bottom_row,clip = FALSE) +
  plot_layout(heights = c(-1, 1))
ggsave("manuscript/figures/pca_shift_r2.jpg",
       plot = PCA_combo,
       width = 11,      # width in inches
       height = 7,      # height in inches
       dpi = 300,       # resolution
       units = "in")

## mean shift and probability of PC shift
posterior_shift_df %>% 
  group_by(Species) %>% 
  summarise(mean_shift1=mean(dPC1),
            mean_shift2=mean(dPC2),
            prob_shift1=mean(dPC1>0),
            prob_shift2=mean(dPC2>0))->PC_shifts
## average shift across PC1 vs 2
PC_shifts %>% 
  summarise(mean1 = mean(abs(mean_shift1)), mean2 = mean(abs(mean_shift2)))

# posterior density plots of single traits ---------------------------------
## generate E+/E- contrasts
lifehistorypost$R0_diff<-lifehistorypost$R0_ep-lifehistorypost$R0_em
lifehistorypost$G_diff<-lifehistorypost$G_ep-lifehistorypost$G_em
lifehistorypost$pRep_diff<-lifehistorypost$pRep_ep-lifehistorypost$pRep_em
lifehistorypost$La_diff<-lifehistorypost$La_ep-lifehistorypost$La_em
lifehistorypost$matlifexp_diff<-lifehistorypost$matlifexp_ep-lifehistorypost$matlifexp_em
lifehistorypost$meanelexp_diff<-lifehistorypost$meanelexp_ep-lifehistorypost$meanelexp_em
lifehistorypost$entropyd_diff<-lifehistorypost$entropyd_ep-lifehistorypost$entropyd_em
lifehistorypost$entropyk_diff<-lifehistorypost$entropyk_ep-lifehistorypost$entropyk_em
lifehistorypost$longevity_diff<-lifehistorypost$longevity_ep-lifehistorypost$longevity_em
lifehistorypost$firstrepro_diff<-lifehistorypost$firstrepro_ep-lifehistorypost$firstrepro_em
lifehistorypost$shapesurv_diff<-lifehistorypost$shape_surv_ep-lifehistorypost$shape_surv_em
lifehistorypost$shaperep_diff<-lifehistorypost$shape_rep_ep-lifehistorypost$shape_rep_em
lifehistorypost$lambda_diff<-lifehistorypost$lambda_ep-lifehistorypost$lambda_em
lifehistorypost$lambda_vt_diff<-lifehistorypost$lambda_ep_vt-lifehistorypost$lambda_em_vt

##lambda
# 1) Long format + relabel
lifehistory_long <- lifehistorypost %>%
  pivot_longer(
    cols = c(lambda_diff, lambda_vt_diff),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(parameter = recode(parameter,
                            "lambda_diff" = "Host",
                            "lambda_vt_diff" = "Symbiont"))

# compute posterior means
means_df <- lifehistory_long %>%
  group_by(species, parameter) %>%
  summarise(mu = mean(value), .groups = "drop") %>%
  # offset y values so Host and Symbiont don't overlap
  mutate(ypos = ifelse(parameter == "Host", 0.002, 0.1))

facet_tags <- lifehistory_long %>%
  distinct(species) %>%
  arrange(species) %>%
  mutate(tag = LETTERS[row_number()])

host_symbiont_lambdadiff<-ggplot(lifehistory_long,
       aes(x = value, color = parameter, linetype = parameter)) +
  geom_density(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # posterior mean points (slightly separated in y)
  geom_point(data = means_df,
             aes(x = mu, y = ypos, color = parameter, shape = parameter),
             inherit.aes = FALSE, size = 3) +
  facet_wrap(~ species, scales = "free",
             labeller = as_labeller(c(
               "AGPE" = "italic(Agrostis~pernnans)",
               "ELRI" = "italic(Elymus~villosus)",
               "ELVI" = "italic(Elymus~virginicus)",
               "FESU" = "italic(Festuca~subverticillata)",
               "POAL" = "italic(Poa~alsodes)",
               "POAU" = "italic(Poa~autumnalis)",
               "POSY" = "italic(Poa~sylvestris)"), 
               label_parsed)) +
  theme_classic(base_size = 12) +
  labs(
    x = expression("Fitness difference (" * lambda["S+"] - lambda["S-"] * ")"),
    y = "Probability density",
    color = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme(legend.title = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(12, 12, 12, 24)
  )
ggsave("manuscript/figures/host_symbiont_lambdadiff_r2.jpg",
       plot = host_symbiont_lambdadiff,
       width = 8,      # width in inches
       height = 6,      # height in inches
       dpi = 300,       # resolution
       units = "in")


lifehistorypost %>% 
  group_by(species) %>% 
  summarise(mean_diff = mean(lambda_diff),
            prop_positive = mean(lambda_diff>0),
            prop_negative = mean(lambda_diff<0),
            x_pos = quantile(lambda_diff[lambda_diff > 0], 0.7, na.rm = TRUE),  
            x_neg = quantile(lambda_diff[lambda_diff < 0], 0.3, na.rm = TRUE)) %>%
              mutate(direction = ifelse(mean_diff >= 0, "positive", "negative"))->proportions_df

lambdadiffplot<-ggplot(lifehistorypost, aes(x = lambda_diff, y = species, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, color = "black") +
  scale_fill_gradient2(
    low = "tomato",
    mid = "white",
    high = "cornflowerblue",
    midpoint = 0,
    limits = c(-.65, .65),  # clip values outside this range
    oob = scales::squish,  # squish outliers into endpoint colors
    name = "Shift"
  ) +  
  geom_text(
    data = proportions_df,
    aes(
      x = x_neg,
      y = as.numeric(factor(species)) + 0.13,  # position just below the ridge
      label = paste0(round(prop_negative, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  # Add positive side (right tail) proportions
  geom_text(
    data = proportions_df,
    aes(
      x = x_pos,
      y = as.numeric(factor(species)) + 0.13,  # same vertical offset
      label = paste0(round(prop_positive, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal(base_size = 13) + theme(legend.position = "none") +
  labs(x = expression("Endophyte effect on host fitness ("*Delta~lambda*")"), y = NULL)+ 
  geom_vline(xintercept = 0)+ 
  geom_point(
    data = proportions_df,
    aes(x = mean_diff, y = species, color = direction),  # now fill is mapped!
    inherit.aes = FALSE,
    shape = 16,
    size = 4
  ) +
  scale_color_manual(
    values = c("positive" = "cornflowerblue", "negative" = "tomato")
  )+
  scale_y_discrete(labels = c(
    AGPE = "Agrostis perennans",
    ELRI = "Elymus villosus",
    ELVI = "Elymus virginicus",
    FESU = "Festuca subverticillata",
    POAL = "Poa alsodes",
    POAU = "Poa autumnalis",
    POSY = "Poa sylvestris"))+
  theme(
    axis.text.y = element_text(face = "italic", size = 12),
    panel.grid = element_blank()
  )+xlim(-0.5, 0.8)

ggsave("manuscript/figures/lambdadiffplot_r2.jpg",
       plot = lambdadiffplot,
       width = 6,      # width in inches
       height = 8,      # height in inches
       dpi = 300,       # resolution
       units = "in")


# individual life history traits ------------------------------------------
pca.dat %>% mutate(Species2 = recode(
  Species,
  POAL = "Poa alsodes",
  FESU = "Festuca subverticillata",
  AGPE = "Agrostis perennans",
  POSY = "Poa sylvestris",
  POAU = "Poa autumnalis",
  ELVI = "Elymus virginicus",
  ELRI = "Elymus villosus"))->pca.dat

##shape survival
ggplot(pca.dat,aes(x=ShapeSurv,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  geom_vline(xintercept=0)+facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Shape of survival")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/ShapeSurv_r2.jpg",
        width = 6,height = 6,dpi = 300,units = "in")

##shape reproduction
ggplot(pca.dat,aes(x=ShapeRep,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  geom_vline(xintercept=0)+facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Shape of reproduction")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/ShapeRep_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")
##R0
ggplot(pca.dat,aes(x=R0,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("R0")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/R0_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")

##life expectancy
ggplot(pca.dat,aes(x=LifeExpect,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Mean life expectancy (years)")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/LifeExpect_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")

##max long
ggplot(pca.dat,aes(x=Longevity,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Maximum longevity (years)")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/Longevity_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")

##generation time
ggplot(pca.dat,aes(x=GenTime,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Generation time (years)")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/GenTime_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")

##Entropy
ggplot(pca.dat,aes(x=EntropyD,fill=Endo))+
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",
                 alpha = 0.5,bins = 30)+
  facet_wrap(~Species2,scales="free")+theme_classic()+
  xlab("Demetrius' entropy")+ylab("Posterior probability density")+
  theme(legend.position = c(0.4, 0.1),legend.justification = c("center", "bottom"))+
  scale_fill_manual(values = c("S-" = "tomato", "S+" = "cornflowerblue"))
ggsave(filename="manuscript/figures/EntropyD_r2.jpg",
       width = 6,height = 6,dpi = 300,units = "in")

## data frame for life history effects
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(`R0` = mean(R0_diff>0),
            `Generation time` = mean(G_diff>0),
            `Life expectancy` = mean(meanelexp_diff>0),
            `Longevity` = mean(longevity_diff>0),
            `Entropy` = mean(entropyd_diff>0),
            `Reproductive age` = mean(firstrepro_diff>0),
            `Survival shape` = mean(shapesurv_diff>0),
            `Reproduction shape` = mean(shaperep_diff>0))->life_history_effects

life_history_effects %>%
  pivot_longer(cols = -species, names_to = "metric", values_to = "probpos") %>% 
  ggplot(aes(x = species, y = metric, fill = probpos)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", probpos)), size = 3) +  # Add values to each tile
  scale_fill_gradient2(
    low = "tomato",
    mid = "grey",
    high = "cornflowerblue",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Pr(S+ > S-)"
  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Life history trait")+
  scale_x_discrete(labels = c(
    AGPE = "Agrostis perennans",
    ELRI = "Elymus villosus",
    ELVI = "Elymus virginicus",
    FESU = "Festuca subverticillata",
    POAL = "Poa alsodes",
    POAU = "Poa autumnalis",
    POSY = "Poa sylvestris"))-> LHtraits_heatmap

ggsave("manuscript/figures/LHtraits_heatmap_r2.jpg",
       plot = LHtraits_heatmap,
       width = 10,      # width in inches
       height = 8,      # height in inches
       dpi = 300,       # resolution
       units = "in")

## what is the mean delay in flowering age?
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(mean_change = mean(firstrepro_diff)) %>% 
  summarise(sp_mean_change = mean(mean_change))

## what is that as a percentage increase in age?
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(firstrepro_em = mean(firstrepro_em),
            firstrepro_ep = mean(firstrepro_ep),
            percent_change = (firstrepro_ep - firstrepro_em)/firstrepro_em) %>% 
  summarise(mean(percent_change))

## what is the mean extension of generation time?
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(mean_change = mean(G_diff)) %>% 
  summarise(sp_mean_change = mean(mean_change))

## what is that as a percentage increase in G?
lifehistorypost %>% 
  group_by(species) %>% 
  summarise(G_em = mean(G_em),
            G_ep = mean(G_ep),
            percent_change = (G_ep - G_em)/G_em) %>% 
  summarise(mean(percent_change))
