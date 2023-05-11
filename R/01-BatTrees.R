library(tidyverse)
library(brms) # Running Bayesian models
library(tidybayes) # Pulling results from, and plotting Bayesian models
library(emmeans) # Calculate marginal means and gather marginal mean draws from brms models
library(ggbeeswarm) # Nicer jittered points than geom_jitter
library(pals) # Get palettes with lots of colors
library(patchwork) # Plot multiple plots together
library(taxize) # Get a phylogenetic distance matrix
library(ape) # Calculate the phylogenetic correlation matrix

# Read in the raw data
raw_data <- read_tsv("Bat Differentiation - Bats.tsv") %>% 
  rename(mating_system = `mating system`) %>% 
  mutate(mating_system = str_remove(mating_system, "\\s*\\([^\\)]+\\)")) %>% 
  mutate(migratory = factor(migratory), 
         mating_system = factor(mating_system),
         location = factor(location)) %>% 
  mutate(mating_system = str_replace(mating_system, "congregatory", "swarm"))

ggplot(data = raw_data, aes(x = migratory, y = global_fst)) +
  geom_violin()

ggplot(data = raw_data, aes(x = migratory, y = gene_diversity)) +
  geom_violin()

ggplot(data = raw_data, aes(x = migratory, y = allelic_richness)) +
  geom_violin()

ggplot(data = raw_data, aes(x = global_fst, y = gene_diversity, colour = region)) +
  geom_point() +
  theme_bw()

ggplot(data = raw_data, aes(x = global_fst, y = allelic_richness)) +
  geom_point() +
  theme_bw()

# Read in the metadata for  migration status
migration_status <- read_tsv("bat_life_history.tsv")

# Join the new migration classifications to the genetic data. Also, replace spaces with underscores to keep brms happy
# Turn mating system into a binary column of harem vs. non-harem species, ignoring the ones we do not have data on
combined_data <- left_join(raw_data, migration_status, by = "species") %>% 
  mutate_if(is.character, str_replace_all, ' ', '_') %>% 
  mutate_if(is.character, str_replace_all, '-', '_') %>% 
  mutate(matingsystem2 = case_when(mating_system == "harem" ~ "harem",
                                   mating_system == "swarm" ~ "non-harem",
                                   mating_system == "congregatory" ~ "non-harem")) %>% 
  mutate(species = as.factor(species))


ggplot(data = combined_data, aes(x = region.x, y = allelic_richness)) +
  geom_point() +
  theme_bw()


# Create an object of colors matched up with bat species, to maintain consistency across plots
batColors <- setNames(c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", 
                                 "#005C31", "#2BCE48", "#FFCC99", "#94FFB5", 
                                 "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405", 
                                 "#FFA8BB", "#426600", "#FF0010", "#5EF1F2"), levels(combined_data$species))
                                 
# Create an object of italicized bat species names to maintain consistency across plots
batNames <- setNames(c(bquote(italic("Artibeus jamaicensis")), bquote(italic("Carollia castanea")), bquote(italic("Eptesicus serotinus")), bquote(italic("Lasionycteris noctivagans")), bquote(italic("Lasiurus cinereus")), bquote(italic("Miniopterus schreibersii")), bquote(italic("Myotis blythii")), bquote(italic("Myotis daubentonii")), bquote(italic("Myotis escalerai")), bquote(italic("Myotis lucifugus")), bquote(italic("Myotis myotis")), bquote(italic("Myotis myotis")*"x"*italic("Myotis blythii")), bquote(italic("Myotis septentrionalis")), bquote(italic("Myotis thysanodes")), bquote(italic("Nyctalus lasiopterus")), bquote(italic("Nyctalus leisleri")), bquote(italic("Rhynchonycteris naso")), bquote(italic("Thyroptera tricolor"))), levels(combined_data$species))

########################################################################################################################
# Use taxize to get classifications for the different species. Two useful links are here
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# https://cran.r-project.org/web/packages/taxize/taxize.pdf
# Look up the species names on NCBI
taxize_classification <- classification(as.character(unique(combined_data$species)), db = 'ncbi')
# Create a phylogenetic tree
species_tree <- class2tree(taxize_classification)
plot(species_tree)

# Create a covariance matrix with ape using the phylo object within the species tree
species_cov <- ape::vcv.phylo(species_tree$phylo)
colnames(species_cov) <- make.names(colnames(species_cov), unique=TRUE)
rownames(species_cov) <- make.names(rownames(species_cov), unique=TRUE)

# The species Lasiurus cinereus is listed as Aeorestes cinereus in the NCBI data, so replace that species name with the NCBI one to make the models work
setdiff(combined_data$species_nounderscore, rownames(species_cov))
setdiff(rownames(species_cov), combined_data$species_nounderscore)

# Make the Lasiurus/Aeorestes cinereus name consistent with the NCBI taxonomy
combined_data <- combined_data %>% 
  mutate(species_nounderscore = str_replace(species, "_", ".")) %>% 
  mutate(species_nounderscore = str_replace(species_nounderscore, "Lasiurus.cinereus", "Aeorestes.cinereus"))
#####################################################################################

# What is bwt in terms of mating system?
bwt_mating_model <- brm(data = combined_data, family = gaussian,
                 data2 = list(species_cov = species_cov),
                 bf(global_fst ~ matingsystem2 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)),
                 iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                 control = list(adapt_delta = 0.99, max_treedepth = 16),
                 backend = "cmdstanr")

summary(bwt_mating_model)

# Look at trace plots
plot(bwt_mating_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(bwt_mating_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(bwt_mating_model), points = T, theme = theme_bw(base_size = 18))

get_variables(bwt_mating_model)

# Gather emmeans draws for the model
bwt_mating_emmeans <- gather_emmeans_draws(emmeans(bwt_mating_model, ~ matingsystem2))

# Plot the model with full eye plots and the raw datapoints used in it
bwt_mating_plot <- ggplot(data = bwt_mating_emmeans, aes(x = matingsystem2, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = matingsystem2, y = global_fst, colour = species), varwidth = TRUE, alpha = 0.6) +
  scale_color_manual(values = batColors, labels = batNames) +
  labs(x = element_blank(), y = element_blank(), color = "Species") +
  scale_x_discrete(labels = c("Harem", "Non-Harem")) +
  ylim(-0.1, 0.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
bwt_mating_plot
ggsave(filename = "bwt_mating_plot.png", plot = bwt_mating_plot, dpi = 900, width = 9.5, height = 5.6)



Hs_mating_model <- brm(data = combined_data, family = gaussian,
                       data2 = list(species_cov = species_cov),
                 bf(gene_diversity ~ matingsystem2 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)),
                 iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                 control = list(adapt_delta = 0.99, max_treedepth = 16),
                 backend = "cmdstanr")

summary(Hs_mating_model)

# Look at trace plots
plot(Hs_mating_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(Hs_mating_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(Hs_mating_model), points = T, theme = theme_bw(base_size = 18))
get_variables(Hs_mating_model)

# Gather emmeans draws for the model
Hs_mating_emmeans <- gather_emmeans_draws(emmeans(Hs_mating_model, ~ matingsystem2))

# Plot the model with full eye plots and the raw datapoints used in it
Hs_mating_plot <- ggplot(data = Hs_mating_emmeans, aes(x = matingsystem2, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = matingsystem2, y = gene_diversity, colour = species), varwidth = TRUE, alpha = 0.6) +
  scale_color_manual(values = batColors, labels = batNames) +
  labs(x = element_blank(), y = element_blank(), color = "Species") +
  scale_x_discrete(labels = c("Harem", "Non-Harem")) +
  ylim(0.0, 1.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
Hs_mating_plot
ggsave(filename = "Hs_mating_plot.png", plot = Hs_mating_plot, dpi = 900, width = 9.5, height = 5.6)

# What is the relationship between allelic richness and mating
allelic_richness_mating_model <- brm(data = combined_data, family = gaussian,
                                     data2 = list(species_cov = species_cov),
                bf(allelic_richness ~ matingsystem2 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))),
                prior = c(prior(normal(0, 5), class = Intercept),
                          prior(normal(0, 5), class = b)),
                iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                control = list(adapt_delta = 0.99, max_treedepth = 16),
                backend = "cmdstanr")

summary(allelic_richness_mating_model)

# Look at trace plots
plot(allelic_richness_mating_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(allelic_richness_mating_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(allelic_richness_mating_model), points = T, theme = theme_bw(base_size = 18))

# Gather emmeans draws for the model
allelic_richness_mating_emmeans <- gather_emmeans_draws(emmeans(allelic_richness_mating_model, ~ matingsystem2))

# Plot the model with full eye plots and the raw datapoints used in it
allelic_richness_mating_plot <- ggplot(data = allelic_richness_mating_emmeans, aes(x = matingsystem2, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = matingsystem2, y = allelic_richness, colour = species), varwidth = TRUE, alpha = 0.6) +
  scale_color_manual(values = batColors, labels = batNames) +
  labs(x = "Mating System", y = element_blank(), color = "Species") +
  scale_x_discrete(labels = c("Harem", "Non-Harem")) +
  ylim(0, 10) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
allelic_richness_mating_plot
ggsave(filename = "allelic_richness_mating_plot.png", plot = allelic_richness_mating_plot, dpi = 900, width = 9.5, height = 5.6)

# What is bwt in terms of migration strategy?
bwt_migration_model <- brm(data = combined_data, family = gaussian,
                           data2 = list(species_cov = species_cov),
                        bf(global_fst ~ migratory_class_3 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(normal(0, 1), class = b)),
                        iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                        control = list(adapt_delta = 0.99, max_treedepth = 16),
                        backend = "cmdstanr")

summary(bwt_migration_model)

# Look at trace plots
plot(bwt_migration_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(bwt_migration_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(bwt_migration_model), points = T, theme = theme_bw(base_size = 18))

get_variables(bwt_migration_model)

# Gather emmeans draws for the model
bwt_migration_emmeans <- gather_emmeans_draws(emmeans(bwt_migration_model, ~ migratory_class_3))

# Plot the model with full eye plots and the raw datapoints used in it
bwt_migration_plot <- ggplot(data = bwt_migration_emmeans, aes(x = migratory_class_3, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = migratory_class_3, y = global_fst, colour = species), varwidth = TRUE, alpha = 0.6) +
  #scale_color_brewer(palette = "Paired", labels = c(bquote(italic("Artibeus jamaicensis")), bquote(italic("Carollia castanea")), bquote(italic("Eptesicus serotinus")), bquote(italic("Miniopterus schreibersii")), bquote(italic("Myotis escalerai")), bquote(italic("Myotis lucifugus")), bquote(italic("Myotis septentrionalis")), bquote(italic("Myotis thysanodes")), bquote(italic("Nyctalus lasiopterus")), bquote(italic("Nyctalus leisleri")), bquote(italic("Rhynchonycteris naso")), bquote(italic("Thyroptera tricolor")))) +
  scale_colour_manual(values = batColors, labels = batNames) +
  labs(x = element_blank(), y = bquote(italic("\u03B2")[WT]), color = "Species") +
  scale_x_discrete(limits = c("non_migratory", "regional_migrant", "long_distance_migrant"), labels = c("Non-Migratory", "Regional", "Long Distance")) +
  ylim(-0.1, 0.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
bwt_migration_plot

ggsave(filename = "bwt_migration_plot.png", plot = bwt_migration_plot, dpi = 900, width = 9.5, height = 5.6)


# What is Hs in terms of migration strategy?
Hs_migration_model <- brm(data = combined_data, family = gaussian,
                          data2 = list(species_cov = species_cov),
                           bf(gene_diversity ~ migratory_class_3 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                           prior = c(prior(normal(0, 1), class = Intercept),
                                     prior(normal(0, 1), class = b)),
                           iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                           control = list(adapt_delta = 0.99, max_treedepth = 16),
                           backend = "cmdstanr")

summary(Hs_migration_model)

# Look at trace plots
plot(Hs_migration_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(Hs_migration_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(Hs_migration_model), points = T, theme = theme_bw(base_size = 18))

get_variables(Hs_migration_model)

# Gather emmeans draws for the model
Hs_migration_emmeans <- gather_emmeans_draws(emmeans(Hs_migration_model, ~ migratory_class_3))

# Plot the model with full eye plots and the raw datapoints used in it
Hs_migration_plot <- ggplot(data = Hs_migration_emmeans, aes(x = migratory_class_3, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = migratory_class_3, y = gene_diversity, colour = species), varwidth = TRUE, alpha = 0.6) +
  scale_color_manual(values = batColors, labels = batNames) +
  labs(x = element_blank(), y = bquote("H"[S]), color = "Species") +
  scale_x_discrete(limits = c("non_migratory", "regional_migrant", "long_distance_migrant"), labels = c("Non-Migratory", "Regional", "Long Distance")) +
  ylim(0.0, 1.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
Hs_migration_plot

ggsave(filename = "Hs_migration_plot.png", plot = Hs_migration_plot, dpi = 900, width = 9.5, height = 5.6)


# What is the relationship between allelic richness and migration strategy?
allelic_richness_migration_model <- brm(data = combined_data, family = gaussian,
                                        data2 = list(species_cov = species_cov),
                                     bf(allelic_richness ~ migratory_class_3 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))),
                                     prior = c(prior(normal(0, 5), class = Intercept),
                                               prior(normal(0, 5), class = b)),
                                     iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                                     control = list(adapt_delta = 0.99, max_treedepth = 16),
                                     backend = "cmdstanr")

summary(allelic_richness_migration_model)

# Look at trace plots
plot(allelic_richness_migration_model, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(allelic_richness_migration_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(allelic_richness_migration_model), points = T, theme = theme_bw(base_size = 18))

# Gather emmeans draws for the model
allelic_richness_migration_emmeans <- gather_emmeans_draws(emmeans(allelic_richness_migration_model, ~ migratory_class_3))

# Plot the model with full eye plots and the raw datapoints used in it
allelic_richness_migration_plot <- ggplot(data = allelic_richness_migration_emmeans, aes(x = migratory_class_3, y = .value)) +
  stat_eye(alpha = 0.6) +
  geom_quasirandom(data = combined_data, aes(x = migratory_class_3, y = allelic_richness, colour = species), varwidth = TRUE, alpha = 0.6) +
  scale_color_manual(values = batColors, labels = batNames) +
  labs(x = "Migratory Classification", y = "Allelic Richness", color = "Species") +
  scale_x_discrete(limits = c("non_migratory", "regional_migrant", "long_distance_migrant"), labels = c("Non-Migratory", "Regional", "Long Distance")) +
  ylim(0, 10) +
  theme_bw(base_size = 12)
allelic_richness_migration_plot
ggsave(filename = "allelic_richness_migration_plot.png", plot = allelic_richness_migration_plot, dpi = 900, width = 9.5, height = 5.6)


# Put all the plots together
combined_plot <- (bwt_migration_plot + bwt_mating_plot) / (Hs_migration_plot + Hs_mating_plot) / (allelic_richness_migration_plot + allelic_richness_mating_plot) +
  plot_layout(guides = 'collect')
# Write out all the plots in one figure. Writing them out as a pdf led to a weird issue with the Bwt label being partly covered by the plots, so I was forced to use another format.
ggsave(filename = "combined_plot.png", plot = combined_plot, dpi = 900, height = 7, width = 9)


############################################################################################################
# Assemble various tables as supplemental tables for the manuscript
# A table of sample sizes, first for samples with Hs measures then with Bwt measures
sample_sizes_Hs <- combined_data %>% 
  filter(!is.na(gene_diversity)) %>% 
  group_by(region.x, species) %>% 
  tally()

write_tsv(x = sample_sizes_Hs, file = "sample_sizes_Hs.txt")


sample_sizes_Bwt <- combined_data %>% 
  filter(!is.na(global_fst)) %>% 
  group_by(region.x, species) %>% 
  tally()

write_tsv(x = sample_sizes_Bwt, file = "sample_sizes_Bwt.txt")


# Use emmeans to pull pairwise highest posterior density interval comparisons for each of the 6 models used.
# First, the migratory models
bwt_migration_pairs <- summary(pairs(emmeans(bwt_migration_model, ~ migratory_class_3)))
Hs_migration_pairs <- summary(pairs(emmeans(Hs_migration_model, ~ migratory_class_3)))
allelic_richness_migration_pairs <- summary(pairs(emmeans(allelic_richness_migration_model, ~ migratory_class_3)))

# Write out the results as tables
write_tsv(x = bwt_migration_pairs, file = "bwt_migration_pairs.txt")
write_tsv(x = Hs_migration_pairs, file = "Hs_migration_pairs.txt")
write_tsv(x = allelic_richness_migration_pairs, file = "allelic_richness_migration_pairs.txt")

# Pull results for the mating models
bwt_mating_pairs <- summary(pairs(emmeans(bwt_mating_model, ~ matingsystem2))) %>% 
  mutate(Measure = "Bwt")
Hs_mating_pairs <- summary(pairs(emmeans(Hs_mating_model, ~ matingsystem2))) %>% 
  mutate(Measure = "Hs")
allelic_richness_mating_pairs <- summary(pairs(emmeans(allelic_richness_mating_model, ~ matingsystem2))) %>% 
  mutate(Measure = "Allelic_Richness")

# Combine the different mating model pairwise results because they're just one line each
mating_pairs_combined <- bind_rows(bwt_mating_pairs, Hs_mating_pairs, allelic_richness_mating_pairs)

# Write out the table
write_tsv(x = mating_pairs_combined, file = "mating_pairs_combined.txt")

####################################################################################
# Run intercept-free models to directly compare groups
# Check if the non-harem maters have higher Bwt than the harem ones
bwt_mating_model_NoIntercept <- brm(data = combined_data, family = gaussian,
                        data2 = list(species_cov = species_cov),
                        bf(global_fst ~ 0 + matingsystem2 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                        prior = c(#prior(normal(0, 1), class = Intercept),
                                  prior(normal(0, 1), class = b)),
                        iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                        control = list(adapt_delta = 0.99, max_treedepth = 16),
                        backend = "cmdstanr")

summary(bwt_mating_model_NoIntercept)

# Look at trace plots
plot(bwt_mating_model_NoIntercept, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(bwt_mating_model_NoIntercept, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(bwt_mating_model_NoIntercept), points = T, theme = theme_bw(base_size = 18))

get_variables(bwt_mating_model_NoIntercept)

# What's the chance that the non-harem maters have a higher Bwt than the harem ones?
hypothesis(bwt_mating_model_NoIntercept, "matingsystem2nonMharem > matingsystem2harem")


# Check if Hs is higher in the long distance migrators than the regional and non-migratory ones
Hs_migration_model_NoIntercept <- brm(data = combined_data, family = gaussian,
                          data2 = list(species_cov = species_cov),
                          bf(gene_diversity ~ 0 + migratory_class_3 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                          prior = c(#prior(normal(0, 1), class = Intercept),
                                    prior(normal(0, 1), class = b)),
                          iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                          control = list(adapt_delta = 0.99, max_treedepth = 16),
                          backend = "cmdstanr")

summary(Hs_migration_model_NoIntercept)

# Look at trace plots
plot(Hs_migration_model_NoIntercept, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(Hs_migration_model_NoIntercept, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(Hs_migration_model_NoIntercept), points = T, theme = theme_bw(base_size = 18))

get_variables(Hs_migration_model_NoIntercept)

# What are the chances that Hs is higher in long distance migrators compared to the other two groups?
hypothesis(Hs_migration_model_NoIntercept, c("migratory_class_3long_distance_migrant > migratory_class_3regional_migrant",
                                             "migratory_class_3long_distance_migrant > migratory_class_3non_migratory"))

hypothesis(Hs_migration_model_NoIntercept, c("migratory_class_3long_distance_migrant > migratory_class_3regional_migrant + migratory_class_3non_migratory",
                                             "migratory_class_3long_distance_migrant > migratory_class_3regional_migrant & migratory_class_3non_migratory"))


# Check if allelic richness is higher in the long distance migrators than the regional and non-migratory ones
allelic_richness_migration_model_NoIntercept <- brm(data = combined_data, family = gaussian,
                                      data2 = list(species_cov = species_cov),
                                      bf(allelic_richness ~ 0 + migratory_class_3 + (1 | species) + (1 | gr(species_nounderscore, cov = species_cov))), 
                                      prior = c(#prior(normal(0, 1), class = Intercept),
                                        prior(normal(0, 5), class = b)),
                                      iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                                      control = list(adapt_delta = 0.99, max_treedepth = 16),
                                      backend = "cmdstanr")

summary(allelic_richness_migration_model_NoIntercept)

# Look at trace plots
plot(allelic_richness_migration_model_NoIntercept, N = 4, ask = FALSE)

# Look at posterior predictive checks.
pp_check(allelic_richness_migration_model_NoIntercept, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot conditional effects in the model
plot(conditional_effects(allelic_richness_migration_model_NoIntercept), points = T, theme = theme_bw(base_size = 18))

get_variables(allelic_richness_migration_model_NoIntercept)

# What are the chances that Hs is higher in long distance migrators compared to the other two groups?
hypothesis(allelic_richness_migration_model_NoIntercept, c("migratory_class_3long_distance_migrant > migratory_class_3regional_migrant",
                                                           "migratory_class_3long_distance_migrant > migratory_class_3non_migratory"))


hypothesis(allelic_richness_migration_model_NoIntercept, c("migratory_class_3long_distance_migrant > migratory_class_3regional_migrant + migratory_class_3non_migratory",
                                                           "migratory_class_3long_distance_migrant > migratory_class_3regional_migrant & migratory_class_3non_migratory"))

