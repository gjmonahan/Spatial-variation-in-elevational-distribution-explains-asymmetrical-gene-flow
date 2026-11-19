library("dplyr")
library("ggplot2")
setwd("C:/Users/gabri/Documents/Thom_Lab_Project/DB_Elevational_data")
#bringing in elevational data and pruning it to species for this study
elevational_data_all_species <- read.csv("1_All_sp.csv")
elevational_data_all_species <- elevational_data_all_species[!is.na(elevational_data_all_species$a.elevation), ]
all_species <- unique(elevational_data_all_species$name)
species_to_keep <- c("Castanozoster thoracica", "Drymophila rubricollis","Dysithamnus xanthopterus", "Heliobletus contaminatus", "Hemitriccus obsoletus", "Microspingus lateralis", "Phylloscartes difficilis", "Piculus aurulentus", "Poospiza lateralis", "Stephanophorus diadematus", "Stephanoxis lalandi", "Synallaxis cinerascens", "Syndactyla rufosuperciliata", "Tangara desmaresti", "Saltator maxillosus")
elevational_data_pruned_species <- elevational_data_all_species %>% filter(name %in% species_to_keep)

#elevational_data_pruned_species %>% summarize(avg_val = mean(name, na.rm = TRUE))
#Remove top and bottom 5% values for each species
elevational_data_pruned_species <- elevational_data_pruned_species %>%
  group_by(name) %>%
  filter(a.elevation > quantile(a.elevation, 0.05) & 
           a.elevation < quantile(a.elevation, 0.95)) %>%
  ungroup()
#dividing each species by tropical vs subtropical latitude (-24.5)
elevational_threshold <- -24.5
northern_elevational_data <- elevational_data_pruned_species[elevational_data_pruned_species$latitude >= elevational_threshold, ]
southern_elevational_data <- elevational_data_pruned_species[elevational_data_pruned_species$latitude < elevational_threshold, ]
#finding mean and standard devitation of each species in north and south)
northern_mean_elevation_by_species <- northern_elevational_data %>% group_by(name) %>% summarize(mean_value = mean(a.elevation))
northern_sd_elevation_by_species <- northern_elevational_data %>% group_by(name) %>% summarize(sd_value = sd(a.elevation))
southern_mean_elevation_by_species <- southern_elevational_data %>% group_by(name) %>% summarize(mean_value = mean(a.elevation))
southern_sd_elevation_by_species <- southern_elevational_data %>% group_by(name) %>% summarize(sd_value = sd(a.elevation))
#merging data frames
northern_mean_sd_species <- merge(northern_mean_elevation_by_species, northern_sd_elevation_by_species, by="name")
southern_mean_sd_species <- merge(southern_mean_elevation_by_species, southern_sd_elevation_by_species, by="name")

write.csv(northern_mean_sd_species, "north_elevation_24.csv", row.names = FALSE)
write.csv(southern_mean_sd_species, "south_elevation_24.csv", row.names = FALSE)

boxplot(a.elevation ~ name, data = northern_elevational_data) + coord_flip()

ggplot(northern_elevational_data, aes(x = name, y= a.elevation)) + geom_boxplot() + coord_flip()


N_elevation_only <- select(northern_elevational_data, name, a.elevation)
N_elevation_only <- N_elevation_only %>% mutate(Region = 'North')
S_elevation_only <- select(southern_elevational_data, name, a.elevation)
S_elevation_only <- S_elevation_only %>% mutate(Region = 'South')

species_elevation_localities <- bind_rows(N_elevation_only, S_elevation_only)
#Pruning to make graph look nicer
north_order <- species_elevation_localities %>%
  filter(Region == "North") %>%
  group_by(name) %>%
  summarise(mean_north = mean(a.elevation, na.rm = TRUE)) %>%
  arrange(mean_north) %>%
  pull(name)

species_elevation_localities$name <- factor(species_elevation_localities$name, levels = north_order)

species_elevation_localities$Region <- factor(
  species_elevation_localities$Region,
  levels = c("South", "North")
)

elevational_distribution_figure <- ggplot(species_elevation_localities, aes(x = name, 
                          y = a.elevation, 
                          fill = Region)) + 
       geom_boxplot(position = position_dodge(width = 1.0)) +
  theme_bw() +
  labs(x = "Taxa",
       y = "Elevation",
       fill = "Region") +
  scale_fill_manual(values = c("North" = "#ff5555ff",
                               "South" = "#00aad4ff")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      coord_flip() + stat_summary(fun = mean,              
                                geom = "point",          
                                shape = 21,              
                                size = 2,              
                                color = "black",
                                position = position_dodge(1.0))          

print(elevational_distribution_figure)
ggsave(filename = "elevational_distribution_figure.pdf", plot = elevational_distribution_figure, width = 10, height = 10)

#2 degrees within gap
elevational_data_all_species <- read.csv("1_All_sp.csv")
elevational_data_all_species <- elevational_data_all_species[!is.na(elevational_data_all_species$a.elevation), ]
all_species <- unique(elevational_data_all_species$name)
species_to_keep <- c("Castanozoster thoracica", "Drymophila rubricollis","Dysithamnus xanthopterus", "Heliobletus contaminatus", "Hemitriccus obsoletus", "Microspingus lateralis", "Phylloscartes difficilis", "Piculus aurulentus", "Poospiza lateralis", "Stephanophorus diadematus", "Stephanoxis lalandi", "Synallaxis cinerascens", "Syndactyla rufosuperciliata", "Tangara desmaresti","Saltator maxillosus")
elevational_data_pruned_species <- elevational_data_all_species %>% filter(name %in% species_to_keep)

elevational_data_pruned_species<- elevational_data_pruned_species %>%
  group_by(name) %>%
  filter(a.elevation > quantile(a.elevation, 0.05) & 
           a.elevation < quantile(a.elevation, 0.95)) %>%
  ungroup()

#diving each species by tropical vs subtropical latitude (-24.5)
gap_latitude <- -24.5
north_latitude <- -22.5
south_latitude <- -26.5
north_massif <- elevational_data_pruned_species %>%
  filter(latitude > gap_latitude & latitude < north_latitude)
south_massif <- elevational_data_pruned_species %>%
  filter(latitude < gap_latitude & latitude > south_latitude)

###########
north_massif_mean <- north_massif %>% group_by(name) %>% summarize(mean_value = mean(a.elevation))
north_massif_sd <- north_massif %>% group_by(name) %>% summarize(sd_value = sd(a.elevation))
south_massif_mean <- south_massif %>% group_by(name) %>% summarize(mean_value = mean(a.elevation))
south_massif_sd <- south_massif %>% group_by(name) %>% summarize(sd_value = sd(a.elevation))
#merging data frames
north_massif_mean_sd <- merge(north_massif_mean, north_massif_sd, by="name")
south_massif_mean_sd <- merge(south_massif_mean, south_massif_sd, by="name")


write.csv(north_massif_mean_sd, "north_massif_elevation_23_1_degrees.csv", row.names = FALSE)
write.csv(south_massif_mean_sd, "south_massif_elevation_23_1_degrees.csv", row.names = FALSE)