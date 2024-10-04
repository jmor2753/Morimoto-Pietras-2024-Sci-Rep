#### Morimoto & Pietras (2024) Sci Rep
#### Oct 2024

###### Script for the manuscript Morimoto and Pietras (2024) 
###### NOte*: The analysis are not presented in the order of the manuscript, but they have been signposted here accordingly.
set.seed(12309)
## packages
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)
library(seqinr)
library(ape)
library(phangorn)
library(tidyr)
library(purrr)
library(BiocManager)
library(Biostrings)
library(patchwork)
library(pracma)


## AMino acid datasets
amino_acids_codon <- data.frame(aa_letter = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
                                numb_codons = c(4, 6, 2, 3, 2, 2, 2, 4, 2, 3, 6, 2, 1, 2, 4, 6, 4, 1, 2, 4))

### costs
# Define the data
amino_acids_costs <- data.frame(aa_letter =  c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), 
                                costATP = c(11.7, 24.7, 12.7, 15.3, 52, 11.7, 38.3, 32.3, 30.3, 27.3, 34.3, 14.7, 20.3, 16.3, 27.3, 11.7, 18.7, 23.3, 74.3, 50),
                                decay_invtime = c(1, 30, 9, 5, 4, 1, 14, 2, 8, 2, 13, 10, 3, 8, 4, 6, 6, 2, 12, 7),
                                cost_ATPtime = c(12, 741, 114, 77, 208, 12, 536, 65, 242, 55, 446, 147, 61, 130, 109, 70, 112, 47, 892, 350)) %>%
  right_join(., amino_acids_codon) %>%
  mutate(cost_ATPtime_codon = cost_ATPtime/numb_codons)



### amino acid information
# Your existing data frame
# Creating a single data frame with all columns, including full amino acid names
aa_data <- data.frame(
  'aa_letter' = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                  'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'),
  'aa_name' = c('Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid', 'Phenylalanine', 'Glycine', 'Histidine', 'Isoleucine', 
                'Lysine', 'Leucine', 'Methionine', 'Asparagine', 'Proline', 'Glutamine', 'Arginine', 'Serine', 'Threonine', 
                'Valine', 'Tryptophan', 'Tyrosine'),
  'mol_weight' = c(89.1, 121.2, 133.1, 147.1, 165.2, 75.1, 155.2, 131.2, 146.2, 131.2, 149.2, 132.1, 115.1, 146.2, 174.2, 105.1, 119.1, 117.1, 204.2, 181.2))


# Create data frame
amino_acids <- right_join(aa_data, amino_acids_codon, by = join_by('aa_letter')) %>%
  right_join(., amino_acids_costs)
amino_acids 




















### Analysis 1: Testing if amino acid profiles vary
## Output path
exome_processed_path <- "... your path..."


### Storing all file names
exome_all <- list.files(path = exome_processed_path, 
                        pattern="\\.csv$")

### Listing all files in the folder and pasting the file name
exome_datasets_all <- lapply(paste(exome_processed_path ,
                                   exome_all,
                                   sep = "/"), 
                             read.csv, 
                             header = TRUE)

## joining aa information
exome_all_full <- bind_rows(exome_datasets_all) %>%
  select(-X.1, -X) %>%
  right_join(., amino_acids, by = join_by("aa_letter")) %>%
  mutate()


# getting the total sum of aa and proportions per species
total_count_overall_sp <- exome_all_full %>%
  group_by(species) %>%
  summarise(total_aa_overall_sp = sum(aa_total_count)) %>%
  left_join(., exome_all_full) %>%
  mutate(relativfreq_aa_sp = aa_total_count/total_aa_overall_sp,
         relativfreq_aa_sp_codon = relativfreq_aa_sp/numb_codons)### HERE

### getting the mean per superkingdom
total_count_overall_superkingdom <- total_count_overall_sp %>%
  group_by(superkingdom, aa_letter) %>%
  summarise(se_superkingdom = sd(relativfreq_aa_sp),
            mean_superkingdom = mean(relativfreq_aa_sp),
            se_superkingdom_codon = sd(relativfreq_aa_sp_codon),
            mean_superkingdom_codon = mean(relativfreq_aa_sp_codon)) %>%
  left_join(., total_count_overall_sp)


#### (a) codon vs frequency
codonfreq_model <- lmerTest::lmer(relativfreq_aa_sp ~  
                                    numb_codons*superkingdom+ (1|species), data = total_count_overall_superkingdom)
summary(codonfreq_model)
anova(codonfreq_model, test = "F")

## Codon frequency positively correlated with number of redundant codons


## modelling differences in amino acid frequency
## model mean frequency
mean_freq_aa_superkingdom <- lmerTest::lmer((relativfreq_aa_sp_codon)^0.3 ~ superkingdom * aa_letter + 
                                              (1|species), data = total_count_overall_superkingdom)
anova(mean_freq_aa_superkingdom)
summary(mean_freq_aa_superkingdom)




## Fig 1a

bar_superkingdom_plot <- ggplot(data = total_count_overall_superkingdom, aes(x = reorder(paste(aa_name, paste0("(", aa_letter, ")"), sep = " "),-mean_superkingdom_codon), 
                                                                             y = mean_superkingdom_codon, 
                                                                             fill = superkingdom
)) + 
  geom_bar(stat ="identity", position= position_dodge(), col = "black", linewidth = 0.5) +
  facet_wrap(~superkingdom, ncol = 4) + 
  geom_errorbar(aes(ymin=mean_superkingdom_codon - se_superkingdom_codon, 
                    ymax=mean_superkingdom_codon + se_superkingdom_codon), 
                width=.0, linewidth = 0.5,
                position= position_dodge(0.9)) + 
  theme(legend.position = "bottom") + 
  scale_fill_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  xlab('Amino acid') + 
  ylab('Average frequency') +
  theme_linedraw() + 
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust=0.95,vjust=0.2),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) 
bar_superkingdom_plot







### Analysis 2: Amino acid usage rank diversity 

## EXome analysis _fitting regression to calculate inflection points 

library(mgcv)
library(gratia) # for the derivatives
library(dplyr)
library(tidyr)
library(stringr)
library(msa)
library(Biostrings)
library(purrr)




## ** (2) Obtaining consensus sequence and analyzing proportion of sequences
exome_wide <- total_count_overall_sp %>%
  group_by(species) %>%
  arrange(desc(relativfreq_aa_sp_codon), .by_group = TRUE) %>%
  nest(data = -species) %>%
  mutate(string_pattern_codon = map(data, function(.x){
    as.factor(as.character(AAString(paste(as.character(.x$aa_letter), collapse = ''))))
  })) %>%
  unnest(cols = c(data, string_pattern_codon))

## unique string per species
exome_wide_unique_msa <- exome_wide %>%
  dplyr::select(species, string_pattern_codon, family:superkingdom) %>%
  unique() 

## split the character string and assign value from 1(most frequent) to 20 (least frequent), to each species
cost_aa_df <- exome_wide_unique_msa %>%
  nest(data = -species) %>%
  mutate(aa_order = map(data, function(.x){
    order_aa_list <- data.frame(aa_order = unlist(str_split(.x$string_pattern_codon, "")),
                                numeric_aa_order = 1:20)
  })) %>%
  unnest(c(data, aa_order)) %>%
  right_join(., total_count_overall_superkingdom[c("aa_letter", "relativfreq_aa_sp_codon", "species")], 
             by = join_by("aa_order" == "aa_letter",
                          "species")) 

## Fit a gam model from the most frerq to the least freq, then estimate the gradient (derivative) of this fitted smooth curve,
cost_aa_df_deriv <- cost_aa_df %>%
  nest(data = -species) %>%
  mutate(aa_deriv = map(data, function(.x){
    aa_derivdt <- derivatives(gam(relativfreq_aa_sp_codon ~ s(numeric_aa_order), data = .x), order = 1)
  })) 

## loop that unlist the above data. This was not possible without a loop. This loop is slow but runs.
full_cost_aa_df <- data.frame()
for(i in 1:nrow(cost_aa_df_deriv)){
  int_aa_loop <- cost_aa_df_deriv$aa_deriv[[i]]
  int_aa_loop$species <- cost_aa_df_deriv$species[i]
  full_cost_aa_df <- bind_rows(full_cost_aa_df,int_aa_loop)
} # slow: about 5 to 10 minutes to run in my laptop. Sorry but no other way. Happy to take suggestions.

## filter all aa that overlap zero
zero_sequence_sp <- full_cost_aa_df %>%
  mutate(round_zeros = round(numeric_aa_order)) %>%
  dplyr::select(species, round_zeros) %>%
  unique() 

### here we are analysing the proportions of each amino acid in each of the positions
test_zeros <- full_cost_aa_df %>%
  mutate(round_zeros = round(numeric_aa_order)) %>%
  dplyr::select(species, round_zeros) %>%
  unique() %>%
  right_join(., cost_aa_df, 
             by = join_by("round_zeros" == "numeric_aa_order",
                          "species")) %>%
  unique()



test_zeros_sum <- test_zeros %>%
  ungroup() %>%
  group_by(aa_order, round_zeros, superkingdom) %>%
  tally() %>%
  arrange(desc(superkingdom))


divergent_palette <- colorRampPalette(c("khaki", "grey60", "turquoise", "orange1", "orchid2", "royalblue", "firebrick4", "green2"))(20)

prop_aa_position <- ggplot(test_zeros_sum, aes(fill=aa_order, y=n, x=round_zeros)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~superkingdom, scales = "free_y") + 
  scale_fill_manual("Amino acid", values = divergent_palette) +
  xlab('Rank\n(most to least freq)') + 
  ylab('Proportion') +
  theme_linedraw() + 
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) +
  scale_x_continuous(breaks = seq(1, 20, 2)) + 
  guides(fill=guide_legend(ncol=2)) 
prop_aa_position




## model
prop_all_model <- glm(n ~ round_zeros*aa_order*superkingdom, data = test_zeros_sum, family = "quasipoisson")
prop_all_model <- glm(cbind(n, 20) ~ round_zeros*aa_order*superkingdom, data = test_zeros_sum, family = "quasibinomial")
summary(prop_all_model)
anova(prop_all_model, test= "F")





#### Edge effect global
### must calculate diversity and similarity per site.

diveristy_aa <- test_zeros_sum %>% group_by(superkingdom, round_zeros) %>%
  unique() %>%
  summarise(relativ_diversity = n()/20,
            absolute_diversity = n())



diversity_matrix_wide <- test_zeros_sum %>% ungroup() %>%
  group_by(superkingdom, round_zeros) %>%
  unique() %>%
  pivot_wider(values_from = n,
              names_from = aa_order,
              values_fill = 0) %>%
  mutate(superkingdom = as.factor(superkingdom))


## estimating alpha diversity by rank
div_metrics2 <- diversity_matrix_wide %>%
  ungroup() %>%
  nest(data = c(-superkingdom)) %>%
  dplyr::mutate(alpha = map(data, function(.x){
    shannon <- data.frame(position = 1:20,
                          value = tabula::heterogeneity(.x %>%
                                                          arrange(round_zeros) %>%
                                                          dplyr::select(-round_zeros), method = "shannon")@.Data,
                          index = "Shannon",
                          type = "heterogeneity")
    
    count <- data.frame(position = 1:20,
                        value = tabula::richness(.x %>%
                                                   arrange(round_zeros) %>%
                                                   dplyr::select(-round_zeros), method = "count")@.Data,
                        index = "Count",
                        type = "richness")
    dplyr::bind_rows(shannon,
                     count)
  })) %>%
  tidyr::unnest(cols = c(alpha, superkingdom))



index_plot <- ggplot(data = div_metrics2, 
                     aes(x = position, y = value, 
                         col = superkingdom,
                         group = superkingdom)) +
  geom_point() + 
  geom_line() + 
  facet_grid(index~., scales = "free_y") + 
  scale_x_continuous(breaks=seq(1, 20, 2)) + 
  theme(legend.position = "right") + 
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  xlab('Rank') + 
  ylab('Amino acid diversity') +
  theme_linedraw() + 
  theme(legend.position = "bottom", 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(vjust = 0),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) 

index_plot 



## model diversity 
model_div_proteome_Shannon <- lm(value ~ position*superkingdom + 
                                   I(position^2)*superkingdom, data = div_metrics2 %>%
                                   filter(index == "Shannon"))

summary(model_div_proteome_Shannon)
anova(model_div_proteome_Shannon)

model_div_proteome_Count<- lm(value ~ position*superkingdom + 
                                I(position^2)*superkingdom, data = div_metrics2 %>%
                                filter(index == "Count"))

summary(model_div_proteome_Count)
anova(model_div_proteome_Count)







zero_cost_total <- cost_aa_df %>%
  dplyr::select(species:numeric_aa_order, relativfreq_aa_sp_codon) %>%
  inner_join(., zero_sequence_sp, by = join_by("species", "numeric_aa_order" == "round_zeros")) 


tableS1 <- zero_cost_total %>% 
  dplyr::select(superkingdom, phylum, 
                class, order, species, string_position = numeric_aa_order, 
                string_pattern_codon, aa_letter = aa_order,
                relative_freq_codon = relativfreq_aa_sp_codon)





### Analysis 3: Testing temperature effects on the proteome

curveType_analysis_superkingdom <- total_count_overall_superkingdom %>%
  dplyr::select(superkingdom, species, aa_letter, relativfreq_aa_sp_codon) %>%
  unique() %>%
  arrange(desc(relativfreq_aa_sp_codon)) %>%
  group_by(species) %>%
  mutate(rank = row_number(relativfreq_aa_sp_codon)) %>%
  mutate(rank = abs(20-rank) + 1) %>%
  arrange(desc(relativfreq_aa_sp_codon))




### Effect of growth temperature
#### loading extreme halophilic bacteria data 
extremedb <- read.csv("...yoour path ",
                      header = TRUE,
                      strip.white = TRUE)

extremedb_exome <- curveType_analysis_superkingdom[which(curveType_analysis_superkingdom$species %in% extremedb$Species),]



thermodb <- read.csv("your path",
                     header = TRUE,
                     strip.white = TRUE,
                     na.strings = "")

## data from Paper:https://doi.org/10.1016/j.jtbi.2016.08.011
lehmanndt <- data.frame(
  species = c("Arabidopsis thaliana", "Caenorhabditis elegans", "Galus gallus", "Drosophila melanogaster", "Homo sapiens", "Mus musculus", "Saccharomyces cerevisiae"),
  avg_optTemp = c(21.5, 25, 41.2, 25, 37, 37, 30),
  superkingdom = rep("Eukaryota", 7)
)

## joining lehmann et al data
thermodb <- full_join(thermodb, lehmanndt)
## removing rows with all NA entries
thermodb_exome <- curveType_analysis_superkingdom[which(curveType_analysis_superkingdom$species %in% thermodb$species),] %>% na.omit()
thermodb_selected <- thermodb[which(thermodb$species %in% curveType_analysis_superkingdom$species),] %>%
  dplyr::select(species, ecosystem, environment, oxygen_type, avg_optTemp) %>%
  filter(!duplicated(.))


thermodb_final <- thermodb_exome %>%
  inner_join(., thermodb_selected, by = join_by("species"), multiple = "first", unmatched = "drop")


## discretising
thermodb_final$temp_d <- cut(thermodb_final$avg_optTemp, breaks = c(0, 70,110))


##
thermodb_final_prop <- thermodb_final %>%
  group_by(superkingdom, aa_letter, rank, temp_d) %>%
  tally() %>%
  na.omit()


thermo_prop_aa <- ggplot(thermodb_final_prop, aes(fill=aa_letter, y=n, x=rank)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(superkingdom~temp_d) + 
  scale_fill_manual("Aminoacid", values = divergent_palette) +
  xlab('Rank\n(most to least freq)') + 
  ylab('Proportion') +
  theme_linedraw() + 
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "grey80")) +
  scale_x_continuous(breaks = seq(1, 20, 2))  +
  guides(fill=guide_legend(ncol=2)) 

thermo_prop_aa




### must calculate diversity and similarity per site.
thermodb_final_matrix_wide_all <- thermodb_final_prop %>%
  group_by(superkingdom, rank, temp_d, aa_letter) %>%
  summarise(n = n()) %>%
  unique() %>%
  group_by(superkingdom, temp_d) %>%
  unique() %>%
  pivot_wider(values_from = n,
              names_from = aa_letter,
              values_fill = 0) %>%
  mutate(superkingdom = as.factor(superkingdom)) %>%
  ungroup()




thermodb_final %>%
  select(species, temp_d) %>%
  unique() %>%
  group_by(temp_d) %>%
  summarise(n())



thermodb_diversity_metrics <-  thermodb_final_matrix_wide_all %>%
  ungroup() %>%
  nest(data = -c(superkingdom, temp_d)) %>%
  dplyr::mutate(alpha_class = map(data, function(.x){
    shannon <- data.frame(position = 1:20,
                          value = tabula::heterogeneity(.x %>%
                                                          arrange(rank) %>%
                                                          dplyr::select(-rank), method = "shannon")@.Data,
                          index = "Shannon",
                          type = "heterogeneity")
    count <- data.frame(position = 1:20,
                        value = tabula::richness(.x %>%
                                                   arrange(rank) %>%
                                                   dplyr::select(-rank), method = "count")@.Data,
                        index = "Count",
                        type = "richness")
    dplyr::bind_rows(shannon,
                     count)
  })) %>%
  tidyr::unnest(cols = c(alpha_class, superkingdom, temp_d))



thermodb_index_plot_class <- ggplot(data = thermodb_diversity_metrics, 
                                    aes(x = position, y = value, 
                                        col = superkingdom,
                                        group = superkingdom)) +
  geom_point() + 
  geom_line() + 
  facet_grid(index~temp_d, scales = "free_y") + 
  scale_x_continuous(breaks=seq(1, 20, 2)) + 
  theme(legend.position = "right") + 
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  xlab('Rank') + 
  ylab('Index value') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(vjust = 1),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "grey80"))

SupplFile4 <- thermo_prop_aa/thermodb_index_plot_class 
SupplFile4








### model diversity temp
model_div_temperature <- lm(value ~ I(position^2)*temp_d*superkingdom, data = subset(thermodb_diversity_metrics,
                                                                                     index == "Shannon"))
anova(model_div_temperature )


model_div_temperatureCount <- lm(value ~ I(position^2)*temp_d*superkingdom, data = subset(thermodb_diversity_metrics,
                                                                                          index == "Count"))
anova(model_div_temperatureCount)




## Model
thermodb_finaldf <- thermodb_final %>% dplyr::select(relativfreq_aa_sp_codon, 
                                                     avg_optTemp,
                                                     superkingdom,
                                                     aa_letter,
                                                     species) %>%
  na.omit()

thermodb_model <- lmerTest::lmer(log(relativfreq_aa_sp_codon) ~ log(avg_optTemp)*superkingdom*aa_letter +
                                   (1|species), data = thermodb_finaldf)

thermodb_model_notemp <- lmerTest::lmer(log(relativfreq_aa_sp_codon) ~ superkingdom*aa_letter +
                                          (1|species), data = thermodb_finaldf)

stats::anova(thermodb_model, thermodb_model_notemp )

MuMIn::r.squaredGLMM(thermodb_model)
MuMIn::r.squaredGLMM(thermodb_model_notemp)


anova(thermodb_model, test = "F")


thermodb_final_plot_df <- right_join(thermodb_final, amino_acids, by = join_by("aa_letter")) %>%
  mutate(aa_name_symbol = paste(aa_name,
                                paste0("(", aa_letter, ")")))

Fig1b <- ggplot(thermodb_final_plot_df, aes(x = log(avg_optTemp), y = log(relativfreq_aa_sp_codon),fill = superkingdom, col = superkingdom)) +
  facet_wrap(~aa_name_symbol, ncol = 5)+
  geom_point(alpha = 0.08, size =1.5, pch =21) + 
  geom_smooth(method = "lm", se = FALSE, formula = 'y~poly(x, 1)', fullrange = FALSE) +
  theme(legend.position = "right") + 
  scale_fill_manual('Domain:', values = c("darkolivegreen1", "pink1", "royalblue1")) +
  scale_colour_manual('Domain:', values = c("olivedrab4", "orchid3", "royalblue3")) +
  xlab('Average optimal growth temperature (ln)') + 
  ylab('Average frequency (ln)') +
  theme_linedraw() + 
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 11, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) 

Fig1b












##### Analysis 4: PDB analysis and secondary structures

pdb_output_full <- read.csv("your path",
                            header = TRUE,
                            strip.white = TRUE)
sse_unique_species <- unique(pdb_output_full$species)


## adding taxonomy for those we know
#### Taxonomy of PDB output
simplify_freq_aa <- total_count_overall_sp %>%
  dplyr::select(species, aa = aa_letter, aa_full, relativefreq_aa, superkingdom, numb_codons, relativfreq_aa_sp, relativfreq_aa_sp_codon)


## assigning taxonomy for those that already exist in the data
sse_full_dt1 <- inner_join(pdb_output_full, simplify_freq_aa, by = c("aa", "species"), keep = FALSE)

## checking whatever is left
sse_full_missing <- anti_join(pdb_output_full, simplify_freq_aa, by = c("aa", "species"))

sse_full_missing$superkingdom <- NA

## Viral ones
virus_insex <- str_detect(sse_full_missing$species, "virus")
sse_full_missing$superkingdom[virus_insex] <- "Viruses"
sse_full_viruses <- sse_full_missing %>% filter(superkingdom == "Viruses")

### joining viruses with the rest of the data
sse_full_dt2 <- bind_rows(sse_full_dt1, sse_full_viruses)





## taxonomy for the rest
sse_full_missing_2 <- sse_full_missing %>%
  filter(is.na(superkingdom))

sse_species_tax <- unique(sse_full_missing_2$species)



sse_species_tax2 <- ifelse(str_count(sse_species_tax, "\\w+") > 2, 
                           word(sse_species_tax, 1, 2),
                           sse_species_tax)

sse_species_tax2 <- unique(sse_species_tax2)



## combining to create a final version
tax_out_clean <- read.csv("your path",
                          header = TRUE)





sse_full_missing_2$species <- ifelse(str_count(sse_full_missing_2$species, "\\w+") > 2, 
                                     word(sse_full_missing_2$species, 1, 2),
                                     sse_full_missing_2$species)


sse_full_dt3_final <- sse_full_missing_2 %>%
  select(-superkingdom) %>%
  inner_join(., tax_out_clean, by = c("species"), relationship = "many-to-many") %>%
  select(aa, count, pdb, sse, species, total_aa, aa_prop, superkingdom) %>%
  bind_rows(., sse_full_dt2 %>%
              select(aa, count, pdb, sse, species, total_aa, aa_prop, superkingdom)) %>%
  unique()



sse_full_dt4_final <- sse_full_dt3_final %>% 
  full_join(., amino_acids, by = join_by("aa" == "aa_letter")) %>%
  mutate(aa_prop_sse_codon = aa_prop/numb_codons)


### calculating proportions
sse_prop_aa_superkingdom <- sse_full_dt4_final %>%
  group_by(aa, superkingdom, sse) %>%
  summarise(mean_sse_aa_codon = mean(aa_prop_sse_codon),
            sd_sse_aa_codon = sd(aa_prop_sse_codon),
            mean_sse_aa = mean(aa_prop),
            sd_sse_aa = sd(aa_prop)) %>%
  ungroup()





### Creating plot by rank
rank_sse_analysis_superkingdom <- sse_full_dt4_final %>%
  dplyr::select(superkingdom, species, aa, sse) %>%
  unique() %>% 
  right_join(., sse_full_dt4_final) %>%
  dplyr::select(superkingdom, species, aa, sse, aa_prop_sse_codon) %>%
  group_by(species, sse, aa) %>%
  mutate(mean_sse_codon_aaType = mean(aa_prop_sse_codon)) %>%
  dplyr::select(-aa_prop_sse_codon) %>%
  unique() %>%
  arrange(desc(mean_sse_codon_aaType), .by_group = TRUE) %>%
  ungroup() %>%
  group_by(species, sse) %>%
  mutate(rank = row_number(mean_sse_codon_aaType)) %>%
  mutate(rank = abs(20-rank) + 1) %>%
  group_by(species) %>%
  arrange(desc(mean_sse_codon_aaType)) %>%
  ungroup()



prop_sse_aa <- rank_sse_analysis_superkingdom %>%
  ungroup() %>%
  group_by(aa, rank, superkingdom, sse) %>%
  tally() %>%
  arrange(desc(superkingdom))


prop_aa_sse_table <- prop_sse_aa %>% 
  pivot_wider(names_from = rank,
              values_from = n,
              values_fill = 0) 

chisq.test(as.matrix(prop_aa_sse_table %>%
                       filter(sse == "Beta") %>%
                       ungroup() %>%
                       select(`1`:`20`)),simulate.p.value = TRUE, B = 500)


## testing that amino acid prop in each sse differ
prop_treatment_test_sseAlpha <- chisq.test(as.matrix(prop_aa_sse_table %>%
                                                       filter(sse == "Alpha") %>%
                                                       ungroup() %>%
                                                       select(`1`:`20`)),simulate.p.value = TRUE, B = 5000)

prop_treatment_test_sseBeta <- chisq.test(as.matrix(prop_aa_sse_table %>%
                                                      filter(sse == "Beta") %>%
                                                      ungroup() %>%
                                                      select(`1`:`20`)),simulate.p.value = TRUE, B = 5000)

prop_treatment_test_sseAlpha
prop_treatment_test_sseBeta







### DIVERSITY
prop_aa_sse_table_diversity <- prop_sse_aa %>% 
  pivot_wider(names_from = aa,
              values_from = n,
              values_fill = 0) 


sse_diversity_rank_df <-  prop_aa_sse_table_diversity %>%
  ungroup() %>%
  nest(data = -c(superkingdom, sse,)) %>%
  dplyr::mutate(alpha_sse = map(data, function(.x){
    shannon <- data.frame(position = 1:20,
                          value = tabula::heterogeneity(.x %>%
                                                          arrange(rank) %>%
                                                          dplyr::select(-rank), method = "shannon")@.Data,
                          index = "Shannon",
                          type = "heterogeneity")
    count <- data.frame(position = 1:20,
                        value = tabula::richness(.x %>%
                                                   arrange(rank) %>%
                                                   dplyr::select(-rank), method = "count")@.Data,
                        index = "Count",
                        type = "richness")
    dplyr::bind_rows(shannon,
                     count)
  })) %>%
  tidyr::unnest(cols = c(alpha_sse, superkingdom, sse))


Fig3b <- ggplot(data = sse_diversity_rank_df, 
                aes(x = position, y = value, 
                    col = superkingdom,
                    group = superkingdom)) +
  geom_point() + 
  geom_line() + 
  facet_grid(index~sse, scales = "free_y") + 
  scale_x_continuous(breaks=seq(1, 20, 2)) + 
  theme(legend.position = "right") + 
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  xlab('Rank') + 
  ylab('Amino acid diversity') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(vjust = 1, size = 10),
        axis.text.y = element_text(vjust = 1, size = 10),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) 

Fig3b



## model SSE
model_sse_diversity <- lm(value ~ I(position^2)*sse*superkingdom, data = subset(sse_diversity_rank_df, index == "Shannon"))
anova(model_sse_diversity)
summary(model_sse_diversity)


model_sse_diversityCount <- lm(value ~ I(position^2)*sse*superkingdom, data = subset(sse_diversity_rank_df, index == "Count"))
anova(model_sse_diversityCount)
summary(model_sse_diversityCount)








## relationship frequency and sse
sse_full_dt5_final <- sse_full_dt4_final %>%
  select(superkingdom, aa, sse, species, aa_prop_sse_codon) %>%
  right_join(., total_count_overall_sp %>%
               select(species, aa_letter, superkingdom, relativfreq_aa_sp_codon),
             by = join_by("aa" == "aa_letter",
                          "superkingdom", 
                          "species")) %>%
  na.omit() %>%
  group_by(superkingdom, aa, sse) %>%
  summarise(mean_freq = mean(relativfreq_aa_sp_codon),
            sd_freq = sd(relativfreq_aa_sp_codon),
            mean_sse = mean(aa_prop_sse_codon),
            sd_sse = sd(aa_prop_sse_codon))





## relationship between average freq in the proteome and in sse
Fig3a<- ggplot(sse_full_dt5_final,
               aes(x = mean_freq,
                   y = mean_sse,
                   col = superkingdom)) + 
  facet_grid(sse~superkingdom) +
  geom_text(aes(label = aa)) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab('Average frequency\n(Proteome)') + 
  ylab('Average frequency\n(SSE)') +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "lightblue1"),
        legend.position = "right")

Fig3a


## models
alpha_model <- lm(mean_freq ~ mean_sse*superkingdom, data = subset(sse_full_dt5_final, sse == "Alpha"))
anova(alpha_model)
summary(alpha_model)

beta_model <- lm(mean_freq ~ mean_sse*superkingdom, data = subset(sse_full_dt5_final, sse == "Beta"))
anova(beta_model)
summary(beta_model)
##












##### ** TEsting for the relationship between frequency propensity calculated by CChou et al and frequency observed here
### Amino acid propensity
# Define the data

### SOurce: Chou PY, Fasman GD: Conformational parameters for amino acids in helical, beta-sheet, and random coil regions calculated from proteins. Biochemistry 1974, 13(2):211–222. 10.1021/bi00699a001

# Create the data frame
propensity_aa_sse <- data.frame(amino_acids = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", 
                                                "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"),
                                aa_letter = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", 
                                              "F", "P", "S", "T", "W", "Y", "V"),
                                freq_alpha = c(0.522, 0.282, 0.263, 0.351, 0.278, 0.421, 0.549, 0.19, 0.446, 0.358, 
                                               0.48, 0.383, 0.429, 0.402, 0.212, 0.282, 0.295, 0.409, 0.22, 0.409),
                                freq_beta = c(0.167, 0.154, 0.113, 0.137, 0.222, 0.211, 0.044, 0.138, 0.122, 0.274, 
                                              0.209, 0.126, 0.286, 0.219, 0.62, 0.124, 0.205, 0.203, 0.22, 0.282))



propensity_aa_sse_long <- propensity_aa_sse %>%
  pivot_longer(cols = freq_alpha:freq_beta,
               names_to = "SSE",
               values_to = "Propensity")

propensity_aa_sse_long$SSE <- ifelse(propensity_aa_sse_long$SSE == "freq_alpha", "Alpha", "Beta")


SupplFile3 <- ggplot(propensity_aa_sse_long, aes(x = Propensity, col = SSE, fill = SSE)) +
  geom_density(alpha = 0.2) + 
  facet_grid(~SSE) +
  theme_linedraw() + 
  theme(panel.grid = element_blank()) +
  scale_fill_manual('Secondary structure', 
                    values = c("orange", "turquoise1")) +
  scale_colour_manual('Secondary structure', 
                      values = c("orange", "turquoise2"))
SupplFile3










### Simulations 
set.seed(57849)

## numb of rounds
R <- 100

sim_sse <- c()
for(w in 1:R) {
  for(i in seq(6, 66, by = 8)){
    alpha_sim <- data.frame(aa = sample(propensity_aa_sse$aa_letter, size = i,
                                        prob = propensity_aa_sse$freq_alpha,
                                        replace = TRUE),
                            sse = "Alpha",
                            iteration = w,
                            size = i)
    beta_sim <- data.frame(aa = sample(propensity_aa_sse$aa_letter, size = i,
                                       prob = propensity_aa_sse$freq_beta,
                                       replace = TRUE),
                           sse = "Beta",
                           iteration = w,
                           size = i)
    sim_sse <- bind_rows(sim_sse,
                         alpha_sim,
                         beta_sim)
  }
  
}



sim_sse_prop <- sim_sse %>% group_by(aa, iteration, sse, size) %>%
  tally() %>%
  ungroup() %>%
  group_by(iteration, sse, size) %>%
  mutate(rank = row_number(n)) %>%
  mutate(rank = abs(20-rank) + 1) %>%
  ungroup() %>%
  arrange(rank)




### diversity of simulated data
sim_sse_prop <- sim_sse %>% group_by(aa, size, sse) %>%
  tally() %>%
  ungroup() %>%
  group_by(sse, size) %>%
  mutate(rank = row_number(n)) %>%
  mutate(rank = abs(20-rank) + 1) %>%
  ungroup() %>%
  arrange(rank)


### assembling with different proportions
### Loop for simulating virtual proteins with varying ratios of sse and sizes of sse
sim_sse_prop_nested <- sim_sse_prop %>%
  nest(data = -c(sse,size))


sim_mixture <- c()
sim_mixture_internal <- c()
for(s in 1:30){
  for(i in seq(2, 60, by = 20)){
    mix1090 <- data.frame(sse = sample(unique(sim_sse_prop_nested$sse), size = i, prob = c(0.1, 0.9), replace = TRUE),
                          mixture = "10/90",
                          rep = s,
                          length_sim = i)
    mix5050 <- data.frame(sse = sample(unique(sim_sse_prop_nested$sse), size = i, prob = c(0.5, 0.5), replace = TRUE),
                          mixture = "50/50",
                          rep = s,
                          length_sim = i)
    mix7525 <- data.frame(sse = sample(unique(sim_sse_prop_nested$sse), size = i, prob = c(0.75, 0.25), replace = TRUE),
                          mixture = "75/25",
                          rep = s,
                          length_sim = i)
    mix2575 <- data.frame(sse = sample(unique(sim_sse_prop_nested$sse), size = i, prob = c(0.25, 0.75), replace = TRUE),
                          mixture = "25/75",
                          rep = s,
                          length_sim = i)
    mix9010 <- data.frame(sse = sample(unique(sim_sse_prop_nested$sse), size = i, prob = c(0.1, 0.9), replace = TRUE),
                          mixture = "90/10",
                          rep = s,
                          length_sim = i)
    sim_mixture_internal <- bind_rows(sim_mixture_internal,
                                      mix5050,
                                      mix7525,
                                      mix2575,
                                      mix1090,
                                      mix9010)
  }
  sim_mixture <- bind_rows(sim_mixture, sim_mixture_internal)
}




sim_mixture_prop <- sim_mixture %>% right_join(., sim_sse_prop_nested %>%
                                                 unnest(cols = c(size, sse, data)),
                                               relationship = "many-to-many")


sim_mixture_prop_final <- sim_mixture_prop %>%
  mutate(length_d = cut(length_sim, breaks=c(0, 20, 40, 70),
                        labels = c("Short: 2-20", "Medium: 21-40", "Long: 41-66"))) %>%
  group_by(mixture, rep, aa, rank, length_d) %>%
  tally()


sim_prop_aa_position_sse_SIMULATED <- ggplot(sim_mixture_prop_final, aes(fill=aa, y=n, x=rank)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(length_d~mixture) + 
  scale_fill_manual("Aminoacid", values = divergent_palette) +
  xlab('Rank\n(most to least freq)') + 
  ylab('Proportion') +
  theme_linedraw() + 
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) +
  scale_x_continuous(breaks = seq(1, 20, 2))

sim_prop_aa_position_sse_SIMULATED




sim_mixture_prop_final_table <- sim_mixture_prop_final %>%
  group_by(mixture, rep, rank, aa, length_d) %>%
  pivot_wider(names_from = aa,
              values_from = n,
              values_fill = 0) %>%
  ungroup()



## diversity sse
sse_diversity_rank_df_simulated_mixtures <-  sim_mixture_prop_final_table %>%
  ungroup() %>%
  nest(data = -c(mixture, rep, length_d)) %>%
  dplyr::mutate(alpha_sse_sim = map(data, function(.x){
    shannon <- data.frame(position = sort(.x$rank),
                          value = tabula::heterogeneity(.x %>%
                                                          arrange(rank) %>%
                                                          dplyr::select(-rank), method = "shannon")@.Data,
                          index = "Shannon",
                          type = "heterogeneity")
    count <- data.frame(position = sort(.x$rank),
                        value = tabula::richness(.x %>%
                                                   arrange(rank) %>%
                                                   dplyr::select(-rank), method = "count")@.Data,
                        index = "Count",
                        type = "richness")
    dplyr::bind_rows(shannon,
                     count)
  })) %>%
  tidyr::unnest(cols = c(alpha_sse_sim, rep, length_d))


sse_diversity_rank_df_simulated_mixtures_means <- sse_diversity_rank_df_simulated_mixtures %>%
  #mutate(length_d = cut(length_sim, breaks=c(0, 20, 40, 70),
  #    labels = c("short", "medium", "long"))) %>%
  group_by(position, index, length_d, mixture) %>%
  summarise(mean_value = mean(value),
            sd_value = sd(value))


index_plot_sse_simulation_SIMULATED <- ggplot(data = sse_diversity_rank_df_simulated_mixtures_means, 
                                              aes(x = position, y = mean_value, col = length_d)) +
  geom_point() + 
  geom_errorbar(aes(ymax = mean_value+sd_value, ymin = mean_value-sd_value), width = 0.1) + 
  geom_line() + 
  facet_grid(index~mixture, scales = "free_y") + 
  scale_x_continuous(breaks=seq(1, 20, 2)) + 
  theme(legend.position = "right") + 
  scale_colour_manual('Length SSE:', values = c("royalblue1", "goldenrod2", "orangered1"),
                      breaks = c("Short: 2-20", "Medium: 21-40", "Long: 41-66"),
                      labels = c("Short: 2-20", "Medium: 21-40", "Long: 41-66")) +
  xlab('Rank') + 
  ylab('Amino acid diversity') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(vjust = 1, size = 10),
        axis.text.y = element_text(vjust = 1, size = 10),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) 

index_plot_sse_simulation_SIMULATED



### model
sse_diversity_rank_df_simulated_mixtures_modeldt <- sse_diversity_rank_df_simulated_mixtures %>%
  ungroup() %>%
  select(mixture, rep, length_d, position, value, index)



modelsim_data_diversity <- lm(value ~ I(position^2)*length_d*mixture, data = subset(sse_diversity_rank_df_simulated_mixtures_modeldt,
                                                                                    index == "Shannon"))

anova(modelsim_data_diversity)
summary(modelsim_data_diversity)


modelsim_data_diversityCount <- lm(value ~ I(position^2)*length_d*mixture, data = subset(sse_diversity_rank_df_simulated_mixtures_modeldt,
                                                                                         index == "Count"))
anova(modelsim_data_diversityCount)
summary(modelsim_data_diversity)









#### ANalysis 5: comparing observed vs simulated edge effect
### Comparing diversity metrics between observed and simulated data

### observed diversity
div_metrics_comparing_dt <- div_metrics2 %>%
  select(-data, -type) %>%
  mutate(type = "Observed",
         mean_value = value) %>%
  select(-value)





div_metrics_comparing_dt2 <- replicate(25, div_metrics_comparing_dt, simplify = FALSE) %>%
  imap_dfr(~ .x %>% 
             mutate(rep = .y))


## simulated diversity



### by each position
mean_position_div_sim <- sse_diversity_rank_df_simulated_mixtures %>%
  ungroup() %>%
  group_by(rep, position, index) %>%
  summarise(mean_value = mean(value)) %>%
  ungroup() %>%
  mutate(type = "Simulated (Propensity)")



diversity_compare_obsexp_df <- div_metrics_comparing_dt2 %>%
  select(-superkingdom) %>%
  ungroup() %>%
  bind_rows(., mean_position_div_sim) %>%
  group_by(position, index, type) %>%
  summarise(mean_val = mean(mean_value),
            sd_val = sd(mean_value))

ggplot(diversity_compare_obsexp_df , aes(x = position, y = mean_val, col = type)) + 
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax =mean_val + sd_val, ymin = mean_val - sd_val), width = 0.1) + 
  facet_grid(index~., scale ="free_y") +
  scale_colour_manual('Type:', values = c("black", "purple", "firebrick3"),
                      breaks = c("Observed", "Simulated (Propensity)", "Simulated (from data)"),
                      labels = c("Observed", "Simulated (Propensity)", "Simulated (Data)")) +
  xlab('Rank') + 
  ylab('Amino acid diversity') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(vjust = 1, size = 10),
        axis.text.y = element_text(vjust = 1, size = 10),
        strip.text = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "lightblue1")) +
  scale_x_continuous(breaks=seq(1, 20, 1))



set.seed(1420)
model_diversity_Shannon <- lm(mean_val ~I(position^2)*type, data = diversity_compare_obsexp_df %>%
                                filter(index == "Shannon"))
anova(model_diversity_Shannon)
summary(model_diversity_Shannon)

model_diversity_Count <- lm(mean_val ~  
                              I(position^2)*type, data = diversity_compare_obsexp_df %>%
                              filter(index == "Count"))
anova(model_diversity_Count)
summary(model_diversity_Count)











######### expected based on evolutionary order
# https://pdf.sciencedirectassets.com/271116/1-s2.0-S0378111900X01722/1-s2.0-S0378111900004765/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOn%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQDGjDum9JdTdEvosauDhn3g1Tlf0PuvmbgpcWj4cI98TwIhANt63QqZFXbjHkQWxPWaNJjUVd6kFsx6KcbKns8BDtfXKrMFCDIQBRoMMDU5MDAzNTQ2ODY1IgxHM3SSjgeZ9x1PsQ4qkAX1WCPFFOs5I7e%2BbkFhhrZVQEq%2FAiMxXF9C5MNaD9Zi%2FG9WM%2BTf7zKyCLU4J4rMjOmuAIkSPlWyUloKmLyM90ALSoKGt56m8nUaoz0H%2B2NOY36JHg46VCSBti8kJTcQeLKn8qU7uPNqShi6cwfC64O7NRjX6yYMbzSfr48qEn1UCkRUuUZr74lEJeEr4w7vtSIP4ef%2BcBdqKm%2FjiwpqOK2IvHyCQFR3av8kNQDnaXzSf1zYKJAJihcb57A9eOmKyNKkwwJ6jSEzkmlqVDIrq0RJBqh1564PgEZU%2BVlqVgOqaFb%2Fv4PE%2FwmG5EzXAFc2gQGKZ7n5kPIpm0ZaHv9%2F4cAjMiiaXdUnE%2B3%2Fz3w43we%2FANy49Bf6pPK2ed9AE232T6Pr3AV%2F%2FGqEcsmgUTpeRMdthiYQdP9M8wHsNxdWdeZ%2FxtQUrbc5z193041%2B7aXZ0gIaQ%2FyorVMAEOsTT01AY88wRu3bRsG1Kj8DRSrfDVQC4A3CSASti70DrOZbTPOlmlRVIO9I5L4yVC03jfIEXFCj2Jv30iynRh8igJ8raGW4n%2F1SkECZMl3iHvfXKVbZOKSH%2BWR64tUyEfiCTN4E%2B5y6ClpmePopVo96u7uO5zpPcnqE52zLnHpC36aXK9uF62PpIq%2BVfrkGkup0v4Lrv4VCSegeRWVPtc4MTVbCE%2Bx%2BVak9zY908QAhkh882uyEmZPFfJQMoEwCbI2DoA0endcLeFiZU4A9pNytNKp8WSIG4gFf9SeouTzZxCpXNDiIUwEe8aLK13US4GaiighCqiDPnJTUALUh%2BD0ZFk0ODuRteeQahpsNIfr6wPGGHA7n8stNN%2FvWv2GEZbmwzNXZuKHVSVgF7uZTXvHZ3CXRaVc9RTDj4I%2BxBjqwAXbdS6uAXplMNQ9Fe66bPdGUYZCg%2B2I6QmUrNrRrM%2FlDyrb%2FOVL1ofhZqAAb5CA6zqiNafu0gaZgmjVH7zoB9oJiuJqSHwpibLUDw4NRSIJjiBRRZlq%2BVfIPuYpbWKq8x3XNmm5yBqodlMZ9WawbIDQq%2BhOcgKlxNjiVcG0ViFF9YLpa329vqHR4CC4pgKVREh4NskFDUDbBknID3XrTlXdhi1yCqjU1mIEZALLwobdV&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240420T174837Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYTOORVMXW%2F20240420%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=6f8fd7a35193037a179e85634240f2c9670136718474ef0f3a2f014aa3c22a64&hash=a7f4f05ed4d27ff42b787024e80818373fb3ba59b7ad9e271f33c883f369c96f&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0378111900004765&tid=spdf-fa99ecd0-8bfe-40a0-9305-5018e80b17f6&sid=147df87c7143c842e78b45074e7ebbec0466gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=02015758505008585604&rr=8776f75ef9883533&cc=de
aa_evolorder <- data.frame(
  aa_letter = c("G", "A", "V", "D", "S", "E", "P", "L", "T", "I", "N", "R", "K", "Q", "C", "F", "H", "M", "Y", "W"),
  rawdata_evol = c(4.4, 4.9, 6.9, 7.2, 7.9, 8.2, 8.3, 9.4, 10.1, 11.2, 11.8, 12.0, 12.0, 12.4, 12.4, 13.0, 13.3, 14.0, 14.7, 15.8),
  filtdata_evol = c(2.9, 2.9, 6.6, 7.0, 7.2, 7.5, 7.7, 9.5, 9.8, 11.5, 12.2, 12.3, 13.0, 13.0, 14.3, 14.9, 15.1, 15.4, 15.6, 16.7)
)

aa_evolorder




average_rank <- test_zeros_sum %>%
  ungroup() %>%
  select(aa_letter = aa_order,
         rank = round_zeros, 
         superkingdom, 
         n) %>%
  group_by(superkingdom, rank, aa_letter) %>%
  summarise(times_rank = sum(n)) %>%
  ungroup() %>%
  group_by(superkingdom, aa_letter) %>%
  summarise(mean_rank = weighted.mean(rank, w = times_rank)) %>%
  ungroup()
average_rank


evol_hist_fulldf <- average_rank %>%
  right_join(., aa_evolorder, by = "aa_letter")

# Models
rawdataevol <- lm(rawdata_evol~ mean_rank*superkingdom, data = evol_hist_fulldf)
summary(rawdataevol)
anova(rawdataevol)
confint(rawdataevol)



filtereddataevol <- lm(filtdata_evol ~ mean_rank*superkingdom, data = evol_hist_fulldf)
summary(filtereddataevol)
anova(filtereddataevol)
confint(filtereddataevol)




### PLotting

## relationship between average freq in the proteome and in sse
Fig4a<- ggplot(evol_hist_fulldf,
               aes(x = mean_rank,
                   y = rawdata_evol,
                   col = superkingdom)) + 
  facet_grid(~superkingdom) +
  geom_text(aes(label = aa_letter)) +
  geom_smooth(method = "lm", se = TRUE) +
  xlab('Average rank') + 
  ylab('Raw order') +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "lightblue1"),
        legend.position = "none")

Fig4a


## relationship between average freq in the proteome and in sse
Fig4b<- ggplot(evol_hist_fulldf,
               aes(x = mean_rank,
                   y = filtdata_evol,
                   col = superkingdom)) + 
  facet_grid(~superkingdom) +
  geom_text(aes(label = aa_letter)) +
  geom_smooth(method = "lm", se = TRUE) +
  xlab('Average rank') + 
  ylab('Filtered order') +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "lightblue1"),
        legend.position = "right")

Fig4b


wrap_plots(Fig4a, Fig4b) + plot_annotation(tag_levels = 'a', tag_suffix = ".")







### by superkingdom
## raw data 
rawdataevol_Archaea <- lm(rawdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Archaea"))
summary(rawdataevol_Archaea)
confint(rawdataevol_Archaea)



rawdataevol_Bacteria <- lm(rawdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Bacteria"))
summary(rawdataevol_Bacteria)
confint(rawdataevol_Bacteria)


rawdataevol_Eukaryota <- lm(rawdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Eukaryota"))
summary(rawdataevol_Eukaryota)
confint(rawdataevol_Eukaryota)


rawdataevol_Viruses <- lm(rawdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Viruses"))
summary(rawdataevol_Viruses)
confint(rawdataevol_Viruses)



### filtered data
filtereddataevol_Archaea <- lm(filtdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Archaea"))
summary(filtereddataevol_Archaea)
confint(filtereddataevol_Archaea)



filtereddataevol_Bacteria <- lm(filtdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Bacteria"))
summary(filtereddataevol_Bacteria)
confint(filtereddataevol_Bacteria)


filtereddataevol_Eukaryota <- lm(filtdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Eukaryota"))
summary(filtereddataevol_Eukaryota)
confint(filtereddataevol_Eukaryota)


filtereddataevol_Viruses <- lm(filtdata_evol~ mean_rank, data = subset(evol_hist_fulldf, superkingdom == "Viruses"))
summary(filtereddataevol_Viruses)
confint(filtereddataevol_Viruses)
















## Principal Component Analysis
PCAdata_sp <- total_count_overall_sp %>%
  select(species, class, superkingdom, aa_full, relativfreq_aa_sp_codon) %>%
  pivot_wider(names_from = aa_full,
              values_from = relativfreq_aa_sp_codon)


# install.packages("devtools")
devtools::install_github("arleyc/PCAtest")

PCAspecies <- prcomp(PCAdata_sp[4:23])
summary(PCAspecies)

library(PCAtest)
PCAspecies_stats <- PCAtest(PCAdata_sp[4:23])


PCAdata_sp_final <- as.data.frame(cbind(PCAdata_sp[1:3],
                                        PCAspecies$x))

PCAdata_sp_final_centroids <- 
  PCAdata_sp_final %>%
  select(PC1, PC2, superkingdom) %>%
  group_by(superkingdom) %>%
  summarise(cPC1 = mean(PC1),
            cPC2 = mean(PC2)) %>%
  ungroup()



## estimating centroid standard error
out_PCA_replacement <- data.frame()
for(i in 1:1000){
  calc <- PCAdata_sp_final %>% group_by(superkingdom) %>%
    sample_n(dim(.)[1], replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    group_by(superkingdom) %>%
    summarise(cPC1 = mean(PC1),
              cPC2 = mean(PC2)) %>%
    ungroup() %>%
    mutate(rep = i)
  
  out_PCA_replacement <- rbind(out_PCA_replacement, calc)
  
  
}
out_PCA_replacement_summarised <- out_PCA_replacement %>%
  group_by(superkingdom) %>%
  summarise(mean_cPC1 = mean(cPC1),
            mean_cPC2 = mean(cPC2),
            sd_cPC1 = sd(cPC1),
            sd_cPC2 = sd(cPC2))

library(ggConvexHull)
PC12_plot <- ggplot(PCAdata_sp_final, aes(x = PC1, 
                                          y = PC2, colour = superkingdom)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = out_PCA_replacement_summarised, mapping = aes(x = mean_cPC1, y = mean_cPC2, fill = superkingdom), pch = 22, size = 4, col = "black") +
  #geom_errorbar(data = out_PCA_replacement_summarised, mapping = aes(y = mean_cPC2, x = mean_cPC1, ymin = mean_cPC2 - sd_cPC2, ymax = mean_cPC2 + sd_cPC2), size = 2,
  #              inherit.aes = FALSE) + 
  geom_convexhull(aes(fill = superkingdom, colour = superkingdom), alpha = 0.1) +
  #facet_wrap(~superkingdom, ncol = 2)
  xlab('PC1 (71.88 %)') + 
  ylab('PC2 (10.32 %)') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "grey80"),
        legend.position = "right") +
  scale_fill_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1"))
PC12_plot



### Using hausdorff distance
library(pracma)
## estimating centroid standard error
out_PCA_replacement_HD <- data.frame()
for(i in 1:500){
  calc <- PCAdata_sp_final %>% group_by(superkingdom) %>%
    sample_n(100, replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    ungroup()
  
  
  BA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  EA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  VV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BB_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  HD_estimates <- rbind(data.frame(HD = BA_comp, Comparison = "Bacteria-Archaea", Type = "Between", rep = i),
                        data.frame(HD = BV_comp, Comparison = "Bacteria-Viruses", Type = "Between", rep = i),
                        data.frame(HD = BE_comp, Comparison = "Bacteria-Eukaryota",Type = "Between", rep = i),
                        data.frame(HD = EA_comp, Comparison = "Eukaryota-Archaea", Type = "Between", rep = i),
                        data.frame(HD = EV_comp, Comparison = "Eukaryota-Viruses", Type = "Between", rep = i),
                        data.frame(HD = AV_comp, Comparison = "Archaea-Viruses", Type = "Between", rep = i),
                        data.frame(HD = AA_comp, Comparison = "Archaea-Archaea", Type = "Within", rep = i),
                        data.frame(HD = BB_comp, Comparison = "Bacteria-Bacteria", Type = "Within", rep = i),
                        data.frame(HD = EE_comp, Comparison = "Eukaryota-Eukaryota", Type = "Within", rep = i),
                        data.frame(HD = VV_comp, Comparison = "Viruses-Viruses", Type = "Within", rep = i))
  
  
  out_PCA_replacement_HD <- rbind(out_PCA_replacement_HD, HD_estimates)
  
  
}




HD_Plot <- ggplot(data = out_PCA_replacement_HD %>% filter(Type == "Between"),
                  aes(x = HD, fill = Comparison), col = "black") + 
  geom_histogram(bins = 45) +
  xlab("Hausdorff distance") + 
  scale_fill_brewer(palette = "Dark2")+ 
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black")) + 
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", col = "red") + 
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", col = "black")

HD_Plot 
wrap_plots(PC12_plot, HD_Plot) + plot_annotation(tag_levels = 'a', tag_suffix = ".")

















### PCA Analysis without standardisation (Supplementary material)
PCAdata_sp_raw <- total_count_overall_sp %>%
  select(species, class, superkingdom, aa_full, relativfreq_aa_sp) %>%
  pivot_wider(names_from = aa_full,
              values_from = relativfreq_aa_sp)

## Principal Component Analysis
# install.packages("devtools")


PCAspecies_raw <- prcomp(PCAdata_sp_raw[4:23])
summary(PCAspecies_raw)



PCAdata_sp_final_raw <- as.data.frame(cbind(PCAdata_sp_raw[1:3],
                                            PCAspecies_raw$x))

PCAdata_sp_final_centroids_raw <- 
  PCAdata_sp_final_raw %>%
  select(PC1, PC2, superkingdom) %>%
  group_by(superkingdom) %>%
  summarise(cPC1 = mean(PC1),
            cPC2 = mean(PC2)) %>%
  ungroup()



## estimating centroid standard error
out_PCA_replacement_raw <- data.frame()
for(i in 1:1000){
  calc_raw <- PCAdata_sp_final_raw %>% group_by(superkingdom) %>%
    sample_n(dim(.)[1], replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    group_by(superkingdom) %>%
    summarise(cPC1 = mean(PC1),
              cPC2 = mean(PC2)) %>%
    ungroup() %>%
    mutate(rep = i)
  
  out_PCA_replacement_raw <- rbind(out_PCA_replacement_raw, calc_raw)
  
  
}
out_PCA_replacement_summarised_raw <- out_PCA_replacement_raw %>%
  group_by(superkingdom) %>%
  summarise(mean_cPC1 = mean(cPC1),
            mean_cPC2 = mean(cPC2),
            sd_cPC1 = sd(cPC1),
            sd_cPC2 = sd(cPC2))

library(ggConvexHull)
PC12_plot_raw <- ggplot(PCAdata_sp_final_raw, aes(x = PC1, 
                                                  y = PC2, colour = superkingdom)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = out_PCA_replacement_summarised_raw, mapping = aes(x = mean_cPC1, y = mean_cPC2, fill = superkingdom), pch = 22, size = 4, col = "black") +
  #geom_errorbar(data = out_PCA_replacement_summarised, mapping = aes(y = mean_cPC2, x = mean_cPC1, ymin = mean_cPC2 - sd_cPC2, ymax = mean_cPC2 + sd_cPC2), size = 2,
  #              inherit.aes = FALSE) + 
  geom_convexhull(aes(fill = superkingdom, colour = superkingdom), alpha = 0.1) +
  #facet_wrap(~superkingdom, ncol = 2)
  xlab('PC1 (71.88 %)') + 
  ylab('PC2 (10.32 %)') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "grey80"),
        legend.position = "right") +
  scale_fill_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1"))
PC12_plot_raw



### Using hausdorff sitance
## estimating centroid standard error
out_PCA_replacement_HD_raw <- data.frame()
for(i in 1:500){
  calc<- PCAdata_sp_final_raw %>% group_by(superkingdom) %>%
    sample_n(100, replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    ungroup()
  
  
  BA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  EA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  VV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BB_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  HD_estimates_raw <- rbind(data.frame(HD = BA_comp, Comparison = "Bacteria-Archaea", Type = "Between", rep = i),
                            data.frame(HD = BV_comp, Comparison = "Bacteria-Viruses", Type = "Between", rep = i),
                            data.frame(HD = BE_comp, Comparison = "Bacteria-Eukaryota",Type = "Between", rep = i),
                            data.frame(HD = EA_comp, Comparison = "Eukaryota-Archaea", Type = "Between", rep = i),
                            data.frame(HD = EV_comp, Comparison = "Eukaryota-Viruses", Type = "Between", rep = i),
                            data.frame(HD = AV_comp, Comparison = "Archaea-Viruses", Type = "Between", rep = i),
                            data.frame(HD = AA_comp, Comparison = "Archaea-Archaea", Type = "Within", rep = i),
                            data.frame(HD = BB_comp, Comparison = "Bacteria-Bacteria", Type = "Within", rep = i),
                            data.frame(HD = EE_comp, Comparison = "Eukaryota-Eukaryota", Type = "Within", rep = i),
                            data.frame(HD = VV_comp, Comparison = "Viruses-Viruses", Type = "Within", rep = i))
  
  
  out_PCA_replacement_HD_raw <- rbind(out_PCA_replacement_HD_raw, HD_estimates_raw)
  
  
}
library(pracma)

HD_Plot_raw <- ggplot(data = out_PCA_replacement_HD_raw %>% filter(Type == "Between"),
                      aes(x = HD, fill = Comparison), col = "black") + 
  geom_histogram(bins = 45) +
  xlab("Hausdorff distance") + 
  scale_fill_brewer(palette = "Dark2")+ 
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black")) + 
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", col = "red") + 
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", col = "black")

HD_Plot_raw 
final_plot_raw <- wrap_plots(PC12_plot_raw, HD_Plot_raw) + plot_annotation(title = "No standardisation")













### PCA Analysis standardizing by Molecular Weight (Supplementary material)
PCAdata_sp_MW <- total_count_overall_sp %>%
  mutate(aa_mw = relativfreq_aa_sp/mol_weight) %>%
  select(species, class, superkingdom, aa_full, aa_mw) %>%
  pivot_wider(names_from = aa_full,
              values_from = aa_mw)

## Principal Component Analysis
# install.packages("devtools")


PCAspecies_MW <- prcomp(PCAdata_sp_MW[4:23])
summary(PCAspecies_MW)



PCAdata_sp_final_MW <- as.data.frame(cbind(PCAdata_sp_MW[1:3],
                                           PCAspecies_MW$x))

PCAdata_sp_final_centroids_MW <- 
  PCAdata_sp_final_MW %>%
  select(PC1, PC2, superkingdom) %>%
  group_by(superkingdom) %>%
  summarise(cPC1 = mean(PC1),
            cPC2 = mean(PC2)) %>%
  ungroup()



## estimating centroid standard error
out_PCA_replacement_MW <- data.frame()
for(i in 1:1000){
  calc_MW <- PCAdata_sp_final_MW %>% group_by(superkingdom) %>%
    sample_n(dim(.)[1], replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    group_by(superkingdom) %>%
    summarise(cPC1 = mean(PC1),
              cPC2 = mean(PC2)) %>%
    ungroup() %>%
    mutate(rep = i)
  
  out_PCA_replacement_MW <- rbind(out_PCA_replacement_MW, calc_MW)
  
  
}
out_PCA_replacement_summarised_MW <- out_PCA_replacement_MW %>%
  group_by(superkingdom) %>%
  summarise(mean_cPC1 = mean(cPC1),
            mean_cPC2 = mean(cPC2),
            sd_cPC1 = sd(cPC1),
            sd_cPC2 = sd(cPC2))


PC12_plot_MW <- ggplot(PCAdata_sp_final_MW, aes(x = PC1, 
                                                y = PC2, colour = superkingdom)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = out_PCA_replacement_summarised_MW, mapping = aes(x = mean_cPC1, y = mean_cPC2, fill = superkingdom), pch = 22, size = 4, col = "black") +
  #geom_errorbar(data = out_PCA_replacement_summarised, mapping = aes(y = mean_cPC2, x = mean_cPC1, ymin = mean_cPC2 - sd_cPC2, ymax = mean_cPC2 + sd_cPC2), size = 2,
  #              inherit.aes = FALSE) + 
  geom_convexhull(aes(fill = superkingdom, colour = superkingdom), alpha = 0.1) +
  #facet_wrap(~superkingdom, ncol = 2)
  xlab('PC1 (71.88 %)') + 
  ylab('PC2 (10.32 %)') +
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "grey80"),
        legend.position = "right") +
  scale_fill_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1")) +
  scale_colour_manual('Domain:', values = c("olivedrab3", "orchid3", "dodgerblue3", "darkorange1"))
PC12_plot_MW



### Using hausdorff sitance
## estimating centroid standard error
out_PCA_replacement_HD_MW <- data.frame()
for(i in 1:500){
  calc<- PCAdata_sp_final_MW %>% group_by(superkingdom) %>%
    sample_n(100, replace = TRUE) %>%
    select(PC1, PC2, superkingdom) %>%
    ungroup()
  
  
  BA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  EA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Archaea") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  AA_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == c("Archaea")) %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  VV_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Viruses") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  EE_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Eukaryota") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  BB_comp <- hausdorff_dist(calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix(),
                            calc  %>%
                              filter(superkingdom == "Bacteria") %>%
                              select(PC1, PC2) %>%
                              as.matrix())
  
  HD_estimates_MW <- rbind(data.frame(HD = BA_comp, Comparison = "Bacteria-Archaea", Type = "Between", rep = i),
                           data.frame(HD = BV_comp, Comparison = "Bacteria-Viruses", Type = "Between", rep = i),
                           data.frame(HD = BE_comp, Comparison = "Bacteria-Eukaryota",Type = "Between", rep = i),
                           data.frame(HD = EA_comp, Comparison = "Eukaryota-Archaea", Type = "Between", rep = i),
                           data.frame(HD = EV_comp, Comparison = "Eukaryota-Viruses", Type = "Between", rep = i),
                           data.frame(HD = AV_comp, Comparison = "Archaea-Viruses", Type = "Between", rep = i),
                           data.frame(HD = AA_comp, Comparison = "Archaea-Archaea", Type = "Within", rep = i),
                           data.frame(HD = BB_comp, Comparison = "Bacteria-Bacteria", Type = "Within", rep = i),
                           data.frame(HD = EE_comp, Comparison = "Eukaryota-Eukaryota", Type = "Within", rep = i),
                           data.frame(HD = VV_comp, Comparison = "Viruses-Viruses", Type = "Within", rep = i))
  
  
  out_PCA_replacement_HD_MW <- rbind(out_PCA_replacement_HD_MW, HD_estimates_MW)
  
  
}


HD_Plot_MW <- ggplot(data = out_PCA_replacement_HD_MW %>% filter(Type == "Between"),
                     aes(x = HD, fill = Comparison), col = "black") + 
  geom_histogram(bins = 45) +
  xlab("Hausdorff distance") + 
  scale_fill_brewer(palette = "Dark2")+ 
  theme_linedraw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, color = "black")) + 
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", col = "red") + 
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid", col = "black")

HD_Plot_MW 

final_plot_MW <- wrap_plots(PC12_plot_MW, HD_Plot_MW) + plot_annotation(title = "Mol Weight standardisation")

final_plot_raw / final_plot_MW







