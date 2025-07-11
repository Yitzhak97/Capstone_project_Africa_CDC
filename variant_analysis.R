#get working directory
getwd() #"/Users/sequencingplatform/Documents/Linux_Basics_to_Mastery_training/Module-15"

#load libraries
library(dplyr)
library(openxlsx)
library(janitor)
library(stringi)
library(stringr)
library(ggplot2)
library(data.table)

#load hq_allvariants.tsv
data_allvariants <- fread("vcf/hq_allvariants.tsv")
data_allvariants <- na.omit(data_allvariants)

#Calculate and plot the total number of high quality variants per sample

high_qual_vars <- data_allvariants %>%
  group_by(Sample) %>%
  summarise(count = n(), .groups = "drop") 
plot1 <- high_qual_vars %>%
ggplot(aes(x = Sample, y = count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "High Quality Variants per Sample",
       x = "Sample",
       y = "Variant Count") +
  theme_classic()
plot1

#Merge REF and ALT column to a new column mutation
data_allvariants <- data_allvariants %>%
  mutate(mutation = paste0(REF,ALT))

#Determine the type of each variant (Transition:A<->C, OR G<->C, or Transversion).
data_allvariants <- data_allvariants %>%
  mutate(Type = case_when(mutation %in% c("AT","TA","GC","CG") ~ "Transition",
                          mutation %in% c("AG", "GA", "AC", "CA", "TC", "CT", "TG", "GT") ~ "Transversion"))

is_transition <- function(ref, alt) {
  (ref == "A" & alt == "G") |
    (ref == "G" & alt == "A") |
    (ref == "C" & alt == "T") |
    (ref == "T" & alt == "C")
}

data_allvariants <- data_allvariants %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%  # keep only SNPs
  mutate(Type = ifelse(is_transition(REF, ALT), "Transition", "Transversion"))


#Analyze and plot the distribution of variant types per sample both count and proportions
Var_dist <- data_allvariants %>%
  na.omit() %>%
  group_by(Sample, Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prop = count/sum(count)*100)
plot2 <- Var_dist %>%
  ggplot(aes(x = Sample, prop, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of Ts vs. Tv variants per sample.",
    x = "Sample",
    y = "Proportion"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size =14),
    axis.text = element_text(size =10, angle = 45, hjust = 1)
    
  )

plot2   

#Identify common variants across samples or unique variants

## Create a unique identifier for each variant
data_allvariants$variant_id <- paste(data_allvariants$CHROM,
                                     data_allvariants$POS,
                                     data_allvariants$REF,
                                     data_allvariants$ALT,
                                     sep = "_")

## Keep one row per sample per variant (remove duplicates)
variant_sample_table <- data_allvariants %>%
  distinct(Sample, variant_id)

##Count how many samples each variant appears in
variant_counts <- variant_sample_table %>%
  group_by(variant_id) %>%
  summarise(sample_count = n(), .groups = "drop")

## Identify unique variants (appear in only one sample)
unique_variants <- variant_counts %>%
  filter(sample_count == 1)

## Identify common variants (appear in all samples)
common_variants <- variant_counts %>%
  filter(sample_count == total_samples)

## Join back to see which sample each unique variant belongs to
unique_variants_with_samples <- unique_variants %>%
  inner_join(variant_sample_table, by = "variant_id")


