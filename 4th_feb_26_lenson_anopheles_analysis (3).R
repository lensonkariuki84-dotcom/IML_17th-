# =============================================================================
# MALARIA VECTOR CHARACTERIZATION IN KENYA (2016-2024)
# Complete Data Cleaning, EDA & Analysis Script
# =============================================================================

# Load required packages
library(tidyverse)
library(janitor)
library(lubridate)
library(scales)
library(viridis)
library(patchwork)

# -----------------------------------------------------------------------------
# 1.1 LOAD AND INITIAL CLEANING
# -----------------------------------------------------------------------------

# Load data
anopheledata_nov25 <- read_csv("cleaned_mosquito_data.csv", show_col_types = FALSE)

# Check original structure
glimpse(anopheledata_nov25)
names(anopheledata_nov25)

# Clean column names and remove empty column
anopheles_clean <- anopheledata_nov25 %>%
  clean_names() %>%
  dplyr::select(-x17) %>%
  dplyr::rename(
    specimen_id = specimenid,
    year = yearofcollection,
    abdomen_status = abdomenstatus,
    house_no = houseno,
    collection_method = methodcollection,
    collection_place = place_of_collection,
    sporozoite_status = sporozoitestatus,
    bloodmeal_status = bloodmeasstatus
  )

# Check for duplicates
cat("Duplicate specimen IDs:", sum(duplicated(anopheles_clean$specimen_id)), "\n")

# -----------------------------------------------------------------------------
# 1.2 STANDARDIZE TEXT - TITLE CASE
# -----------------------------------------------------------------------------

# Function to clean location names (handles special characters)
clean_location_names <- function(x) {
  x %>%
    str_trim() %>%
    str_squish() %>%
    # Handle ng'a patterns
    str_replace_all("(?i)ng\\?'a", "nga") %>%
    str_replace_all("(?i)ng'a", "nga") %>%
    str_replace_all("(?i)ng`a", "nga") %>%
    # Remove remaining special characters
    str_replace_all("[''`'\"\\?]", "") %>%
    str_replace_all("[^[:alnum:][:space:]-]", "") %>%
    str_squish() %>%
    str_to_title()
}

# Apply to location columns
anopheles_clean <- anopheles_clean %>%
  mutate(
    county = clean_location_names(county),
    subcounty = clean_location_names(subcounty),
    ward = clean_location_names(ward),
    village = clean_location_names(village)
  )

# Fix specific county issues
anopheles_clean <- anopheles_clean %>%
  mutate(
    county = str_replace(county, "Murangaa", "Muranga"),
    county = str_replace(county, "Homabay", "Homa Bay")
  )

# Standardize other text columns
anopheles_clean <- anopheles_clean %>%
  mutate(
    across(c(species, subspecies, sex, abdomen_status, 
             collection_method, collection_place, 
             sporozoite_status, bloodmeal_status),
           ~ str_trim(.) %>% str_squish() %>% str_to_title())
  )

# Check counties after cleaning
print(sort(unique(anopheles_clean$county)))

# -----------------------------------------------------------------------------
# 1.3 STANDARDIZE MONTHS AND CREATE SEASONS
# -----------------------------------------------------------------------------

# Check current month values
print(unique(anopheles_clean$month))

# Standardize month names
month_mapping <- c(
  "January" = 1, "February" = 2, "March" = 3, "April" = 4,
  "May" = 5, "June" = 6, "July" = 7, "August" = 8,
  "September" = 9, "October" = 10, "November" = 11, "December" = 12,
  "Jan" = 1, "Feb" = 2, "Mar" = 3, "Apr" = 4,
  "Jun" = 6, "Jul" = 7, "Aug" = 8, "Sep" = 9, "Sept" = 9,
  "Oct" = 10, "Nov" = 11, "Dec" = 12
)

# Kenya's seasons
assign_season <- function(month_num) {
  case_when(
    month_num %in% c(1, 2) ~ "Short Dry (Jan-Feb)",
    month_num %in% c(3, 4, 5) ~ "Long Rains (Mar-May)",
    month_num %in% c(6, 7, 8, 9) ~ "Long Dry (Jun-Sep)",
    month_num %in% c(10, 11, 12) ~ "Short Rains (Oct-Dec)",
    TRUE ~ NA_character_
  )
}

anopheles_clean <- anopheles_clean %>%
  mutate(
    month_clean = str_to_title(str_trim(month)),
    month_num = month_mapping[month_clean],
    month_abbr = month.abb[month_num],
    season = assign_season(month_num),
    season = factor(season, levels = c(
      "Short Dry (Jan-Feb)", 
      "Long Rains (Mar-May)", 
      "Long Dry (Jun-Sep)", 
      "Short Rains (Oct-Dec)"
    ))
  )

# Check
cat("\nMonths after cleaning:\n")
anopheles_clean %>% 
  dplyr::count(month_clean, month_num, month_abbr) %>%
  print(n = 20)

# -----------------------------------------------------------------------------
# 1.4 ASSIGN EPIDEMIOLOGICAL ZONES
# -----------------------------------------------------------------------------

assign_malaria_burden <- function(county) {
  county <- str_to_title(str_trim(county))
  
  # No transmission
  no_transmission <- c("Nairobi")
  
  # Very low burden
  very_low <- c("Embu", "Kirinyaga", "Kitui", "Laikipia", "Lamu",
                "Machakos", "Makueni", "Mandera", "Muranga", 
                "Nyandarua", "Nyeri", "Wajir")
  
  # Low burden
  low <- c("Bomet", "Garissa", "Isiolo", "Kajiado", "Kiambu",
           "Marsabit", "Meru", "Nakuru", "Taita Taveta", "Tharaka Nithi")
  
  # Moderate burden
  moderate <- c("Elgeyo Marakwet", "Mombasa", "Narok", "Nyamira",
                "Samburu", "Tana River", "Uasin Gishu")
  
  # Moderate to high burden
  moderate_high <- c("Baringo", "Kericho", "Kilifi", "Kisii", "Nandi", 
                     "Trans Nzoia")
  
  # High burden
  high <- c("Bungoma", "Homa Bay", "Kwale", "Turkana", "West Pokot")
  
  # Very high burden
  very_high <- c("Busia", "Kakamega", "Kisumu", "Migori", "Siaya", "Vihiga")
  
  case_when(
    county %in% no_transmission ~ "No Transmission",
    county %in% very_low ~ "Very Low Burden",
    county %in% low ~ "Low Burden",
    county %in% moderate ~ "Moderate Burden",
    county %in% moderate_high ~ "Moderate-High Burden",
    county %in% high ~ "High Burden",
    county %in% very_high ~ "Very High Burden",
    TRUE ~ "Unclassified"
  )
}

# Apply to data
anopheles_clean <- anopheles_clean %>%
  mutate(
    malaria_burden = assign_malaria_burden(county),
    malaria_burden = factor(malaria_burden, levels = c(
      "No Transmission",
      "Very Low Burden",
      "Low Burden", 
      "Moderate Burden",
      "Moderate-High Burden",
      "High Burden",
      "Very High Burden",
      "Unclassified"
    ))
  )

# Check assignments
anopheles_clean %>%
  dplyr::count(malaria_burden, county) %>%
  print(n = 50)

# Check for unclassified counties
anopheles_clean %>%
  filter(malaria_burden == "Unclassified") %>%
  distinct(county)

# -----------------------------------------------------------------------------
# 1.5 CLEAN SPECIES AND SUBSPECIES
# -----------------------------------------------------------------------------

# Check current values
print(unique(anopheles_clean$species))
print(unique(anopheles_clean$subspecies))

# Clean species
anopheles_clean <- anopheles_clean %>%
  mutate(
    species_clean = case_when(
      str_detect(species, regex("gambiae", ignore_case = TRUE)) ~ "An. gambiae",
      str_detect(species, regex("funestus", ignore_case = TRUE)) ~ "An. funestus",
      str_detect(species, regex("unknown", ignore_case = TRUE)) ~ "Unknown",
      is.na(species) ~ NA_character_,
      TRUE ~ species
    )
  )

# Clean subspecies
anopheles_clean <- anopheles_clean %>%
  mutate(
    subspecies_clean = case_when(
      str_detect(subspecies, regex("arabiensis", ignore_case = TRUE)) ~ "An. arabiensis",
      str_detect(subspecies, regex("gambiae\\s*s\\.?s", ignore_case = TRUE)) ~ "An. gambiae s.s.",
      str_detect(subspecies, regex("funestus\\s*s\\.?s", ignore_case = TRUE)) ~ "An. funestus s.s.",
      str_detect(subspecies, regex("coluzzi", ignore_case = TRUE)) ~ "An. coluzzii",
      str_detect(subspecies, regex("leesoni", ignore_case = TRUE)) ~ "An. leesoni",
      str_detect(subspecies, regex("rivulorum", ignore_case = TRUE)) ~ "An. rivulorum",
      str_detect(subspecies, regex("parensis", ignore_case = TRUE)) ~ "An. parensis",
      str_detect(subspecies, regex("merus", ignore_case = TRUE)) ~ "An. merus",
      str_detect(subspecies, regex("vaneedeni", ignore_case = TRUE)) ~ "An. vaneedeni",
      str_detect(subspecies, regex("stephensi", ignore_case = TRUE)) ~ "An. stephensi",
      str_detect(subspecies, regex("not\\s*amplif", ignore_case = TRUE)) ~ "Not amplified",
      is.na(subspecies) ~ NA_character_,
      TRUE ~ subspecies
    )
  )

# Create combined species/subspecies column for detailed analysis
anopheles_clean <- anopheles_clean %>%
  mutate(
    species_detail = case_when(
      !is.na(subspecies_clean) & subspecies_clean != "Not amplified" ~ subspecies_clean,
      TRUE ~ species_clean
    )
  )

# Check cleaned species
print(sort(unique(anopheles_clean$species_clean)))

# Check cleaned subspecies
print(sort(unique(anopheles_clean$subspecies_clean)))

cat("\nSpecies detail:\n")
anopheles_clean %>%
  dplyr::count(species_clean, subspecies_clean, species_detail, sort = TRUE) %>%
  print(n = 20)

# -----------------------------------------------------------------------------
# 1.6 CLEAN SPOROZOITE STATUS
# -----------------------------------------------------------------------------

print(unique(anopheles_clean$sporozoite_status))

anopheles_clean <- anopheles_clean %>%
  mutate(
    sporozoite_clean = case_when(
      str_detect(sporozoite_status, regex("^pos", ignore_case = TRUE)) ~ "Positive",
      str_detect(sporozoite_status, regex("^neg", ignore_case = TRUE)) ~ "Negative",
      str_detect(sporozoite_status, regex("not\\s*done|^np$", ignore_case = TRUE)) ~ "Not Done",
      str_detect(sporozoite_status, regex("not\\s*applicable|^---$", ignore_case = TRUE)) ~ NA_character_,
      is.na(sporozoite_status) ~ NA_character_,
      TRUE ~ NA_character_
    ),
    sporozoite_positive = case_when(
      sporozoite_clean == "Positive" ~ 1L,
      sporozoite_clean == "Negative" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Check the result
cat("Sporozoite status (cleaned):\n")
anopheles_clean %>%
  dplyr::count(sporozoite_status, sporozoite_clean, sporozoite_positive) %>%
  print(n = 15)

# -----------------------------------------------------------------------------
# 1.7 CLEAN BLOODMEAL STATUS - CREATE CATEGORIES
# -----------------------------------------------------------------------------

print(unique(anopheles_clean$bloodmeal_status))

anopheles_clean <- anopheles_clean %>%
  mutate(
    # Check if human is present in bloodmeal
    has_human = str_detect(bloodmeal_status, regex("human", ignore_case = TRUE)),
    
    # Count number of hosts (comma separated)
    host_count = str_count(bloodmeal_status, ",") + 1,
    
    # Create bloodmeal categories
    bloodmeal_category = case_when(
      # Keep these as they are
      str_detect(bloodmeal_status, regex("^not\\s*applicable$", ignore_case = TRUE)) ~ "Not Applicable",
      str_detect(bloodmeal_status, regex("^not\\s*done$", ignore_case = TRUE)) ~ "Not Done",
      str_detect(bloodmeal_status, regex("^none$", ignore_case = TRUE)) ~ "None",
      is.na(bloodmeal_status) ~ NA_character_,
      
      # Three or more hosts = Mixed
      host_count >= 3 ~ "Mixed",
      
      # Two hosts - keep the combination as is
      host_count == 2 ~ str_to_title(bloodmeal_status),
      
      # Single host - keep as is
      host_count == 1 ~ str_to_title(bloodmeal_status),
      
      TRUE ~ NA_character_
    ),
    
    # Create human bloodmeal indicator (for HBI calculation)
    human_bloodmeal = case_when(
      bloodmeal_category %in% c("Not Applicable", "Not Done", "None") ~ NA_integer_,
      is.na(bloodmeal_category) ~ NA_integer_,
      has_human ~ 1L,
      !has_human ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Check the result
cat("Bloodmeal categories:\n")
anopheles_clean %>%
  dplyr::count(bloodmeal_category, sort = TRUE) %>%
  print(n = 35)

cat("\nDetailed bloodmeal mapping:\n")
anopheles_clean %>%
  dplyr::count(bloodmeal_status, bloodmeal_category, has_human, human_bloodmeal) %>%
  print(n = 35)

# -----------------------------------------------------------------------------
# 1.8 CLEAN ABDOMEN STATUS
# -----------------------------------------------------------------------------
print(unique(anopheles_clean$abdomen_status))

anopheles_clean <- anopheles_clean %>%
  mutate(
    abdomen_clean = case_when(
      str_detect(abdomen_status, regex("^unfed", ignore_case = TRUE)) ~ "Unfed",
      str_detect(abdomen_status, regex("^fed", ignore_case = TRUE)) ~ "Fed",
      str_detect(abdomen_status, regex("gravid", ignore_case = TRUE)) ~ "Gravid",
      str_detect(abdomen_status, regex("half", ignore_case = TRUE)) ~ "Half-gravid",
      str_detect(abdomen_status, regex("not\\s*(applicable|known)", ignore_case = TRUE)) ~ NA_character_,
      is.na(abdomen_status) ~ NA_character_,
      TRUE ~ abdomen_status
    )
  )

cat("\nAbdomen status (cleaned):\n")
anopheles_clean %>%
  dplyr::count(abdomen_status, abdomen_clean)

# -----------------------------------------------------------------------------
# 1.9 CREATE FINAL ANALYSIS DATASET
# -----------------------------------------------------------------------------

anopheles_final <- anopheles_clean %>%
  dplyr::select(
    # Identifiers
    specimen_id,
    
    # Temporal
    year,
    month_num,
    month_abbr,
    season,
    
    # Spatial
    county,
    subcounty,
    ward,
    village,
    malaria_burden,
    
    # Collection info
    house_no,
    collection_method,
    collection_place,
    
    # Taxonomy - BOTH levels
    species = species_clean,
    subspecies = subspecies_clean,
    species_detail,
    
    # Mosquito characteristics
    sex,
    abdomen_status = abdomen_clean,
    
    # Sporozoite outcomes
    sporozoite_status = sporozoite_clean,
    sporozoite_positive,
    
    # Bloodmeal outcomes
    bloodmeal_raw = bloodmeal_status,
    bloodmeal_category,
    # primary_host,
    has_human,
    human_bloodmeal
  )

# Filter to females only (for vector analysis)
anopheles_final <- anopheles_final %>%
  filter(sex == "Female" | is.na(sex))

# -----------------------------------------------------------------------------
# 1.10 DATA QUALITY SUMMARY
# -----------------------------------------------------------------------------
cat("Total records (females):", nrow(anopheles_final), "\n")
cat("Year range:", min(anopheles_final$year, na.rm = TRUE), "-", 
    max(anopheles_final$year, na.rm = TRUE), "\n")
cat("Unique counties:", n_distinct(anopheles_final$county, na.rm = TRUE), "\n")
cat("Unique villages:", n_distinct(anopheles_final$village, na.rm = TRUE), "\n\n")

# Missing data summary
missing_summary <- anopheles_final %>%
  dplyr::summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing") %>%
  mutate(
    Total = nrow(anopheles_final),
    Percent_Missing = round(Missing / Total * 100, 1)
  ) %>%
  arrange(desc(Percent_Missing))

cat("Missing Data Summary:\n")
print(missing_summary, n = 25)

# Final structure
cat("\nFinal dataset structure:\n")
glimpse(anopheles_final)

# Save cleaned data
# write_csv(anopheles_final, "anopheles_final_cleaned.csv")
cat("\nCleaned data saved to: anopheles_final_cleaned.csv\n")


# =============================================================================
# EXPLORATORY DATA ANALYSIS
# =============================================================================

# Set theme
theme_set(theme_minimal(base_size = 12) +
            theme(
              plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom",
              panel.grid.minor = element_blank()
            ))

# Color palettes
species_colors <- c(
  "An. gambiae" = "#E41A1C",
  "An. funestus" = "#377EB8",
  "Unknown" = "#999999"
)

subspecies_colors <- c(
  "An. arabiensis" = "#E41A1C",
  "An. gambiae s.s." = "#FF7F00",
  "An. funestus s.s." = "#377EB8",
  "An. leesoni" = "#4DAF4A",
  "An. rivulorum" = "#984EA3",
  "An. merus" = "#FFFF33",
  "An. parensis" = "#A65628",
  "An. coluzzii" = "#F781BF",
  "An. stephensi" = "#66C2A5",
  "An. vaneedeni" = "#8DA0CB",
  "Not amplified" = "#CCCCCC"
)

zone_colors <- c(
  "No Transmission" = "#1a9850",
  "Very Low Burden" = "#91cf60",
  "Low Burden" = "#d9ef8b",
  "Moderate Burden" = "#fee08b",
  "Moderate-High Burden" = "#fc8d59",
  "High Burden" = "#d73027",
  "Very High Burden" = "#a50026",
  "Unclassified" = "#999999"
)

# -----------------------------------------------------------------------------
# 2.1 OVERALL SUMMARIES
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("           EXPLORATORY DATA ANALYSIS\n")
cat("============================================================\n\n")

# Basic counts
cat("Total specimens:", nrow(anopheles_final), "\n")
cat("Years:", paste(sort(unique(anopheles_final$year)), collapse = ", "), "\n")
cat("Counties:", n_distinct(anopheles_final$county), "\n")
cat("Epidemiological zones:", n_distinct(anopheles_final$malaria_burden), "\n\n")

# Species summary
cat("=== SPECIES COMPOSITION ===\n")
anopheles_final %>%
  dplyr::count(species, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\n=== SUBSPECIES COMPOSITION ===\n")
anopheles_final %>%
  dplyr::count(subspecies, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\n=== DETAILED SPECIES COMPOSITION ===\n")
anopheles_final %>%
  dplyr::count(species, subspecies, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  print(n = 20)

# -----------------------------------------------------------------------------
# 2.2 VISUALIZATION - OVERVIEW
# -----------------------------------------------------------------------------

# Plot 1: Species composition
p1 <- anopheles_final %>%
  filter(!is.na(species)) %>%
  dplyr::count(species) %>%
  mutate(species = fct_reorder(species, n)) %>%
  ggplot(aes(x = n, y = species, fill = species)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(format(n, big.mark = ","), " (", round(n/sum(n)*100, 1), "%)")),
            hjust = -0.05, size = 4) +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title = "Species Composition",
    subtitle = paste("N =", format(nrow(anopheles_final), big.mark = ",")),
    x = "Number of Specimens",
    y = NULL
  )
p1

# Plot 2: Subspecies composition
p2 <- anopheles_final %>%
  filter(!is.na(subspecies)) %>%
  dplyr::count(subspecies) %>%
  mutate(subspecies = fct_reorder(subspecies, n)) %>%
  ggplot(aes(x = n, y = subspecies, fill = subspecies)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(format(n, big.mark = ","), " (", round(n/sum(n)*100, 1), "%)")),
            hjust = -0.05, size = 3.5) +
  scale_fill_manual(values = subspecies_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(
    title = "Subspecies Composition",
    x = "Number of Specimens",
    y = NULL
  )
p2

# Combine
p1 + p2 + plot_annotation(title = "Anopheles Species and Subspecies Composition")

# Plot 3: Annual distribution
p3 <- anopheles_final %>%
  dplyr::count(year) %>%
  ggplot(aes(x = factor(year), y = n)) +
  geom_col(fill = "#2166ac", alpha = 0.8) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Annual Sample Distribution", x = "Year", y = "Number of Specimens")
p3

# Plot 4: Monthly distribution
p4 <- anopheles_final %>%
  filter(!is.na(month_num)) %>%
  dplyr::count(month_num) %>%
  mutate(month_name = month.abb[month_num]) %>%
  ggplot(aes(x = reorder(month_name, month_num), y = n)) +
  geom_col(fill = "#2166ac", alpha = 0.8) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Monthly Sample Distribution", x = "Month", y = "Number of Specimens")
p4

# Alternative: Use scale_x_discrete with breaks
p_monthly_all <- anopheles_final %>%
  filter(!is.na(month_num)) %>%
  dplyr::count(month_num) %>%
  complete(month_num = 1:12, fill = list(n = 0)) %>%  # Fill missing months with 0
  mutate(month_name = factor(month.abb[month_num], levels = month.abb)) %>%
  ggplot(aes(x = month_name, y = n)) +
  geom_col(fill = "#2166ac", alpha = 0.8) +
  geom_text(aes(label = ifelse(n > 0, format(n, big.mark = ","), "")), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Monthly Sample Distribution",
    x = "Month",
    y = "Number of Specimens"
  )
p_monthly_all



# Plot 5: Seasonal distribution
p5 <- anopheles_final %>%
  filter(!is.na(season)) %>%
  dplyr::count(season) %>%
  ggplot(aes(x = season, y = n, fill = season)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.5, size = 3.5) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Seasonal Distribution", x = NULL, y = "Number of Specimens") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p5


# Plot 6: Epidemiological zones
p6 <- anopheles_final %>%
  filter(!is.na(malaria_burden)) %>%
  dplyr::count(malaria_burden) %>%
  mutate(malaria_burden = fct_reorder(malaria_burden, n)) %>%
  ggplot(aes(x = n, y = malaria_burden, fill = malaria_burden)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = zone_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Samples by Malaria Burden Zone", x = "Number of Specimens", y = NULL)
p6

# Combine temporal and spatial plots
(p3 + p4) / (p5 + p6) + 
  plot_annotation(title = "Temporal and Spatial Distribution of Samples")

# -----------------------------------------------------------------------------
# 2.3 SUBSPECIES BY SPECIES COMPLEX
# -----------------------------------------------------------------------------

# Subspecies within each species complex
p_subspecies_detail <- anopheles_final %>%
  filter(!is.na(subspecies), subspecies != "Not amplified") %>%
  dplyr::count(species, subspecies) %>%
  ggplot(aes(x = n, y = reorder(subspecies, n), fill = species)) +
  geom_col() +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3) +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~species, scales = "free_y", ncol = 1) +
  labs(
    title = "Subspecies Composition within Species Complexes",
    x = "Number of Specimens",
    y = NULL,
    fill = "Complex"
  ) +
  theme(legend.position = "none")

print(p_subspecies_detail)

# -----------------------------------------------------------------------------
# 2.4 SPECIES BY ZONE AND SEASON
# -----------------------------------------------------------------------------

# Species by zone
species_zone <- anopheles_final %>%
  filter(!is.na(species), !is.na(malaria_burden)) %>%
  dplyr::count(malaria_burden, species) %>%
  group_by(malaria_burden) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_species_zone <- species_zone %>%
  ggplot(aes(x = malaria_burden, y = prop, fill = species)) +
  geom_col(position = "fill", color = "white") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = species_colors) +
  labs(
    title = "Species Composition by Malaria Burden Zone",
    x = NULL, y = "Proportion", fill = "Species"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_species_zone

# Subspecies by zone
subspecies_zone <- anopheles_final %>%
  filter(!is.na(subspecies), subspecies != "Not amplified", !is.na(malaria_burden)) %>%
  dplyr::count(malaria_burden, subspecies) %>%
  group_by(malaria_burden) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_subspecies_zone <- subspecies_zone %>%
  ggplot(aes(x = malaria_burden, y = prop, fill = subspecies)) +
  geom_col(position = "fill", color = "white") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = subspecies_colors) +
  labs(
    title = "Subspecies Composition by Malaria Burden Zone",
    x = NULL, y = "Proportion", fill = "Subspecies"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_subspecies_zone

p_species_zone / p_subspecies_zone

# -----------------------------------------------------------------------------
# 2.5 COUNTY-LEVEL SUMMARY
# -----------------------------------------------------------------------------

county_summary <- anopheles_final %>%
  dplyr::count(county, malaria_burden, sort = TRUE)

p_county <- county_summary %>%
  head(38) %>%
  mutate(county = fct_reorder(county, n)) %>%
  ggplot(aes(x = n, y = county, fill = malaria_burden)) +
  geom_col() +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3) +
  scale_fill_manual(values = zone_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Top 20 Counties by Sample Size",
    x = "Number of Specimens",
    y = NULL,
    fill = "Epi Zone"
  )

print(p_county)


# Clean collection methods
anopheles_final <- anopheles_final %>%
  mutate(
    collection_method_clean = case_when(
      str_detect(collection_method, regex("light\\s*trap", ignore_case = TRUE)) ~ "CDC Light Trap",
      str_detect(collection_method, regex("pyrethrum", ignore_case = TRUE)) ~ "Pyrethrum Spray Catch",
      str_detect(collection_method, regex("mouth\\s*aspiration|^aspiration$", ignore_case = TRUE)) ~ "Mouth Aspiration",
      str_detect(collection_method, regex("prokopack", ignore_case = TRUE)) ~ "Prokopack Aspiration",
      str_detect(collection_method, regex("not\\s*known|^0$", ignore_case = TRUE)) ~ NA_character_,
      is.na(collection_method) ~ NA_character_,
      TRUE ~ NA_character_
    )
  )

# Check reclassification
anopheles_final %>%
  dplyr::count(collection_method, collection_method_clean) %>%
  print(n = 20)

# Summary
collection_summary <- anopheles_final %>%
  filter(!is.na(collection_method_clean)) %>%
  dplyr::count(collection_method_clean, collection_place) %>%
  arrange(desc(n))

print(collection_summary)

# Plot
p_collection <- anopheles_final %>%
  filter(!is.na(collection_method_clean)) %>%
  dplyr::count(collection_method_clean) %>%
  mutate(
    percent = round(n / sum(n) * 100, 1),
    collection_method_clean = fct_reorder(collection_method_clean, n)
  ) %>%
  ggplot(aes(x = n, y = collection_method_clean)) +
  geom_col(fill = "#4575b4") +
  geom_text(aes(label = paste0(format(n, big.mark = ","), " (", percent, "%)")), 
            hjust = -0.1, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Collection Methods",
    subtitle = paste("N =", format(sum(!is.na(anopheles_final$collection_method_clean)), big.mark = ",")),
    x = "Number of Specimens",
    y = NULL
  )

print(p_collection)

# =============================================================================
# OBJECTIVE 1: SPECIES COMPOSITION & SPATIO-TEMPORAL PATTERNS
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 SPECIES-LEVEL ANALYSIS
# -----------------------------------------------------------------------------

# Overall species composition
species_comp <- anopheles_final %>%
  filter(!is.na(species)) %>%
  dplyr::count(species, name = "n") %>%
  mutate(
    percent = round(n / sum(n) * 100, 1),
    cum_percent = cumsum(percent)
  )

print(species_comp)

# Species by epidemiological zone
species_by_zone <- anopheles_final %>%
  filter(!is.na(species), !is.na(malaria_burden)) %>%
  dplyr::count(malaria_burden, species) %>%
  group_by(malaria_burden) %>%
  mutate(
    zone_total = sum(n),
    percent = round(n / zone_total * 100, 1)
  ) %>%
  ungroup()

species_by_zone %>%
  dplyr::select(malaria_burden, species, n, percent) %>%
  pivot_wider(names_from = species, values_from = c(n, percent), values_fill = 0) %>%
  print()

# Species by season
species_by_season <- anopheles_final %>%
  filter(!is.na(species), !is.na(season)) %>%
  dplyr::count(season, species) %>%
  group_by(season) %>%
  mutate(
    season_total = sum(n),
    percent = round(n / season_total * 100, 1)
  ) %>%
  ungroup()

print(species_by_season)

# Species by year
species_by_year <- anopheles_final %>%
  filter(!is.na(species)) %>%
  dplyr::count(year, species) %>%
  group_by(year) %>%
  mutate(
    year_total = sum(n),
    percent = round(n / year_total * 100, 1)
  ) %>%
  ungroup()

print(species_by_year)

# -----------------------------------------------------------------------------
# 1.2 SUBSPECIES-LEVEL ANALYSIS
# -----------------------------------------------------------------------------

# Overall subspecies composition
subspecies_comp <- anopheles_final %>%
  filter(!is.na(subspecies)) %>%
  dplyr::count(subspecies, name = "n") %>%
  mutate(
    percent = round(n / sum(n) * 100, 1),
    cum_percent = cumsum(percent)
  ) %>%
  arrange(desc(n))

print(subspecies_comp)

# Subspecies by epidemiological zone
subspecies_by_zone <- anopheles_final %>%
  filter(!is.na(subspecies), !is.na(malaria_burden)) %>%
  dplyr::count(malaria_burden, subspecies) %>%
  group_by(malaria_burden) %>%
  mutate(
    zone_total = sum(n),
    percent = round(n / zone_total * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(malaria_burden, desc(n))

print(subspecies_by_zone, n = 40)

# Subspecies by season
cat("\n=== SUBSPECIES BY SEASON ===\n")
subspecies_by_season <- anopheles_final %>%
  filter(!is.na(subspecies), !is.na(season)) %>%
  dplyr::count(season, subspecies) %>%
  group_by(season) %>%
  mutate(
    season_total = sum(n),
    percent = round(n / season_total * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(season, desc(n))

print(subspecies_by_season, n = 40)

# -----------------------------------------------------------------------------
# 1.3 VISUALIZATIONS - SPECIES LEVEL
# -----------------------------------------------------------------------------

# Monthly patterns by species
p_monthly_species <- anopheles_final %>%
  filter(!is.na(month_num), !is.na(species)) %>%
  dplyr::count(month_num, species) %>%
  mutate(month_name = month.abb[month_num]) %>%
  ggplot(aes(x = reorder(month_name, month_num), y = n, fill = species)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = species_colors) +
  labs(
    title = "Monthly Species Composition",
    x = "Month", y = "Number of Specimens", fill = "Species"
  )
p_monthly_species

# stacked barplot
p_species_zone <- species_by_zone %>%
  ggplot(aes(x = malaria_burden, y = percent, fill = species)) +
  geom_col(position = "fill", color = "white", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  # scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = species_colors) +
  labs(
    title = "Species Composition by Malaria Burden Zone",
    x = NULL, y = "Proportion", fill = "Species"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_species_zone


# Annual trends
p_annual_species <- species_by_year %>%
  ggplot(aes(x = factor(year), y = n, fill = species)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = species_colors) +
  labs(
    title = "Annual Species Abundance",
    x = "Year", y = "Number of Specimens", fill = "Species"
  )

p_annual_species

(p_monthly_species + p_species_zone) / p_annual_species

# -----------------------------------------------------------------------------
# 1.4 VISUALIZATIONS - SUBSPECIES LEVEL
# -----------------------------------------------------------------------------

# Monthly patterns by subspecies
p_monthly_subspecies <- anopheles_final %>%
  filter(!is.na(month_num), !is.na(subspecies), subspecies != "Not amplified") %>%
  dplyr::count(month_num, subspecies) %>%
  mutate(month_name = month.abb[month_num]) %>%
  ggplot(aes(x = reorder(month_name, month_num), y = n, fill = subspecies)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = subspecies_colors) +
  labs(
    title = "Monthly Subspecies Composition",
    x = "Month", y = "Number of Specimens", fill = "Subspecies"
  ) +
  theme(legend.position = "right")
p_monthly_subspecies

# Stacked bar chart for subspecies
p_subspecies_zone <- subspecies_by_zone %>%
  ggplot(aes(x = malaria_burden, y = percent, fill = subspecies)) +
  geom_col(position = "fill", color = "white", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  # scale_fill_brewer(palette = "Paired") +
  scale_fill_manual(values = subspecies_colors) +
  labs(
    title = "Subspecies Composition by Malaria Burden Zone",
    x = NULL, y = "Proportion", fill = "Subspecies"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
p_subspecies_zone

p_monthly_subspecies / p_subspecies_zone


# Subspecies by Year - Stacked Bar Plot
subspecies_by_year <- anopheles_final %>%
  filter(!is.na(subspecies), subspecies != "Not Amplified", !is.na(year)) %>%
  dplyr::count(year, subspecies) %>%
  group_by(year) %>%
  mutate(
    year_total = sum(n),
    prop = n / year_total
  ) %>%
  ungroup()

# Check data
print(subspecies_by_year, n = 30)


# Stacked bar plot with manual colors
p_subspecies_year <- subspecies_by_year %>%
  ggplot(aes(x = factor(year), y = prop, fill = subspecies)) +
  geom_col(position = "fill", color = "white", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = subspecies_colors) +
  labs(
    title = "Subspecies Composition by Year of Collection",
    subtitle = paste("N =", format(sum(subspecies_by_year$n), big.mark = ",")),
    x = "Year",
    y = "Proportion",
    fill = "Subspecies"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

print(p_subspecies_year)


# Subspecies by County - Stacked Bar Plot
subspecies_by_county <- anopheles_final %>%
  filter(!is.na(subspecies), subspecies != "Not Amplified", !is.na(county)) %>%
  dplyr::count(county, subspecies) %>%
  group_by(county) %>%
  mutate(
    county_total = sum(n),
    prop = n / county_total
  ) %>%
  ungroup() %>%
  filter(county_total >= 20)  # Filter counties with sufficient samples

# Check data
print(subspecies_by_county, n = 30)

# Order counties by total samples
county_order <- subspecies_by_county %>%
  group_by(county) %>%
  dplyr::summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  pull(county)

# Horizontal stacked bar for counties with manual colors
p_subspecies_county_h <- subspecies_by_county %>%
  mutate(county = factor(county, levels = rev(county_order))) %>%
  ggplot(aes(x = prop, y = county, fill = subspecies)) +
  geom_col(position = "fill", color = "white", width = 0.8) +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = subspecies_colors) +
  labs(
    title = "Subspecies Composition by County",
    subtitle = "Counties with ≥20 samples",
    x = "Proportion",
    y = NULL,
    fill = "Subspecies"
  ) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_subspecies_county_h)

# -----------------------------------------------------------------------------
# 1.5 SEASONAL PATTERNS WITH CONFIDENCE
# -----------------------------------------------------------------------------

# Seasonal abundance by species
seasonal_species_summary <- anopheles_final %>%
  filter(!is.na(season), !is.na(species)) %>%
  group_by(season, species) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(season) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    se = sqrt(prop * (1 - prop) / total),
    ci_lower = pmax(0, prop - 1.96 * se),
    ci_upper = pmin(1, prop + 1.96 * se)
  )

p_seasonal_prop <- seasonal_species_summary %>%
  ggplot(aes(x = season, y = prop, fill = species)) +
  geom_col(position = "dodge", color = "white") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = species_colors) +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Seasonal Species Proportions with 95% CI",
    x = NULL, y = "Proportion", fill = "Species"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_seasonal_prop)

# -----------------------------------------------------------------------------
# 1.6 OBJECTIVE 1 SUMMARY TABLES
# -----------------------------------------------------------------------------

# Summary table: Species by Zone and Season
obj1_summary_species <- anopheles_final %>%
  filter(!is.na(malaria_burden), !is.na(season), !is.na(species)) %>%
  dplyr::count(malaria_burden, season, species, name = "abundance") %>%
  arrange(malaria_burden, season, desc(abundance))

cat("\n=== OBJECTIVE 1 SUMMARY: SPECIES BY ZONE & SEASON ===\n")
print(obj1_summary_species, n = 50)

# Summary table: Subspecies by Zone and Season  
obj1_summary_subspecies <- anopheles_final %>%
  filter(!is.na(malaria_burden), !is.na(season), !is.na(subspecies), 
         subspecies != "Not amplified") %>%
  dplyr::count(malaria_burden, season, subspecies, name = "abundance") %>%
  arrange(malaria_burden, season, desc(abundance))

cat("\n=== OBJECTIVE 1 SUMMARY: SUBSPECIES BY ZONE & SEASON ===\n")
print(obj1_summary_subspecies, n = 50)


# =============================================================================
# OBJECTIVE 2: SPOROZOITE INFECTION RATES
# =============================================================================
# Filter to specimens with sporozoite data
spz_data <- anopheles_final %>%
  filter(!is.na(sporozoite_status))

cat("Specimens with sporozoite data:", nrow(spz_data), "\n")
cat("Positive:", sum(spz_data$sporozoite_positive, na.rm = TRUE), "\n")
cat("Negative:", sum(spz_data$sporozoite_positive == 0, na.rm = TRUE), "\n\n")

# Function to calculate sporozoite rate with CI
calc_spz_rate <- function(data) {
  data %>%
    dplyr::summarise(
      n_tested = n(),
      n_positive = sum(sporozoite_positive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      spz_rate = n_positive / n_tested,
      se = sqrt(spz_rate * (1 - spz_rate) / n_tested),
      ci_lower = pmax(0, spz_rate - 1.96 * se),
      ci_upper = pmin(1, spz_rate + 1.96 * se),
      spz_rate_pct = round(spz_rate * 100, 2),
      ci_lower_pct = round(ci_lower * 100, 2),
      ci_upper_pct = round(ci_upper * 100, 2)
    )
}

# -----------------------------------------------------------------------------
# 2.1 OVERALL SPOROZOITE RATE
# -----------------------------------------------------------------------------

overall_spz <- spz_data %>% calc_spz_rate()

cat("Tested:", overall_spz$n_tested, "\n")
cat("Positive:", overall_spz$n_positive, "\n")
cat("Rate:", overall_spz$spz_rate_pct, "% (95% CI:", 
    overall_spz$ci_lower_pct, "-", overall_spz$ci_upper_pct, "%)\n\n")

# -----------------------------------------------------------------------------
# 2.2 SPOROZOITE RATE BY SPECIES
# -----------------------------------------------------------------------------

spz_by_species <- spz_data %>%
  filter(!is.na(species)) %>%
  group_by(species) %>%
  calc_spz_rate() %>%
  arrange(desc(spz_rate))

print(spz_by_species %>% 
        dplyr::select(species, n_tested, n_positive, spz_rate_pct, ci_lower_pct, ci_upper_pct))

# -----------------------------------------------------------------------------
# 2.3 SPOROZOITE RATE BY SUBSPECIES
# -----------------------------------------------------------------------------

spz_by_subspecies <- spz_data %>%
  filter(!is.na(subspecies)) %>%
  group_by(subspecies) %>%
  calc_spz_rate() %>%
  filter(n_tested >= 1) %>%
  arrange(desc(spz_rate))

print(spz_by_subspecies %>% 
        dplyr::select(subspecies, n_tested, n_positive, spz_rate_pct, ci_lower_pct, ci_upper_pct))

# Visualization
p_spz_subspecies <- spz_by_subspecies %>%
  filter(n_tested >= 1) %>%
  mutate(subspecies = fct_reorder(subspecies, spz_rate)) %>%
  ggplot(aes(x = spz_rate * 100, y = subspecies)) +
  geom_point(aes(size = n_tested), color = "#d73027") +
  geom_errorbarh(aes(xmin = ci_lower * 100, xmax = ci_upper * 100), height = 0.3) +
  geom_text(aes(label = paste0("n=", n_tested)), hjust = -0.3, size = 3) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    title = "Sporozoite Rates by Subspecies",
    subtitle = "Error bars = 95% CI",
    x = "Sporozoite Rate (%)", y = NULL
  ) +
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_spz_subspecies)


p_spz_subspecies <- spz_by_subspecies %>%
  filter(n_tested >= 1) %>%
  mutate(subspecies = fct_reorder(subspecies, spz_rate)) %>%
  ggplot(aes(x = spz_rate * 100, y = subspecies)) +
  geom_col(fill = "#d73027", alpha = 0.3, width = 0.7) +
  geom_point(aes(size = n_tested), color = "#d73027") +
  geom_errorbarh(aes(xmin = ci_lower * 100, xmax = ci_upper * 100), 
                 height = 0.2, color = "black", linewidth = 0.8) +
  geom_text(aes(label = paste0(spz_rate_pct, "% (n=", n_tested, ")")), 
            hjust = -0.1, size = 3) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title = "Sporozoite Rates by Subspecies",
    subtitle = "Error bars = 95% CI; Point size = sample size",
    x = "Sporozoite Rate (%)", 
    y = NULL
  ) +
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_spz_subspecies)

# -----------------------------------------------------------------------------
# 2.4 SPOROZOITE RATE BY EPIDEMIOLOGICAL ZONE
# -----------------------------------------------------------------------------

spz_by_zone <- spz_data %>%
  filter(!is.na(malaria_burden), malaria_burden != "Unclassified") %>%
  group_by(malaria_burden) %>%
  calc_spz_rate() %>%
  arrange(desc(spz_rate))

print(spz_by_zone %>% 
        dplyr::select(malaria_burden, n_tested, n_positive, spz_rate_pct, ci_lower_pct, ci_upper_pct))

p_spz_zone <- spz_by_zone %>%
  mutate(malaria_burden = fct_reorder(malaria_burden, spz_rate)) %>%
  ggplot(aes(x = spz_rate * 100, y = malaria_burden, fill = malaria_burden)) +
  geom_col(show.legend = FALSE) +
  geom_errorbarh(aes(xmin = ci_lower * 100, xmax = ci_upper * 100), height = 0.3) +
  geom_text(aes(label = paste0(spz_rate_pct, "% (n=", n_tested, ")")), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = zone_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(title = "Sporozoite Rates by Malaria Burden Zone", x = "Sporozoite Rate (%)", y = NULL)

print(p_spz_zone)

# -----------------------------------------------------------------------------
# 2.5 SPOROZOITE RATE BY SEASON
# -----------------------------------------------------------------------------

spz_by_season <- spz_data %>%
  filter(!is.na(season)) %>%
  group_by(season) %>%
  calc_spz_rate()

print(spz_by_season %>% 
        dplyr::select(season, n_tested, n_positive, spz_rate_pct, ci_lower_pct, ci_upper_pct))

# -----------------------------------------------------------------------------
# 2.6 SPOROZOITE RATE: SUBSPECIES × ZONE
# -----------------------------------------------------------------------------

spz_subspecies_zone <- spz_data %>%
  filter(!is.na(subspecies),
         !is.na(malaria_burden), malaria_burden != "Unclassified") %>%
  group_by(malaria_burden, subspecies) %>%
  calc_spz_rate() %>%
  filter(n_tested >= 5) %>%
  arrange(malaria_burden, desc(spz_rate))

print(spz_subspecies_zone %>% 
        dplyr::select(malaria_burden, subspecies, n_tested, n_positive, spz_rate_pct), n = 40)

# Stacked bar plot - Subspecies by Zone showing sporozoite rates
p_spz_zone_stacked <- spz_subspecies_zone %>%
  ggplot(aes(x = malaria_burden, y = spz_rate * 100, fill = subspecies)) +
  geom_col(position = "dodge", color = "white", width = 0.8) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Sporozoite Rates by Subspecies and Malaria Burden Zone",
    subtitle = "Error bars = 95% CI",
    x = NULL, 
    y = "Sporozoite Rate (%)", 
    fill = "Subspecies"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_spz_zone_stacked)


# # Stacked bar plot - Subspecies stacked within each zone
# p_spz_zone_stacked <- spz_subspecies_zone %>%
#   ggplot(aes(x = malaria_burden, y = spz_rate * 100, fill = subspecies)) +
#   geom_col(position = "stack", color = "white", width = 0.8) +
#   scale_fill_brewer(palette = "Paired") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     title = "Sporozoite Rates by Subspecies and Malaria Burden Zone",
#     x = NULL, 
#     y = "Sporozoite Rate (%)", 
#     fill = "Subspecies"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# print(p_spz_zone_stacked)

# Calculate sporozoite rates by year
spz_by_year <- anopheles_final %>%
  filter(!is.na(sporozoite_status), sporozoite_status %in% c("Positive", "Negative")) %>%
  group_by(year) %>%
  dplyr::summarise(
    n_tested = n(),
    n_positive = sum(sporozoite_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    spz_rate = n_positive / n_tested,
    se = sqrt(spz_rate * (1 - spz_rate) / n_tested),
    ci_lower = pmax(0, spz_rate - 1.96 * se),
    ci_upper = pmin(1, spz_rate + 1.96 * se),
    spz_rate_pct = round(spz_rate * 100, 2),
    ci_lower_pct = round(ci_lower * 100, 2),
    ci_upper_pct = round(ci_upper * 100, 2)
  )

# Check data
print(spz_by_year)


# Bar plot with error bars and trend line
p_spz_year <- spz_by_year %>%
  ggplot(aes(x = factor(year), y = spz_rate * 100)) +
  geom_col(fill = "#d73027", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100),
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_line(aes(group = 1), color = "#2166ac", linewidth = 1.2) +
  geom_point(color = "#2166ac", size = 3) +
  geom_text(aes(label = paste0(spz_rate_pct, "%\n(n=", n_tested, ")")),
            vjust = -0.5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Sporozoite Prevalence by Year",
    subtitle = "Error bars = 95% CI; Blue line = trend",
    x = "Year",
    y = "Sporozoite Rate (%)"
  ) +
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_spz_year)


# Bar plot with smoothed trend line
p_spz_year <- spz_by_year %>%
  mutate(year_num = as.numeric(as.character(year))) %>%
  ggplot(aes(x = factor(year), y = spz_rate * 100)) +
  geom_col(fill = "#d73027", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100),
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_smooth(aes(x = as.numeric(factor(year)), y = spz_rate * 100),
              method = "loess", se = TRUE, color = "#2166ac", 
              fill = "#2166ac", alpha = 0.2, linewidth = 1.2) +
  geom_text(aes(label = paste0(spz_rate_pct, "%\n(n=", n_tested, ")")),
            vjust = -0.5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Sporozoite Prevalence by Year",
    subtitle = "Error bars = 95% CI; Blue line = LOESS trend with 95% CI",
    x = "Year",
    y = "Sporozoite Rate (%)"
  ) +
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_spz_year)


# Sporozoite rates by year and subspecies
spz_year_subspecies <- anopheles_final %>%
  filter(!is.na(sporozoite_status), sporozoite_status %in% c("Positive", "Negative"),
         !is.na(subspecies), subspecies != "Not Amplified") %>%
  group_by(year, subspecies) %>%
  dplyr::summarise(
    n_tested = n(),
    n_positive = sum(sporozoite_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    spz_rate = n_positive / n_tested,
    se = sqrt(spz_rate * (1 - spz_rate) / n_tested),
    ci_lower = pmax(0, spz_rate - 1.96 * se),
    ci_upper = pmin(1, spz_rate + 1.96 * se),
    spz_rate_pct = round(spz_rate * 100, 2)
  ) %>%
  filter(n_tested >= 10)  # Filter for sufficient samples

# Grouped bar plot
p_spz_year_subspecies <- spz_year_subspecies %>%
  ggplot(aes(x = factor(year), y = spz_rate * 100, fill = subspecies)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower * 100, ymax = ci_upper * 100),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Sporozoite Prevalence by Year and Subspecies",
    subtitle = "Error bars = 95% CI",
    x = "Year",
    y = "Sporozoite Rate (%)",
    fill = "Subspecies"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_spz_year_subspecies)

# =============================================================================
# OBJECTIVE 3: BLOOD MEAL SOURCES & HOST PREFERENCE
# =============================================================================
# Filter to specimens with bloodmeal data
bloodmeal_data <- anopheles_final %>%
  filter(!is.na(bloodmeal_category))

cat("Specimens with bloodmeal data:", nrow(bloodmeal_data), "\n\n")

# -----------------------------------------------------------------------------
# 3.1 BLOODMEAL CATEGORY COMPOSITION
# -----------------------------------------------------------------------------
# Filter bloodmeal data for analysis
bloodmeal_data <- anopheles_final %>%
  filter(!is.na(bloodmeal_category),
         !bloodmeal_category %in% c("Not Applicable", "Not Done", "None", "NA"))

# Check what's left
bloodmeal_data %>%
  dplyr::count(bloodmeal_category, sort = TRUE)


cat("=== BLOODMEAL CATEGORY COMPOSITION ===\n")
bloodmeal_summary <- bloodmeal_data %>%
  dplyr::count(bloodmeal_category, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

print(bloodmeal_summary)

# Detailed host composition
cat("\n=== DETAILED HOST COMPOSITION ===\n")
host_detail <- bloodmeal_data %>%
  dplyr::count(bloodmeal_raw, bloodmeal_category, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

print(host_detail, n = 30)

# Visualization
host_colors <- c(
  "Human" = "#d73027",
  "Bovine" = "#4575b4",
  "Goat" = "#91bfdb",
  "Dog" = "#fc8d59",
  "Cat" = "#fee090",
  "Chicken" = "#f46d43",
  "None" = "#999999",
  "Goat, Bovine" = "#8c510a",
  "Goat, Human" = "#bf812d",
  "Human, Bovine" = "#dfc27d",
  "Cat, Human" = "#c7eae5",
  "Human, Dog" = "#80cdc1",
  "Cat, Goat" = "#35978f",
  "Cat, Dog" = "#01665e",
  "Mixed" = "#984ea3"
)



p_bloodmeal_cat <- bloodmeal_summary %>%
  mutate(bloodmeal_category = fct_reorder(bloodmeal_category, n)) %>%
  ggplot(aes(x = n, y = bloodmeal_category, fill = bloodmeal_category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", percent, "%)")), hjust = -0.1, size = 4) +
  scale_fill_manual(values = host_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title = "Blood Meal Categories",
    subtitle = paste("N =", nrow(bloodmeal_data)),
    x = "Number of Blood Meals", y = NULL
  )

print(p_bloodmeal_cat)



# Create bloodmeal data with reclassified categories
bloodmeal_data <- anopheles_final %>%
  filter(!is.na(bloodmeal_category),
         !bloodmeal_category %in% c("Not Applicable", "Not Done", "None", "NA")) %>%
  mutate(
    bloodmeal_simple = case_when(
      bloodmeal_category == "Human" ~ "Human",
      bloodmeal_category == "Bovine" ~ "Bovine",
      bloodmeal_category == "Goat" ~ "Goat",
      bloodmeal_category == "Dog" ~ "Dog",
      bloodmeal_category == "Cat" ~ "Cat",
      bloodmeal_category == "Chicken" ~ "Chicken",
      str_detect(bloodmeal_category, ",") ~ "Mixed",  # Any with comma = Mixed
      bloodmeal_category == "Mixed" ~ "Mixed",
      TRUE ~ bloodmeal_category
    )
  )

# Check reclassification
bloodmeal_data %>%
  dplyr::count(bloodmeal_simple, sort = TRUE)

# Create summary
bloodmeal_summary <- bloodmeal_data %>%
  dplyr::count(bloodmeal_simple, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

print(bloodmeal_summary)

# Updated color palette
host_colors <- c(
  "Human" = "#d73027",
  "Bovine" = "#4575b4",
  "Goat" = "#91bfdb",
  "Dog" = "#fc8d59",
  "Cat" = "#fee090",
  "Chicken" = "#f46d43",
  # "None" = "#999999",
  "Mixed" = "#984ea3"
)

# Plot
p_bloodmeal_cat_new <- bloodmeal_summary %>%
  mutate(bloodmeal_simple = fct_reorder(bloodmeal_simple, n)) %>%
  ggplot(aes(x = n, y = bloodmeal_simple, fill = bloodmeal_simple)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", percent, "%)")), hjust = -0.1, size = 4) +
  scale_fill_manual(values = host_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title = "Blood Meal Host Composition",
    subtitle = paste("N =", format(nrow(bloodmeal_data), big.mark = ",")),
    x = "Number of Blood Meals", 
    y = NULL
  )

print(p_bloodmeal_cat_new)


# -----------------------------------------------------------------------------
# 3.2 HUMAN BLOOD INDEX (HBI)
# -----------------------------------------------------------------------------

# Function to calculate HBI with CI
calc_hbi <- function(data) {
  data %>%
    dplyr::summarise(
      n_total = n(),
      n_human = sum(human_bloodmeal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      hbi = n_human / n_total,
      se = sqrt(hbi * (1 - hbi) / n_total),
      ci_lower = pmax(0, hbi - 1.96 * se),
      ci_upper = pmin(1, hbi + 1.96 * se),
      hbi_pct = round(hbi * 100, 1),
      ci_lower_pct = round(ci_lower * 100, 1),
      ci_upper_pct = round(ci_upper * 100, 1)
    )
}

# Overall HBI
overall_hbi <- bloodmeal_data %>% calc_hbi()

cat("\n=== OVERALL HUMAN BLOOD INDEX ===\n")
cat("Total blood meals:", overall_hbi$n_total, "\n")
cat("Human blood meals:", overall_hbi$n_human, "\n")
cat("HBI:", overall_hbi$hbi_pct, "% (95% CI:", 
    overall_hbi$ci_lower_pct, "-", overall_hbi$ci_upper_pct, "%)\n\n")

# -----------------------------------------------------------------------------
# 3.3 HBI BY SPECIES
# -----------------------------------------------------------------------------

hbi_by_species <- bloodmeal_data %>%
  filter(!is.na(species)) %>%
  group_by(species) %>%
  calc_hbi() %>%
  arrange(desc(hbi))

print(hbi_by_species %>% 
        dplyr::select(species, n_total, n_human, hbi_pct, ci_lower_pct, ci_upper_pct))

# -----------------------------------------------------------------------------
# 3.4 HBI BY SUBSPECIES
# -----------------------------------------------------------------------------

hbi_by_subspecies <- bloodmeal_data %>%
  filter(!is.na(subspecies)) %>%
  group_by(subspecies) %>%
  calc_hbi() %>%
  filter(n_total >= 1) %>%
  arrange(desc(hbi))

print(hbi_by_subspecies %>% 
        dplyr::select(subspecies, n_total, n_human, hbi_pct, ci_lower_pct, ci_upper_pct))

# Visualization
p_hbi_subspecies <- hbi_by_subspecies %>%
  filter(n_total >= 1) %>%
  mutate(subspecies = fct_reorder(subspecies, hbi)) %>%
  ggplot(aes(x = hbi * 100, y = subspecies)) +
  geom_point(aes(size = n_total), color = "#d73027") +
  geom_errorbarh(aes(xmin = ci_lower * 100, xmax = ci_upper * 100), height = 0.3) +
  geom_vline(xintercept = overall_hbi$hbi * 100, linetype = "dashed", color = "gray50") +
  scale_size_continuous(range = c(2, 8), name = "N blood meals") +
  labs(
    title = "Human Blood Index by Subspecies",
    subtitle = "Dashed line = overall HBI",
    x = "HBI (%)", y = NULL
  ) +
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_hbi_subspecies)

# -----------------------------------------------------------------------------
# 3.5 HBI BY EPIDEMIOLOGICAL ZONE
# -----------------------------------------------------------------------------

hbi_by_zone <- bloodmeal_data %>%
  filter(!is.na(malaria_burden), malaria_burden != "Unclassified") %>%
  group_by(malaria_burden) %>%
  calc_hbi() %>%
  arrange(desc(hbi))

print(hbi_by_zone %>% 
        dplyr::select(malaria_burden, n_total, n_human, hbi_pct, ci_lower_pct, ci_upper_pct))

p_hbi_zone <- hbi_by_zone %>%
  mutate(malaria_burden = fct_reorder(malaria_burden, hbi)) %>%
  ggplot(aes(x = hbi * 100, y = malaria_burden, fill = malaria_burden)) +
  geom_col(show.legend = FALSE) +
  geom_errorbarh(aes(xmin = ci_lower * 100, xmax = ci_upper * 100), height = 0.3) +
  geom_text(aes(label = paste0(hbi_pct, "% (n=", n_total, ")")), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = zone_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(title = "Human Blood Index by Epidemiological Zone", x = "HBI (%)", y = NULL)

print(p_hbi_zone)

# -----------------------------------------------------------------------------
# 3.6 HBI BY SEASON
# -----------------------------------------------------------------------------

hbi_by_season <- bloodmeal_data %>%
  filter(!is.na(season)) %>%
  group_by(season) %>%
  calc_hbi()

   print(hbi_by_season %>% 
        dplyr::select(season, n_total, n_human, hbi_pct, ci_lower_pct, ci_upper_pct))

# -----------------------------------------------------------------------------
# 3.7 HBI: SUBSPECIES × ZONE
# -----------------------------------------------------------------------------

hbi_subspecies_zone <- bloodmeal_data %>%
  filter(!is.na(subspecies),
         !is.na(malaria_burden), malaria_burden != "Unclassified") %>%
  group_by(malaria_burden, subspecies) %>%
  calc_hbi() %>%
  filter(n_total >= 5) %>%
  arrange(malaria_burden, desc(hbi))

print(hbi_subspecies_zone %>% 
        dplyr::select(malaria_burden, subspecies, n_total, n_human, hbi_pct), n = 40)


# Stacked bar plot - HBI by subspecies across zones
# p_hbi_stacked <- hbi_subspecies_zone %>%
#   ggplot(aes(x = malaria_burden, y = hbi * 100, fill = subspecies)) +
#   geom_col(position = "stack", color = "white", width = 0.8) +
#   scale_fill_brewer(palette = "Paired") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     title = "Human Blood Index by Subspecies and Malaria Burden Zone",
#     x = NULL, 
#     y = "HBI (%)", 
#     fill = "Subspecies"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# print(p_hbi_stacked)

# Grouped bar plot - HBI by subspecies across zones
p_hbi_grouped <- hbi_subspecies_zone %>%
  ggplot(aes(x = malaria_burden, y = hbi * 100, fill = subspecies)) +
  geom_col(position = "dodge", color = "white", width = 0.8) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Human Blood Index by Subspecies and Malaria Burden Zone",
    x = NULL, 
    y = "HBI (%)", 
    fill = "Subspecies"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_hbi_grouped)



# -----------------------------------------------------------------------------
# 3.8 HOST COMPOSITION BY SUBSPECIES
# -----------------------------------------------------------------------------
cat("\n=== HOST COMPOSITION BY SUBSPECIES ===\n")
# host_by_subspecies <- bloodmeal_data %>%
#   filter(!is.na(subspecies)) %>%
#   dplyr::count(subspecies, bloodmeal_category) %>%
#   group_by(subspecies) %>%
#   mutate(prop = n / sum(n), total = sum(n)) %>%
#   ungroup() %>%
#   filter(total >= 1)
# 
# print(host_by_subspecies %>% arrange(subspecies, desc(n)), n = 40)
# 
# # Visualization
# p_host_subspecies <- host_by_subspecies %>%
#   ggplot(aes(x = subspecies, y = prop, fill = bloodmeal_category)) +
#   geom_col(position = "fill", color = "white") +
#   scale_y_continuous(labels = percent) +
#   scale_fill_manual(values = host_colors) +
#   labs(
#     title = "Blood Meal Host Composition by Subspecies",
#     x = NULL, y = "Proportion", fill = "Category"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(p_host_subspecies)



# Reclassify bloodmeal categories in host_by_subspecies
host_by_subspecies <- anopheles_final %>%
  filter(!is.na(subspecies),
         !is.na(bloodmeal_category),
         !bloodmeal_category %in% c("Not Applicable", "Not Done", "None")) %>%
  mutate(
    bloodmeal_simple = case_when(
      bloodmeal_category == "Human" ~ "Human",
      bloodmeal_category == "Bovine" ~ "Bovine",
      bloodmeal_category == "Goat" ~ "Goat",
      bloodmeal_category == "Dog" ~ "Dog",
      bloodmeal_category == "Cat" ~ "Cat",
      bloodmeal_category == "Chicken" ~ "Chicken",
      str_detect(bloodmeal_category, ",") ~ "Mixed",
      bloodmeal_category == "Mixed" ~ "Mixed",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::count(subspecies, bloodmeal_simple) %>%
  group_by(subspecies) %>%
  mutate(
    total = sum(n),
    prop = n / total
  ) %>%
  ungroup() %>%
  filter(total >= 1)

# Check categories
host_by_subspecies %>%
  distinct(bloodmeal_simple)

# Updated color palette (only 8 colors needed)
host_colors <- c(
  "Human" = "#d73027",
  "Bovine" = "#4575b4",
  "Goat" = "#91bfdb",
  "Dog" = "#fc8d59",
  "Cat" = "#fee090",
  "Chicken" = "#f46d43",
  # "None" = "#999999",
  "Mixed" = "#984ea3"
)

# Plot
p_host_subspecies <- host_by_subspecies %>%
  ggplot(aes(x = subspecies, y = prop, fill = bloodmeal_simple)) +
  geom_col(position = "fill", color = "white") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = host_colors) +
  labs(
    title = "Blood Meal Host Composition by Subspecies",
    x = NULL, 
    y = "Proportion", 
    fill = "Host"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_host_subspecies)




# =============================================================================
# Machine Learning Approach
# =============================================================================
# =============================================================================
# CLIMATE-VECTOR DATA PREPARATION 
# =============================================================================
# =============================================================================
# STEP 1: READ AND CLEAN CLIMATE DATA
# =============================================================================

# Read climate data
ndvi <- read_csv("climate_data.csv")
temp_precip <- read_csv("temp_precip2.csv")

# Extract NDVI from ndvi dataset
ndvi_clean <- ndvi %>%
  dplyr::select(county, year, month, ndvi) %>%
  mutate(county = str_to_title(str_trim(county)))

# Extract temperature and precipitation from temp_precip dataset
temp_precip_clean <- temp_precip %>%
  dplyr::select(county, year, month, mean_temp, total_precip) %>%
  mutate(county = str_to_title(str_trim(county)))

# Check
cat("NDVI records:", nrow(ndvi_clean), "\n")
cat("NDVI years:", paste(range(ndvi_clean$year), collapse = " - "), "\n")
cat("Temp/Precip records:", nrow(temp_precip_clean), "\n")
cat("Temp/Precip years:", paste(range(temp_precip_clean$year), collapse = " - "), "\n")
cat("Anopheles years:", paste(range(anopheles_final$year, na.rm = TRUE), collapse = " - "), "\n")

# =============================================================================
# STEP 2: FIX COUNTY NAMES IN CLIMATE DATA
# =============================================================================

# Fix county names in NDVI data
ndvi_clean <- ndvi_clean %>%
  mutate(
    county = case_when(
      county == "Homabay" ~ "Homa Bay",
      county == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      county == "Tharaka-Nithi" ~ "Tharaka Nithi",
      county == "Trans-Nzoia" ~ "Trans Nzoia",
      county == "Uasin-Gishu" ~ "Uasin Gishu",
      county == "West-Pokot" ~ "West Pokot",
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Tana-River" ~ "Tana River",
      county == "Homa-Bay" ~ "Homa Bay",
      TRUE ~ county
    )
  )

# Fix county names in temp_precip data
temp_precip_clean <- temp_precip_clean %>%
  mutate(
    county = case_when(
      county == "Homabay" ~ "Homa Bay",
      county == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      county == "Tharaka-Nithi" ~ "Tharaka Nithi",
      county == "Trans-Nzoia" ~ "Trans Nzoia",
      county == "Uasin-Gishu" ~ "Uasin Gishu",
      county == "West-Pokot" ~ "West Pokot",
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Tana-River" ~ "Tana River",
      county == "Homa-Bay" ~ "Homa Bay",
      TRUE ~ county
    )
  )

# =============================================================================
# STEP 3: MERGE NDVI AND TEMP_PRECIP DATA
# =============================================================================

climate_data <- temp_precip_clean %>%
  left_join(ndvi_clean, by = c("county", "year", "month"))

# Check merge
cat("\n=== MERGED CLIMATE DATA ===\n")
glimpse(climate_data)

# Check for missing values
cat("\nMissing values in climate data:\n")
climate_data %>%
  dplyr::summarise(
    n_total = n(),
    missing_temp = sum(is.na(mean_temp)),
    missing_precip = sum(is.na(total_precip)),
    missing_ndvi = sum(is.na(ndvi))
  ) %>%
  print()

# Filter to years matching entomological data (2016-2024)
climate_data <- climate_data %>%
  filter(year >= 2016 & year <= 2024)

cat("\nClimate data records (2016-2024):", nrow(climate_data), "\n")
cat("Unique counties:", n_distinct(climate_data$county), "\n")

# =============================================================================
# STEP 4: FIX COUNTY NAMES IN ANOPHELES DATA
# =============================================================================

anopheles_final <- anopheles_final %>%
  mutate(
    county = case_when(
      county == "Homabay" ~ "Homa Bay",
      county == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      county == "Tharaka-Nithi" ~ "Tharaka Nithi",
      county == "Trans-Nzoia" ~ "Trans Nzoia",
      county == "Uasin-Gishu" ~ "Uasin Gishu",
      county == "West-Pokot" ~ "West Pokot",
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Tana-River" ~ "Tana River",
      county == "Homa-Bay" ~ "Homa Bay",
      TRUE ~ county
    )
  )

# Verify county names match
cat("Counties in anopheles but NOT in climate:\n")
missing_counties <- setdiff(unique(anopheles_final$county), unique(climate_data$county))
print(missing_counties)

cat("\nCounties in climate but NOT in anopheles:\n")
extra_counties <- setdiff(unique(climate_data$county), unique(anopheles_final$county))
print(extra_counties)

# =============================================================================
# STEP 5: CREATE LAGGED VARIABLES FROM CLIMATE DATA
# =============================================================================

climate_with_lags <- climate_data %>%
  arrange(county, year, month) %>%
  group_by(county) %>%
  mutate(
    # Temperature lags
    mean_temp_lag1 = lag(mean_temp, 1),
    mean_temp_lag2 = lag(mean_temp, 2),
    mean_temp_lag3 = lag(mean_temp, 3),
    
    # Precipitation lags
    total_precip_lag1 = lag(total_precip, 1),
    total_precip_lag2 = lag(total_precip, 2),
    total_precip_lag3 = lag(total_precip, 3),
    
    # NDVI lags
    ndvi_lag1 = lag(ndvi, 1),
    ndvi_lag2 = lag(ndvi, 2),
    ndvi_lag3 = lag(ndvi, 3)
  ) %>%
  ungroup()

# Check lags
cat("\n=== CLIMATE DATA WITH LAGS ===\n")
glimpse(climate_with_lags)

# Verify lags are correct (sample check)
cat("\nSample lag verification (Kisumu):\n")
climate_with_lags %>%
  filter(county == "Kisumu", year == 2020) %>%
  dplyr::select(county, year, month, mean_temp, mean_temp_lag1, mean_temp_lag2, mean_temp_lag3) %>%
  print(n = 12)

# =============================================================================
# STEP 6: JOIN CLIMATE TO INDIVIDUAL MOSQUITO RECORDS
# =============================================================================
cat("Individual mosquito records:", nrow(anopheles_final), "\n")

# Join climate data to individual mosquito records
# Use month_num from anopheles_final to match with month in climate data
ml_data <- anopheles_final %>%
  left_join(climate_with_lags, by = c("county", "year", "month_num" = "month"))

# Check merge
cat("Total records after merge:", nrow(ml_data), "\n")
cat("Records WITH climate data:", sum(!is.na(ml_data$mean_temp)), "\n")
cat("Records MISSING climate data:", sum(is.na(ml_data$mean_temp)), "\n")
cat("Percent with climate:", round(mean(!is.na(ml_data$mean_temp)) * 100, 1), "%\n")

# Check which records are missing climate
ml_data %>%
  group_by(year) %>%
  dplyr::summarise(
    total = n(),
    missing_climate = sum(is.na(mean_temp)),
    pct_missing = round(missing_climate / total * 100, 1)
  ) %>%
  print()

# Check by county
ml_data %>%
  filter(is.na(mean_temp)) %>%
  dplyr::count(county, sort = TRUE) %>%
  print(n = 20)

# this is because we had a month labelled as NA in the anopheles data, so no problem


# =============================================================================
# STEP 7: ADD SEASON VARIABLE
# =============================================================================

ml_data <- ml_data %>%
  mutate(
    season = case_when(
      month_num %in% c(1, 2) ~ "Short Dry (Jan-Feb)",
      month_num %in% c(3, 4, 5) ~ "Long Rains (Mar-May)",
      month_num %in% c(6, 7, 8, 9) ~ "Long Dry (Jun-Sep)",
      month_num %in% c(10, 11, 12) ~ "Short Rains (Oct-Dec)"
    ),
    season = factor(season, levels = c(
      "Short Dry (Jan-Feb)",
      "Long Rains (Mar-May)",
      "Long Dry (Jun-Sep)",
      "Short Rains (Oct-Dec)"
    ))
  )

# =============================================================================
# STEP 8: PREPARE OUTCOME VARIABLES
# =============================================================================

ml_data <- ml_data %>%
  mutate(
    # Sporozoite outcome (binary)
    spz_outcome = case_when(
      sporozoite_status == "Positive" ~ "Positive",
      sporozoite_status == "Negative" ~ "Negative",
      TRUE ~ NA_character_
    ),
    spz_outcome = factor(spz_outcome, levels = c("Negative", "Positive")),
    spz_numeric = ifelse(spz_outcome == "Positive", 1, 0),
    
    # Species outcome
    species_outcome = case_when(
      !is.na(species) & species != "Unknown" ~ species,
      TRUE ~ NA_character_
    ),
    species_outcome = factor(species_outcome),
    
    # Subspecies outcome
    subspecies_outcome = case_when(
      !is.na(subspecies)  ~ subspecies,
      TRUE ~ NA_character_
    ),
    subspecies_outcome = factor(subspecies_outcome),
    
    # Blood meal outcome (human vs non-human)
    bloodmeal_outcome = case_when(
      human_bloodmeal == 1 ~ "Human",
      human_bloodmeal == 0 ~ "Non-Human",
      TRUE ~ NA_character_
    ),
    bloodmeal_outcome = factor(bloodmeal_outcome, levels = c("Non-Human", "Human")),
    
    # Malaria burden as factor
    malaria_burden = factor(malaria_burden)
  )

# Save the prepared data
write_csv(ml_data, "ml_data_individual.csv")








# =============================================================================
# EXPLORATORY ANALYSIS: CLIMATE-VECTOR RELATIONSHIPS 
# =============================================================================
library(tidyverse)
library(corrplot)
library(patchwork)
library(caret) # For findCorrelation

# =============================================================================
# DATA SUMMARY
# =============================================================================
# Filter to records with climate data
ml_data_complete <- ml_data %>%
  filter(!is.na(mean_temp), !is.na(total_precip))

cat("Records for analysis (with climate):", nrow(ml_data_complete), "\n\n")

# =============================================================================
# OUTCOME VARIABLE SUMMARIES
# =============================================================================

# Sporozoite status
ml_data_complete %>%
  filter(!is.na(sporozoite_status)) %>%
  dplyr::count(sporozoite_status) %>%
  mutate(percent = round(n / sum(n) * 100, 2)) %>%
  print()

# Species
ml_data_complete %>%
  filter(!is.na(species)) %>%
  dplyr::count(species) %>%
  mutate(percent = round(n / sum(n) * 100, 2)) %>%
  print()

# Subspecies
ml_data_complete %>%
  filter(!is.na(subspecies)) %>%
  dplyr::count(subspecies, sort = TRUE) %>%
  mutate(percent = round(n / sum(n) * 100, 2)) %>%
  print()

# Blood meal
# Blood Meal (Human vs Non-Human)
ml_data_complete %>%
  filter(!is.na(bloodmeal_outcome)) %>%
  dplyr::count(bloodmeal_outcome) %>%
  mutate(percent = round(n / sum(n) * 100, 2)) %>%
  print()

# =============================================================================
# CLIMATE VARIABLE SUMMARIES
# =============================================================================

ml_data_complete %>%
  dplyr::summarise(
    n = n(),
    mean_temp_avg = round(mean(mean_temp, na.rm = TRUE), 2),
    mean_temp_sd = round(sd(mean_temp, na.rm = TRUE), 2),
    mean_temp_range = paste(round(min(mean_temp, na.rm = TRUE), 1), "-", 
                            round(max(mean_temp, na.rm = TRUE), 1)),
    precip_avg = round(mean(total_precip, na.rm = TRUE), 2),
    precip_sd = round(sd(total_precip, na.rm = TRUE), 2),
    ndvi_avg = round(mean(ndvi, na.rm = TRUE), 3),
    ndvi_sd = round(sd(ndvi, na.rm = TRUE), 3)
  ) %>%
  print()

# =============================================================================
# CLIMATE BY SPOROZOITE STATUS
# =============================================================================

ml_data_complete %>%
  filter(sporozoite_status %in% c("Positive", "Negative")) %>%
  group_by(sporozoite_status) %>%
  dplyr::summarise(
    n = n(),
    mean_temp = round(mean(mean_temp, na.rm = TRUE), 2),
    mean_precip = round(mean(total_precip, na.rm = TRUE), 2),
    mean_ndvi = round(mean(ndvi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# CLIMATE BY SPECIES
# =============================================================================

ml_data_complete %>%
  filter(!is.na(species)) %>%
  group_by(species) %>%
  dplyr::summarise(
    n = n(),
    mean_temp = round(mean(mean_temp, na.rm = TRUE), 2),
    mean_precip = round(mean(total_precip, na.rm = TRUE), 2),
    mean_ndvi = round(mean(ndvi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# CLIMATE BY SUBSPECIES
# =============================================================================
ml_data_complete %>%
  filter(!is.na(subspecies)) %>%
  group_by(subspecies) %>%
  dplyr::summarise(
    n = n(),
    mean_temp = round(mean(mean_temp, na.rm = TRUE), 2),
    mean_precip = round(mean(total_precip, na.rm = TRUE), 2),
    mean_ndvi = round(mean(ndvi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n)) %>%
  print()

# =============================================================================
# AGGREGATE DATA FOR CORRELATION ANALYSIS
# =============================================================================

# Create monthly aggregates for correlation analysis
monthly_summary <- ml_data_complete %>%
  group_by(county, year, month_num, season, malaria_burden) %>%
  dplyr::summarise(
    # Counts
    total_count = n(),
    n_positive_spz = sum(spz_numeric, na.rm = TRUE),
    n_tested_spz = sum(!is.na(spz_numeric)),
    spz_rate = ifelse(n_tested_spz > 0, n_positive_spz / n_tested_spz * 100, NA),
    
    n_human_blood = sum(human_bloodmeal, na.rm = TRUE),
    n_bloodmeal = sum(!is.na(human_bloodmeal)),
    hbi = ifelse(n_bloodmeal > 0, n_human_blood / n_bloodmeal * 100, NA),
    
    # Climate (take mean since all records in same month have same climate)
    mean_temp = mean(mean_temp, na.rm = TRUE),
    total_precip = mean(total_precip, na.rm = TRUE),
    ndvi = mean(ndvi, na.rm = TRUE),
    mean_temp_lag1 = mean(mean_temp_lag1, na.rm = TRUE),
    mean_temp_lag2 = mean(mean_temp_lag2, na.rm = TRUE),
    mean_temp_lag3 = mean(mean_temp_lag3, na.rm = TRUE),
    total_precip_lag1 = mean(total_precip_lag1, na.rm = TRUE),
    total_precip_lag2 = mean(total_precip_lag2, na.rm = TRUE),
    total_precip_lag3 = mean(total_precip_lag3, na.rm = TRUE),
    ndvi_lag1 = mean(ndvi_lag1, na.rm = TRUE),
    ndvi_lag2 = mean(ndvi_lag2, na.rm = TRUE),
    ndvi_lag3 = mean(ndvi_lag3, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("Monthly records:", nrow(monthly_summary), "\n")
glimpse(monthly_summary)

# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================
# Prepare correlation data
cor_data <- monthly_summary %>%
  dplyr::select(
    total_count, spz_rate, hbi,
    mean_temp, total_precip, ndvi,
    mean_temp_lag1, mean_temp_lag2, mean_temp_lag3,
    total_precip_lag1, total_precip_lag2, total_precip_lag3,
    ndvi_lag1, ndvi_lag2, ndvi_lag3
  ) %>%
  drop_na()

cat("Records for correlation:", nrow(cor_data), "\n")

# Calculate correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs")

# Visualize correlation matrix
corrplot(cor_matrix, 
         method = "color", 
         type = "upper",
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.7,
         addCoef.col = "black", 
         number.cex = 0.5,
         col = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
         title = "Climate-Vector Correlations (Monthly Aggregated)",
         mar = c(0, 0, 2, 0))
# dev.off()


# --- Remove highly correlated variables ---
#  we define a correlation threshold
# filter out correlated variables.
correlation_threshold <- 0.80 

# Use caret::findCorrelation to identify variables to remove
# This function returns the *indices* of the columns to remove.
high_corr_vars_to_remove_indices <- findCorrelation(cor_matrix, cutoff = correlation_threshold, exact = FALSE)

# Convert indices to variable names
if (length(high_corr_vars_to_remove_indices) > 0) {
  high_corr_vars_to_remove <- colnames(cor_matrix)[high_corr_vars_to_remove_indices]
  cat("Variables identified for removal (absolute correlation >", correlation_threshold, "):\n")
  print(high_corr_vars_to_remove)
} else {
  high_corr_vars_to_remove <- character(0)
  cat("No variables found with absolute correlation above", correlation_threshold, "to remove.\n")
}

# Create a new data frame with the selected variables
# we select all columns that are NOT in the list of variables to remove.
all_original_vars <- colnames(cor_data)
vars_to_keep <- setdiff(all_original_vars, high_corr_vars_to_remove)

cor_data_filtered <- cor_data %>%
  dplyr::select(all_of(vars_to_keep))

cat("\nOriginal number of variables:", ncol(cor_data), "\n")
cat("Number of variables after filtering:", ncol(cor_data_filtered), "\n")
cat("Remaining variables:\n")
print(colnames(cor_data_filtered))

# Re-calculate and visualize the correlation matrix of the filtered data
if (ncol(cor_data_filtered) > 1) {
  cat("\n--- Correlation matrix of filtered variables ---\n")
  cor_matrix_filtered <- cor(cor_data_filtered, use = "complete.obs")
  
  corrplot(cor_matrix_filtered,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.7,
           addCoef.col = "black",
           number.cex = 0.5,
           col = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
           title = paste0("Filtered Climate-Vector Correlations (Cutoff ", correlation_threshold, ")"),
           mar = c(0, 0, 2, 0))
} else if (ncol(cor_data_filtered) == 1) {
  cat("\nOnly one variable remains after filtering; correlation matrix not applicable for a single variable.\n")
} else {
  cat("\nNo variables remain after filtering.\n")
}



# =============================================================================
# CLIMATE-VECTOR RELATIONSHIP PLOTS
# =============================================================================

# Temperature vs Mosquito Count (monthly)
p1 <- monthly_summary %>%
  ggplot(aes(x = mean_temp, y = total_count)) +
  geom_point(aes(color = season), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Temperature vs Mosquito Count",
    x = "Mean Temperature (°C)",
    y = "Monthly Count",
    color = "Season"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p1

# Precipitation vs Mosquito Count
p2 <- monthly_summary %>%
  ggplot(aes(x = total_precip, y = total_count)) +
  geom_point(aes(color = season), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Precipitation vs Mosquito Count",
    x = "Total Precipitation (mm)",
    y = "Monthly Count",
    color = "Season"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p2

# NDVI vs Mosquito Count
p3 <- monthly_summary %>%
  filter(!is.na(ndvi)) %>%
  ggplot(aes(x = ndvi, y = total_count)) +
  geom_point(aes(color = season), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "NDVI vs Mosquito Count",
    x = "NDVI",
    y = "Monthly Count",
    color = "Season"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p3

# Combine plots
(p1 + p2 + p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Climate-Vector Relationships (Monthly Aggregated Data)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# =============================================================================
# CLIMATE VS SPOROZOITE RATE
# =============================================================================

# Temperature vs Sporozoite Rate
p4 <- monthly_summary %>%
  filter(!is.na(spz_rate)) %>%
  ggplot(aes(x = mean_temp, y = spz_rate)) +
  geom_point(aes(color = malaria_burden), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_manual(values = zone_colors) +
  labs(
    title = "Temperature vs Sporozoite Rate",
    x = "Mean Temperature (°C)",
    y = "Sporozoite Rate (%)",
    color = "Malaria Burden"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p4

# Precipitation vs Sporozoite Rate
p5 <- monthly_summary %>%
  filter(!is.na(spz_rate)) %>%
  ggplot(aes(x = total_precip, y = spz_rate)) +
  geom_point(aes(color = malaria_burden), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_manual(values = zone_colors) +
  labs(
    title = "Precipitation vs Sporozoite Rate",
    x = "Total Precipitation (mm)",
    y = "Sporozoite Rate (%)",
    color = "Malaria Burden"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p5

# NDVI vs Sporozoite Rate
p6 <- monthly_summary %>%
  filter(!is.na(spz_rate), !is.na(ndvi)) %>%
  ggplot(aes(x = ndvi, y = spz_rate)) +
  geom_point(aes(color = malaria_burden), alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  scale_color_manual(values = zone_colors) +
  labs(
    title = "NDVI vs Sporozoite Rate",
    x = "NDVI",
    y = "Sporozoite Rate (%)",
    color = "Malaria Burden"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p6

# Combine
(p4 + p5 + p6) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Climate vs Sporozoite Infection Rate",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# =============================================================================
# LAG EFFECT ANALYSIS
# =============================================================================

# Compare current vs lagged precipitation effects on mosquito count
lag_comparison <- monthly_summary %>%
  dplyr::select(total_count, total_precip, total_precip_lag1, 
                total_precip_lag2, total_precip_lag3) %>%
  pivot_longer(cols = starts_with("total_precip"),
               names_to = "lag_type", values_to = "precipitation") %>%
  mutate(
    lag = case_when(
      lag_type == "total_precip" ~ "Current",
      lag_type == "total_precip_lag1" ~ "1 Month Lag",
      lag_type == "total_precip_lag2" ~ "2 Month Lag",
      lag_type == "total_precip_lag3" ~ "3 Month Lag"
    ),
    lag = factor(lag, levels = c("Current", "1 Month Lag", "2 Month Lag", "3 Month Lag"))
  )

p_lag_precip <- lag_comparison %>%
  ggplot(aes(x = precipitation, y = total_count)) +
  geom_point(alpha = 0.5, color = "#4575b4") +
  geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
  facet_wrap(~lag, ncol = 2) +
  labs(
    title = "Effect of Precipitation Lag on Mosquito Count",
    x = "Precipitation (mm)",
    y = "Monthly Count"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold")
  )

print(p_lag_precip)

# Temperature lag effects
lag_temp <- monthly_summary %>%
  dplyr::select(total_count, mean_temp, mean_temp_lag1, 
                mean_temp_lag2, mean_temp_lag3) %>%
  pivot_longer(cols = starts_with("mean_temp"),
               names_to = "lag_type", values_to = "temperature") %>%
  mutate(
    lag = case_when(
      lag_type == "mean_temp" ~ "Current",
      lag_type == "mean_temp_lag1" ~ "1 Month Lag",
      lag_type == "mean_temp_lag2" ~ "2 Month Lag",
      lag_type == "mean_temp_lag3" ~ "3 Month Lag"
    ),
    lag = factor(lag, levels = c("Current", "1 Month Lag", "2 Month Lag", "3 Month Lag"))
  )

p_lag_temp <- lag_temp %>%
  ggplot(aes(x = temperature, y = total_count)) +
  geom_point(alpha = 0.5, color = "#d73027") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~lag, ncol = 2) +
  labs(
    title = "Effect of Temperature Lag on Mosquito Count",
    x = "Temperature (°C)",
    y = "Monthly Count"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(p_lag_temp)

# NDVI lag effects
lag_ndvi <- monthly_summary %>%
  dplyr::select(total_count, ndvi, ndvi_lag1, ndvi_lag2, ndvi_lag3) %>%
  pivot_longer(cols = starts_with("ndvi"),
               names_to = "lag_type", values_to = "ndvi_value") %>%
  mutate(
    lag = case_when(
      lag_type == "ndvi" ~ "Current",
      lag_type == "ndvi_lag1" ~ "1 Month Lag",
      lag_type == "ndvi_lag2" ~ "2 Month Lag",
      lag_type == "ndvi_lag3" ~ "3 Month Lag"
    ),
    lag = factor(lag, levels = c("Current", "1 Month Lag", "2 Month Lag", "3 Month Lag"))
  )

p_lag_ndvi <- lag_ndvi %>%
  filter(!is.na(ndvi_value)) %>%
  ggplot(aes(x = ndvi_value, y = total_count)) +
  geom_point(alpha = 0.5, color = "#1a9850") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~lag, ncol = 2) +
  labs(
    title = "Effect of NDVI Lag on Mosquito Count",
    x = "NDVI",
    y = "Monthly Count"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(p_lag_ndvi)

# =============================================================================
# CLIMATE BY INDIVIDUAL OUTCOMES (BOXPLOTS)
# =============================================================================

# Temperature distribution by sporozoite status
p_temp_spz <- ml_data_complete %>%
  filter(sporozoite_status %in% c("Positive", "Negative")) %>%
  ggplot(aes(x = sporozoite_status, y = mean_temp, fill = sporozoite_status)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Negative" = "#2166ac", "Positive" = "#d73027")) +
  labs(
    title = "Temperature by Sporozoite Status",
    x = "Sporozoite Status",
    y = "Mean Temperature (°C)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p_temp_spz

# Precipitation by sporozoite status
p_precip_spz <- ml_data_complete %>%
  filter(sporozoite_status %in% c("Positive", "Negative")) %>%
  ggplot(aes(x = sporozoite_status, y = total_precip, fill = sporozoite_status)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Negative" = "#2166ac", "Positive" = "#d73027")) +
  labs(
    title = "Precipitation by Sporozoite Status",
    x = "Sporozoite Status",
    y = "Total Precipitation (mm)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p_precip_spz


# NDVI by sporozoite status
p_ndvi_spz <- ml_data_complete %>%
  filter(sporozoite_status %in% c("Positive", "Negative"), !is.na(ndvi)) %>%
  ggplot(aes(x = sporozoite_status, y = ndvi, fill = sporozoite_status)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Negative" = "#2166ac", "Positive" = "#d73027")) +
  labs(
    title = "NDVI by Sporozoite Status",
    x = "Sporozoite Status",
    y = "NDVI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )
p_ndvi_spz

# Combine
(p_temp_spz + p_precip_spz + p_ndvi_spz) +
  plot_annotation(
    title = "Climate Conditions by Sporozoite Infection Status",
    subtitle = "Individual mosquito records",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# =============================================================================
# CLIMATE BY SUBSPECIES
# =============================================================================
# Temperature by subspecies
p_temp_subspecies <- ml_data_complete %>%
  filter(!is.na(subspecies)) %>%
  ggplot(aes(x = reorder(subspecies, mean_temp, FUN = median), y = mean_temp, fill = subspecies)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Temperature Distribution by Subspecies",
    x = NULL,
    y = "Mean Temperature (°C)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  coord_flip()

print(p_temp_subspecies)

# Precipitation by subspecies
p_precip_subspecies <- ml_data_complete %>%
  filter(!is.na(subspecies)) %>%
  ggplot(aes(x = reorder(subspecies, total_precip, FUN = median), y = total_precip, fill = subspecies)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Precipitation Distribution by Subspecies",
    x = NULL,
    y = "Total Precipitation (mm)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  coord_flip()

print(p_precip_subspecies)

# =============================================================================
# STEP 14: STATISTICAL TESTS
# =============================================================================

# T-test for temperature by sporozoite status
cat("T-test: Temperature by Sporozoite Status\n")
t_test_temp <- t.test(mean_temp ~ sporozoite_status, 
                      data = ml_data_complete %>% 
                        filter(sporozoite_status %in% c("Positive", "Negative")))
print(t_test_temp)

# T-test for precipitation by sporozoite status
cat("\nT-test: Precipitation by Sporozoite Status\n")
t_test_precip <- t.test(total_precip ~ sporozoite_status, 
                        data = ml_data_complete %>% 
                          filter(sporozoite_status %in% c("Positive", "Negative")))
print(t_test_precip)

# T-test for NDVI by sporozoite status
cat("\nT-test: NDVI by Sporozoite Status\n")
t_test_ndvi <- t.test(ndvi ~ sporozoite_status, 
                      data = ml_data_complete %>% 
                        filter(sporozoite_status %in% c("Positive", "Negative"),
                               !is.na(ndvi)))
print(t_test_ndvi)





# =============================================================================
# INTERPRETABLE MACHINE LEARNING ANALYSIS FOR MALARIA VECTORS
# (Using Individual Mosquito Data for Sporozoite Prediction)
# =============================================================================

# Load required packages
library(tidyverse)
library(randomForest)
library(xgboost)
library(lightgbm)
library(caret)
library(pROC)
library(PRROC)
library(sf)
library(viridis)
library(patchwork)
library(scales)

# Set seed for reproducibility
set.seed(123)

# Load Kenya county shapefile (adjust path as needed)
kenya_counties <- st_read("C:/LENSON/PHD/PhDData Analysis/IML_Lenson/GADM/gadm41_KEN_1.shp")
# Check shapefile
glimpse(kenya_counties)
plot(kenya_counties["geometry"])

# Standardize county names to match your data
kenya_counties <- kenya_counties %>%
  mutate(county = str_to_title(str_trim(NAME_1)))  

# Check county names match
cat("Counties in shapefile:\n")
print(sort(unique(kenya_counties$county)))

cat("\nCounties in vector data:\n")
print(sort(unique(ml_data$county)))

# Vector of corrections: shapefile names -> vector data names
county_recode <- c(
  "Elgeyo-Marakwet" = "Elgeyo Marakwet",
  "Murang'a" = "Muranga",
  "Tharaka-Nithi" = "Tharaka Nithi"
)

kenya_counties <- kenya_counties %>%
  mutate(county = recode(county, !!!county_recode))


# -----------------------------------------------------------------------------
# 1. PREPARE DATA FOR MACHINE LEARNING
# -----------------------------------------------------------------------------

# Use the individual data (ml_data) prepared earlier
# Filter for sporozoite analysis (only tested specimens)
spz_ml_data <- ml_data %>%
  filter(sporozoite_status %in% c("Positive", "Negative")) %>%
  filter(!is.na(mean_temp), !is.na(total_precip), !is.na(ndvi)) %>%
  mutate(
    # Outcome: Sporozoite Status
    spz_outcome = factor(ifelse(sporozoite_status == "Positive", "Positive", "Negative"),
                         levels = c("Negative", "Positive")),
    spz_numeric = ifelse(spz_outcome == "Positive", 1, 0)
  )

cat("\n=== ML DATA SUMMARY (INDIVIDUAL RECORDS) ===\n")
cat("Total records:", nrow(spz_ml_data), "\n")
cat("Counties:", n_distinct(spz_ml_data$county), "\n")
cat("Positive cases:", sum(spz_ml_data$spz_outcome == "Positive"), "\n")
cat("Negative cases:", sum(spz_ml_data$spz_outcome == "Negative"), "\n")
cat("Positive rate:", round(mean(spz_ml_data$spz_outcome == "Positive") * 100, 2), "%\n")

# Check class balance
cat("\n=== CLASS DISTRIBUTION ===\n")
table(spz_ml_data$spz_outcome)

# -----------------------------------------------------------------------------
# 2. DEFINE FEATURE SETS
# -----------------------------------------------------------------------------

# Climate features (current)
climate_current <- c("mean_temp", "total_precip", "ndvi")

# Climate features (lagged)
climate_lag1 <- c("mean_temp_lag1", "total_precip_lag1", "mean_ndvi_lag1") # Adjust names if needed (e.g., ndvi_lag1)
climate_lag2 <- c("mean_temp_lag2", "total_precip_lag2", "mean_ndvi_lag2")
climate_lag3 <- c("mean_temp_lag3", "total_precip_lag3", "mean_ndvi_lag3")

# Temporal features
temporal <- c("month_num") # Or "month"

# All features combined
all_features <- c(climate_current, climate_lag1, climate_lag2, climate_lag3, temporal)

# Check available features
available_features <- all_features[all_features %in% names(spz_ml_data)]
cat("\n=== FEATURES ===\n")
cat("Total features used:", length(available_features), "\n")
print(available_features)

# -----------------------------------------------------------------------------
# 3. TRAIN/TEST SPLIT (STRATIFIED)
# -----------------------------------------------------------------------------

# Create stratified train/test split (80/20)
train_index <- createDataPartition(spz_ml_data$spz_outcome, p = 0.8, list = FALSE)
train_data <- spz_ml_data[train_index, ]
test_data <- spz_ml_data[-train_index, ]

cat("\n=== TRAIN/TEST SPLIT ===\n")
cat("Training samples:", nrow(train_data), "\n")
cat("Testing samples:", nrow(test_data), "\n")

cat("\nTraining class distribution:\n")
print(table(train_data$spz_outcome))

cat("\nTesting class distribution:\n")
print(table(test_data$spz_outcome))

# Prepare feature matrices
X_train <- train_data %>% dplyr::select(all_of(available_features)) %>% as.data.frame()
y_train <- train_data$spz_outcome
y_train_num <- train_data$spz_numeric

X_test <- test_data %>% dplyr::select(all_of(available_features)) %>% as.data.frame()
y_test <- test_data$spz_outcome
y_test_num <- test_data$spz_numeric

# Handle class imbalance (calculate scale_pos_weight)
scale_pos_weight <- sum(y_train == "Negative") / sum(y_train == "Positive")
cat("Scale positive weight (for imbalance):", round(scale_pos_weight, 2), "\n")

# Remove NAs if any remain
complete_train <- complete.cases(X_train)
complete_test <- complete.cases(X_test)

X_train <- X_train[complete_train, ]
y_train <- y_train[complete_train]
y_train_num <- y_train_num[complete_train]

X_test <- X_test[complete_test, ]
y_test <- y_test[complete_test]
y_test_num <- y_test_num[complete_test]

cat("\nAfter removing NAs:\n")
cat("Training samples:", nrow(X_train), "\n")
cat("Testing samples:", nrow(X_test), "\n")

# =============================================================================
# MODEL TRAINING
# =============================================================================

# -----------------------------------------------------------------------------
# 4.1 RANDOM FOREST
# -----------------------------------------------------------------------------

rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  mtry = floor(sqrt(ncol(X_train))),
  importance = TRUE,
  classwt = c("Negative" = 1, "Positive" = scale_pos_weight), # Handle imbalance
  na.action = na.omit
)

print(rf_model)

# Predictions
rf_pred_prob <- predict(rf_model, X_test, type = "prob")[, "Positive"]
rf_pred_class <- predict(rf_model, X_test)

# -----------------------------------------------------------------------------
# 4.2 XGBOOST
# -----------------------------------------------------------------------------

# Prepare data for XGBoost
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train_num)
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test_num)

# XGBoost parameters
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 6,
  subsample = 0.8,
  colsample_bytree = 0.8,
  scale_pos_weight = scale_pos_weight # Handle imbalance
)

# Train with cross-validation to find optimal rounds
xgb_cv <- xgb.cv(
  params = xgb_params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 20,
  verbose = 0
)

best_nrounds <- xgb_cv$best_iteration
cat("Best number of rounds:", best_nrounds, "\n")

# Train final model
xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)
# Use 100 rounds if the CV fails to provide a number
final_nrounds <- if(length(best_nrounds) > 0) best_nrounds else 100

xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = final_nrounds,
  verbose = 0
)

# Predictions
xgb_pred_prob <- predict(xgb_model, dtest)
xgb_pred_class <- factor(ifelse(xgb_pred_prob > 0.5, "Positive", "Negative"),
                         levels = c("Negative", "Positive"))

# -----------------------------------------------------------------------------
# 4.3 LIGHTGBM
# -----------------------------------------------------------------------------

# Prepare data for LightGBM
lgb_train <- lgb.Dataset(data = as.matrix(X_train), label = y_train_num)
lgb_test <- lgb.Dataset(data = as.matrix(X_test), label = y_test_num, reference = lgb_train)

# LightGBM parameters
lgb_params <- list(
  objective = "binary",
  metric = "auc",
  learning_rate = 0.1,
  num_leaves = 31,
  max_depth = 6,
  scale_pos_weight = scale_pos_weight, # Handle imbalance
  verbose = -1
)

# Train with early stopping
lgb_model <- lgb.train(
  params = lgb_params,
  data = lgb_train,
  nrounds = 500,
  valids = list(test = lgb_test),
  early_stopping_rounds = 20,
  verbose = -1
)

cat("Best iteration:", lgb_model$best_iter, "\n")

# Predictions
lgb_pred_prob <- predict(lgb_model, as.matrix(X_test))
lgb_pred_class <- factor(ifelse(lgb_pred_prob > 0.5, "Positive", "Negative"),
                         levels = c("Negative", "Positive"))

# =============================================================================
# MODEL EVALUATION
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1 METRICS FUNCTION 
# -----------------------------------------------------------------------------

calculate_metrics <- function(actual, predicted_class, predicted_prob, model_name) {
  
  # Use caret::confusionMatrix explicitly
  cm <- caret::confusionMatrix(predicted_class, actual, positive = "Positive")
  
  accuracy <- cm$overall["Accuracy"]
  precision <- cm$byClass["Precision"]
  recall <- cm$byClass["Recall"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  f1 <- cm$byClass["F1"]
  
  # AUC-ROC
  roc_obj <- pROC::roc(as.numeric(actual == "Positive"), predicted_prob, quiet = TRUE)
  auc_roc <- pROC::auc(roc_obj)
  
  # AUC-PR
  pr_obj <- PRROC::pr.curve(scores.class0 = predicted_prob[actual == "Positive"],
                            scores.class1 = predicted_prob[actual == "Negative"],
                            curve = FALSE)
  auc_pr <- pr_obj$auc.integral
  
  tibble(
    Model = model_name,
    Accuracy = round(accuracy, 4),
    # Precision = round(precision, 4),
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    Recall = round(recall, 4),
    F1_Score = round(f1, 4),
    AUC_ROC = round(as.numeric(auc_roc), 4),
    AUC_PR = round(auc_pr, 4)
  )
}

# -----------------------------------------------------------------------------
# 5.2 CALCULATE AND COMPARE METRICS
# -----------------------------------------------------------------------------

rf_metrics <- calculate_metrics(y_test, rf_pred_class, rf_pred_prob, "Random Forest")
xgb_metrics <- calculate_metrics(y_test, xgb_pred_class, xgb_pred_prob, "XGBoost")
lgb_metrics <- calculate_metrics(y_test, lgb_pred_class, lgb_pred_prob, "LightGBM")

all_metrics <- bind_rows(rf_metrics, xgb_metrics, lgb_metrics)

cat("\n=== MODEL COMPARISON ===\n")
print(all_metrics)

# -----------------------------------------------------------------------------
# 5.3 VISUALIZE COMPARISON
# -----------------------------------------------------------------------------

metrics_long <- all_metrics %>%
  pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Accuracy", "Precision", "Recall", 
                                            "Specificity", "F1_Score", "AUC_ROC", "AUC_PR")))

p_metrics <- metrics_long %>%
  ggplot(aes(x = Metric, y = Value, fill = Model)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Value, 2)), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Random Forest" = "#2166ac", 
                               "XGBoost" = "#d73027", 
                               "LightGBM" = "#1a9850")) +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(
    title = "Model Performance Comparison",
    subtitle = "Sporozoite Positivity Prediction (Individual Level)",
    x = NULL, y = "Score"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

print(p_metrics)

# -----------------------------------------------------------------------------
# 5.4 ROC CURVES 
# -----------------------------------------------------------------------------

# Calculate ROC for each model
roc_rf <- pROC::roc(response = y_test, predictor = rf_pred_prob, levels = c("Negative", "Positive"), quiet = TRUE)
roc_xgb <- pROC::roc(response = y_test, predictor = xgb_pred_prob, levels = c("Negative", "Positive"), quiet = TRUE)
roc_lgb <- pROC::roc(response = y_test, predictor = lgb_pred_prob, levels = c("Negative", "Positive"), quiet = TRUE)

# Apply statistical smoothing to each ROC object
smooth_roc_rf <- smooth(roc_rf)
smooth_roc_xgb <- smooth(roc_xgb)
smooth_roc_lgb <- smooth(roc_lgb)

# Create ROC data for plotting using the SMOOTHED coordinates
smooth_roc_data <- bind_rows(
  tibble(
    FPR = 1 - smooth_roc_rf$specificities,
    TPR = smooth_roc_rf$sensitivities,
    # Use the ORIGINAL AUC for the label, as smoothing is for visualization
    Model = paste0("Random Forest (AUC = ", round(pROC::auc(roc_rf), 3), ")")
  ),
  tibble(
    FPR = 1 - smooth_roc_xgb$specificities,
    TPR = smooth_roc_xgb$sensitivities,
    Model = paste0("XGBoost (AUC = ", round(pROC::auc(roc_xgb), 3), ")")
  ),
  tibble(
    FPR = 1 - smooth_roc_lgb$specificities,
    TPR = smooth_roc_lgb$sensitivities,
    Model = paste0("LightGBM (AUC = ", round(pROC::auc(roc_lgb), 3), ")")
  )
)

# Smoothed ROC curve plot
p_roc_smooth <- smooth_roc_data %>%
  ggplot(aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("#2166ac", "#d73027", "#1a9850")) +
  labs(
    title = "ROC Curves - Sporozoite Prediction",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  coord_equal()

print(p_roc_smooth)


# -----------------------------------------------------------------------------
# 4.5 CONFUSION MATRICES
# -----------------------------------------------------------------------------

# Function to plot confusion matrix
plot_confusion_matrix <- function(actual, predicted, model_name, positive_class = "Positive") {
  
  # Use caret::confusionMatrix explicitly and set correct positive class
  cm <- caret::confusionMatrix(predicted, actual, positive = positive_class)
  
  cm_table <- as.data.frame(cm$table)
  names(cm_table) <- c("Predicted", "Actual", "Freq")
  
  # Calculate percentages by actual class
  cm_table <- cm_table %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(Percent = Freq / sum(Freq) * 100) %>%
    dplyr::ungroup()
  
  ggplot(cm_table, aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = paste0(Freq, "\n(", round(Percent, 1), "%)")), 
              color = "white", fontface = "bold", size = 4) +
    scale_fill_gradient(low = "#fee08b", high = "#d73027") +
    labs(
      title = paste("Confusion Matrix -", model_name),
      x = "Actual",
      y = "Predicted",
      fill = "Count"
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    coord_equal()
}
p_cm_rf  <- plot_confusion_matrix(y_test, rf_pred_class,  "Random Forest", positive_class = "Positive")
p_cm_xgb <- plot_confusion_matrix(y_test, xgb_pred_class, "XGBoost",       positive_class = "Positive")
p_cm_lgb <- plot_confusion_matrix(y_test, lgb_pred_class, "LightGBM",      positive_class = "Positive")

(p_cm_rf + p_cm_xgb + p_cm_lgb) +
  plot_annotation(title = "Confusion Matrices - Sporozoite Prediction")

# Combine confusion matrices
(p_cm_rf + p_cm_xgb + p_cm_lgb) +
  plot_annotation(title = "Confusion Matrices - All Models")






# =============================================================================
# CLIMATE–VECTOR TIME SERIES MODELLING (MONTHLY AGGREGATED)
# =============================================================================
library(mgcv)  # For GAMs

# -----------------------------------------------------------------------------
# 3.1 GENERALIZED ADDITIVE MODELS (GAMs)
# -----------------------------------------------------------------------------

# Prepare data: we use monthly_summary created from individual records
# (monthly_summary should have: total_count, mean_temp, total_precip, ndvi, lags, month_num, malaria_burden)
gam_data <- monthly_summary %>%
  filter(!is.na(mean_temp), !is.na(total_precip), !is.na(ndvi)) %>%
  mutate(
    month = month_num,
    malaria_burden = factor(malaria_burden)
  )

cat("\n=== GAM DATA SUMMARY ===\n")
cat("Monthly records:", nrow(gam_data), "\n")
cat("Counties:", n_distinct(gam_data$county), "\n")
cat("Years:", paste(range(gam_data$year, na.rm = TRUE), collapse = " - "), "\n\n")

# GAM with smooth terms for climate variables on monthly mosquito counts
gam_model <- gam(
  total_count ~ 
    s(mean_temp, k = 5) + 
    s(total_precip, k = 5) + 
    s(ndvi, k = 5) +
    s(total_precip_lag1, k = 5) +
    s(total_precip_lag2, k = 5) +
    s(month, bs = "cc", k = min(6, length(unique(gam_data$month)))) +  # Cyclic spline for month
    malaria_burden,
  data = gam_data,
  family = poisson(link = "log"),
  method = "REML"
)

# Model summary
summary(gam_model)

# Plot smooth terms
par(mfrow = c(2, 3))
plot(gam_model, shade = TRUE, pages = 0)
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# 3.2 IDENTIFY OPTIMAL LAGS
# -----------------------------------------------------------------------------

compare_lag_models <- function(data, outcome, climate_var) {
  
  results <- tibble(
    variable = character(),
    lag = integer(),
    AIC = numeric(),
    deviance_explained = numeric()
  )
  
  for (lag_i in 0:3) {
    if (lag_i == 0) {
      var_name <- climate_var
    } else {
      var_name <- paste0(climate_var, "_lag", lag_i)
    }
    
    if (!var_name %in% names(data)) next
    
    model_data <- data %>%
      filter(!is.na(!!sym(var_name)), !is.na(!!sym(outcome)))
    
    if (nrow(model_data) < 20) next
    
    formula_str <- paste0(outcome, " ~ s(", var_name, ", k = 5)")
    
    tryCatch({
      model <- gam(as.formula(formula_str), 
                   data = model_data, 
                   family = poisson(link = "log"))
      
      results <- bind_rows(results, tibble(
        variable = climate_var,
        lag = lag_i,
        AIC = AIC(model),
        deviance_explained = summary(model)$dev.expl * 100
      ))
    }, error = function(e) NULL)
  }
  
  return(results)
}

# Compare lags for each climate variable (on total_count)
lag_results <- bind_rows(
  compare_lag_models(gam_data, "total_count", "mean_temp"),
  compare_lag_models(gam_data, "total_count", "total_precip"),
  compare_lag_models(gam_data, "total_count", "ndvi")
) %>%
  mutate(variable = case_when(
    variable == "mean_temp" ~ "Temperature",
    variable == "total_precip" ~ "Precipitation",
    variable == "ndvi" ~ "NDVI"
  ))

cat("Lag comparison results:\n")
print(lag_results)

# Identify best lag for each climate variable
best_lags <- lag_results %>%
  group_by(variable) %>%
  slice_min(AIC, n = 1) %>%
  ungroup()

cat("\n=== OPTIMAL LAGS FOR ABUNDANCE ===\n")
print(best_lags)

# Visualize lag comparison
p_lag_aic <- lag_results %>%
  ggplot(aes(x = factor(lag), y = AIC, fill = variable)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Temperature" = "#d73027", 
                               "Precipitation" = "#4575b4", 
                               "NDVI" = "#1a9850")) +
  labs(
    title = "Model Fit (AIC) by Lag Period",
    subtitle = "Lower AIC = Better fit",
    x = "Lag (months)",
    y = "AIC",
    fill = "Climate Variable"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_lag_aic)


# =============================================================================
# VARIABLE IMPORTANCE (SPOROZOITE PREDICTION MODELS)
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1 RANDOM FOREST IMPORTANCE
# -----------------------------------------------------------------------------

rf_importance <- rf_model$importance %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(
    Variable_label = case_when(
      Variable == "mean_temp"       ~ "Temperature (current)",
      Variable == "mean_temp_lag1"  ~ "Temperature (1-month lag)",
      Variable == "mean_temp_lag2"  ~ "Temperature (2-month lag)",
      Variable == "mean_temp_lag3"  ~ "Temperature (3-month lag)",
      
      Variable == "total_precip"       ~ "Precipitation (current)",
      Variable == "total_precip_lag1"  ~ "Precipitation (1-month lag)",
      Variable == "total_precip_lag2"  ~ "Precipitation (2-month lag)",
      Variable == "total_precip_lag3"  ~ "Precipitation (3-month lag)",
      
      Variable %in% c("ndvi", "mean_ndvi")        ~ "NDVI (current)",
      Variable %in% c("ndvi_lag1", "mean_ndvi_lag1") ~ "NDVI (1-month lag)",
      Variable %in% c("ndvi_lag2", "mean_ndvi_lag2") ~ "NDVI (2-month lag)",
      Variable %in% c("ndvi_lag3", "mean_ndvi_lag3") ~ "NDVI (3-month lag)",
      
      Variable %in% c("month", "month_num") ~ "Month",
      TRUE ~ Variable
    ),
    Model = "Random Forest"
  )

# -----------------------------------------------------------------------------
# 5.2 XGBOOST IMPORTANCE
# -----------------------------------------------------------------------------

xgb_importance <- xgb.importance(model = xgb_model) %>%
  as.data.frame() %>%
  mutate(
    Variable_label = case_when(
      Feature == "mean_temp"       ~ "Temperature (current)",
      Feature == "mean_temp_lag1"  ~ "Temperature (1-month lag)",
      Feature == "mean_temp_lag2"  ~ "Temperature (2-month lag)",
      Feature == "mean_temp_lag3"  ~ "Temperature (3-month lag)",
      
      Feature == "total_precip"       ~ "Precipitation (current)",
      Feature == "total_precip_lag1"  ~ "Precipitation (1-month lag)",
      Feature == "total_precip_lag2"  ~ "Precipitation (2-month lag)",
      Feature == "total_precip_lag3"  ~ "Precipitation (3-month lag)",
      
      Feature %in% c("ndvi", "mean_ndvi")        ~ "NDVI (current)",
      Feature %in% c("ndvi_lag1", "mean_ndvi_lag1") ~ "NDVI (1-month lag)",
      Feature %in% c("ndvi_lag2", "mean_ndvi_lag2") ~ "NDVI (2-month lag)",
      Feature %in% c("ndvi_lag3", "mean_ndvi_lag3") ~ "NDVI (3-month lag)",
      
      Feature %in% c("month", "month_num") ~ "Month",
      TRUE ~ Feature
    ),
    Model = "XGBoost"
  )

# -----------------------------------------------------------------------------
# 5.3 LIGHTGBM IMPORTANCE
# -----------------------------------------------------------------------------

lgb_importance <- lgb.importance(lgb_model) %>%
  as.data.frame() %>%
  mutate(
    Variable_label = case_when(
      Feature == "mean_temp"       ~ "Temperature (current)",
      Feature == "mean_temp_lag1"  ~ "Temperature (1-month lag)",
      Feature == "mean_temp_lag2"  ~ "Temperature (2-month lag)",
      Feature == "mean_temp_lag3"  ~ "Temperature (3-month lag)",
      
      Feature == "total_precip"       ~ "Precipitation (current)",
      Feature == "total_precip_lag1"  ~ "Precipitation (1-month lag)",
      Feature == "total_precip_lag2"  ~ "Precipitation (2-month lag)",
      Feature == "total_precip_lag3"  ~ "Precipitation (3-month lag)",
      
      Feature %in% c("ndvi", "mean_ndvi")        ~ "NDVI (current)",
      Feature %in% c("ndvi_lag1", "mean_ndvi_lag1") ~ "NDVI (1-month lag)",
      Feature %in% c("ndvi_lag2", "mean_ndvi_lag2") ~ "NDVI (2-month lag)",
      Feature %in% c("ndvi_lag3", "mean_ndvi_lag3") ~ "NDVI (3-month lag)",
      
      Feature %in% c("month", "month_num") ~ "Month",
      TRUE ~ Feature
    ),
    Model = "LightGBM"
  )

# -----------------------------------------------------------------------------
# 5.4 PLOT VARIABLE IMPORTANCE
# -----------------------------------------------------------------------------

# Random Forest importance plot
rf_importance_plot <- rf_importance %>%
  mutate(
    Importance_pct = 100 * MeanDecreaseGini / sum(MeanDecreaseGini),
    Variable_label = fct_reorder(Variable_label, Importance_pct)
  )

# Plot
p_imp_rf <- ggplot(
  rf_importance_plot,
  aes(x = Importance_pct, y = Variable_label)
) +
  geom_col(fill = "#2166ac", alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(Importance_pct, 1), "%")),
    hjust = -0.2,
    size = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.25)),
    labels = label_number(suffix = "%")
  ) +
  labs(
    title = "Random Forest",
    x = "Relative Importance (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

p_imp_rf

# XGBoost importance plot
xgb_importance_plot <- xgb_importance %>%
  mutate(
    Gain_pct = 100 * Gain / sum(Gain),
    Variable_label = fct_reorder(Variable_label, Gain_pct)
  )

# Plot
p_imp_xgb <- ggplot(
  xgb_importance_plot,
  aes(x = Gain_pct, y = Variable_label)
) +
  geom_col(fill = "#d73027", alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(Gain_pct, 1), "%")),
    hjust = -0.2,
    size = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.25)),
    labels = label_number(suffix = "%")
  ) +
  labs(
    title = "XGBoost",
    x = "Relative Gain (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

p_imp_xgb


# LightGBM importance plot
# Prepare LightGBM importance with relative gain
lgb_importance_plot <- lgb_importance %>%
  mutate(
    Gain_pct = 100 * Gain / sum(Gain),
    Variable_label = fct_reorder(Variable_label, Gain_pct)
  )

# Plot
p_imp_lgb <- ggplot(
  lgb_importance_plot,
  aes(x = Gain_pct, y = Variable_label)
) +
  geom_col(fill = "#1a9850", alpha = 0.85) +
  geom_text(
    aes(label = paste0(round(Gain_pct, 1), "%")),
    hjust = -0.2,
    size = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.25)),
    labels = label_number(suffix = "%")
  ) +
  labs(
    title = "LightGBM",
    x = "Relative Gain (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 15,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10)
  )
p_imp_lgb

# Combine
(p_imp_rf + p_imp_xgb + p_imp_lgb) +
  plot_annotation(
    title = "Variable Importance Comparison Across Models",
    subtitle = "Key climate drivers of sporozoite positivity",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  )



# =============================================================================
# COUNTY-LEVEL PREDICTIONS FOR MAPPING (SPOROZOITE RISK)
# =============================================================================

# -----------------------------------------------------------------------------
# 6.1 SELECT BEST MODEL
# -----------------------------------------------------------------------------

# Find best model based on F1 score
best_model_name <- all_metrics %>%
  slice_max(F1_Score, n = 1, with_ties = FALSE) %>%
  pull(Model)

cat("\n=== BEST MODEL:", best_model_name, "===\n")

# -----------------------------------------------------------------------------
# 6.2 GENERATE PREDICTIONS FOR ALL INDIVIDUAL RECORDS
# -----------------------------------------------------------------------------

# Prepare full dataset for prediction 
full_X <- spz_ml_data %>%
  dplyr::select(all_of(available_features)) %>%
  drop_na()

# Corresponding full data frame (same rows, same order)
full_data <- spz_ml_data %>%
  filter(complete.cases(dplyr::select(., all_of(available_features))))

# Generate predictions using best model
if (best_model_name == "Random Forest") {
  pred_prob <- predict(rf_model, full_X, type = "prob")[, "Positive"]
  pred_class <- predict(rf_model, full_X)
} else if (best_model_name == "XGBoost") {
  pred_prob <- predict(xgb_model, xgb.DMatrix(as.matrix(full_X)))
  pred_class <- factor(ifelse(pred_prob > 0.5, "Positive", "Negative"),
                       levels = c("Negative", "Positive"))
} else {
  pred_prob <- predict(lgb_model, as.matrix(full_X))
  pred_class <- factor(ifelse(pred_prob > 0.5, "Positive", "Negative"),
                       levels = c("Negative", "Positive"))
}

full_data <- full_data %>%
  mutate(
    pred_prob = pred_prob,
    pred_class = pred_class
  )

# -----------------------------------------------------------------------------
# 6.3 AGGREGATE BY COUNTY
# -----------------------------------------------------------------------------

county_predictions <- full_data %>%
  group_by(county, malaria_burden) %>%
  dplyr::summarise(
    n_observations = n(),
    
    # Observed sporozoite metrics
    n_positive_spz = sum(spz_outcome == "Positive", na.rm = TRUE),
    spz_rate = ifelse(n_observations > 0, n_positive_spz / n_observations * 100, NA_real_),
    
    # Model-based prediction metrics
    mean_pred_prob = mean(pred_prob, na.rm = TRUE),
    pct_pred_positive = mean(pred_class == "Positive") * 100,
    
    # Climate means
    mean_temp = mean(mean_temp, na.rm = TRUE),
    mean_precip = mean(total_precip, na.rm = TRUE),
    mean_ndvi = mean(ndvi, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Convert to 0–100 scales
    prob_score = mean_pred_prob * 100,
    spz_score = spz_rate,  # Already in percentage
    
    # Composite infection risk score (0–100)
    # Equal weight to observed spz_rate and predicted probability
    risk_score = 0.5 * prob_score + 0.5 * spz_score,
    
    # Risk categories with ranges
    risk_category = case_when(
      risk_score >= 75 ~ "Very High (75–100)",
      risk_score >= 50 ~ "High (50–74)",
      risk_score >= 25 ~ "Moderate (25–49)",
      TRUE ~ "Low (0–24)"
    ),
    risk_category = factor(risk_category,
                           levels = c("Low (0–24)", "Moderate (25–49)",
                                      "High (50–74)", "Very High (75–100)")),
    
    # Sporozoite rate categories
    spz_category = case_when(
      is.na(spz_rate) ~ "No Data",
      spz_rate >= 5   ~ "High (≥5%)",
      spz_rate >= 2   ~ "Moderate (2–5%)",
      spz_rate >= 0.5 ~ "Low (0.5–2%)",
      TRUE            ~ "Very Low (<0.5%)"
    ),
    spz_category = factor(spz_category,
                          levels = c("Very Low (<0.5%)", "Low (0.5–2%)",
                                     "Moderate (2–5%)", "High (≥5%)", "No Data"))
  ) %>%
  arrange(desc(risk_score))

cat("\n=== COUNTY-LEVEL SPOROZOITE RISK PREDICTIONS ===\n")
print(county_predictions, n = nrow(county_predictions))

# =============================================================================
# PUBLICATION-QUALITY MAPS (SPOROZOITE RISK)
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1 JOIN PREDICTIONS TO SHAPEFILE
# -----------------------------------------------------------------------------

kenya_counties <- kenya_counties %>%
  mutate(
    NAME_1 = case_when(
      NAME_1 == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      NAME_1 == "Tharaka-Nithi" ~ "Tharaka Nithi",
      TRUE ~ NAME_1
    )
  ) 

# Ensure county names match; rename NAME_1 to county
kenya_map_data <- kenya_counties %>%
  mutate(county = NAME_1) %>%        # overwrite or create county column
  left_join(county_predictions, by = "county")


# Check for unmatched counties
unmatched <- county_predictions %>%
  filter(!county %in% kenya_map_data$county)

if (nrow(unmatched) > 0) {
  cat("Warning: Unmatched counties in predictions:\n")
  print(unmatched$county)
}

# -----------------------------------------------------------------------------
# 7.2 DEFINE COLOR PALETTES
# -----------------------------------------------------------------------------

# Risk score palette (continuous)
risk_colors_continuous <- c("#1a9850", "#91cf60", "#d9ef8b",
                            "#fee08b", "#fc8d59", "#d73027", "#a50026")

# Risk category palette (discrete with ranges)
risk_colors_discrete <- c(
  "Low (0–24)"         = "#1a9850",
  "Moderate (25–49)"   = "#fee08b",
  "High (50–74)"       = "#fc8d59",
  "Very High (75–100)" = "#d73027"
)

# Sporozoite rate palette
spz_colors <- c(
  "Very Low (<0.5%)" = "#1a9850",
  "Low (0.5–2%)"     = "#a6d96a",
  "Moderate (2–5%)"  = "#fdae61",
  "High (≥5%)"       = "#d73027",
  "No Data"          = "#cccccc"
)

# -----------------------------------------------------------------------------
# 7.3 MAP 1: INFECTION RISK SCORE (CONTINUOUS)
# -----------------------------------------------------------------------------

p_map_risk_continuous <- ggplot(kenya_map_data) +
  geom_sf(aes(fill = risk_score), color = "white", size = 0.3) +
  scale_fill_gradientn(
    colors = risk_colors_continuous,
    values = scales::rescale(c(0, 17.5, 35, 52.5, 70)),  # scale stops for gradient
    limits = c(0, 70),                                    # max at 70
    breaks = c(0, 17.5, 35, 52.5, 70),
    labels = c("0", "25", "50", "75", "100"),            # optional: keep original labels
    na.value = "gray90",
    name = "Infection\nRisk Score"
  ) +
  labs(
    title = "Malaria Vector Infection Risk Score by County",
    subtitle = "Composite of observed sporozoite rate and predicted infectivity"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

print(p_map_risk_continuous)


# -----------------------------------------------------------------------------
# 7.4 MAP 2: INFECTION RISK CATEGORIES (DISCRETE)
# -----------------------------------------------------------------------------

p_map_risk_category <- ggplot(kenya_map_data) +
  geom_sf(aes(fill = risk_category), color = "white", size = 0.3) +
  scale_fill_manual(
    values = risk_colors_discrete,
    na.value = "gray90",
    name = "Infection\nRisk Category"
  ) +
  labs(
    title = "Malaria Vector Infection Risk Categories by County",
    subtitle = "Based on composite infection risk score"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

print(p_map_risk_category)

# -----------------------------------------------------------------------------
# 7.5 MAP 3: OBSERVED SPOROZOITE RATE
# -----------------------------------------------------------------------------

p_map_spz <- ggplot(kenya_map_data) +
  geom_sf(aes(fill = spz_category), color = "white", size = 0.3) +
  scale_fill_manual(
    values = spz_colors,
    na.value = "gray90",
    name = "Observed\nSporozoite Rate"
  ) +
  labs(
    title = "Observed Sporozoite Infection Rate by County",
    subtitle = "Proportion of mosquitoes testing positive"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

print(p_map_spz)

# -----------------------------------------------------------------------------
# 7.6 MAP 4: PREDICTED INFECTIVITY PROBABILITY
# -----------------------------------------------------------------------------

# Create probability bins
kenya_map_data <- kenya_map_data %>%
  mutate(
    prob_bins = cut(mean_pred_prob * 100,
                    breaks = c(0, 20, 40, 60, 80, 100),
                    labels = c("0–20%", "21–40%", "41–60%", "61–80%", "81–100%"),
                    include.lowest = TRUE)
  )

p_map_prob <- ggplot(kenya_map_data) +
  geom_sf(aes(fill = prob_bins), color = "black", size = 0.3) +
  scale_fill_brewer(palette = "YlOrRd", na.value = "gray90", 
                    name = "Predicted\nPositivity") +
  labs(
    title = "Predicted Probability of Sporozoite Positivity by County",
    subtitle = paste("Based on", best_model_name, "model")
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

print(p_map_prob)


# =============================================================================
# RISK MAPPING AND PREDICTIONS
# =============================================================================
# Read climate data
ndvi <- read_csv("climate_data.csv")
temp_precip <- read_csv("temp_precip2.csv")

ndvi_clean <- ndvi %>%
  dplyr::select(county, year, month, ndvi) 

temp_precip_clean <- temp_precip %>%
  dplyr::select(county, year, month, mean_temp, total_precip) 

ndvi_clean <- ndvi_clean %>%
  mutate(
    county = case_when(
      county == "Homabay" ~ "Homa Bay",
      county == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      county == "Tharaka-Nithi" ~ "Tharaka Nithi",
      county == "Trans-Nzoia" ~ "Trans Nzoia",
      county == "Uasin-Gishu" ~ "Uasin Gishu",
      county == "West-Pokot" ~ "West Pokot",
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Tana-River" ~ "Tana River",
      county == "Homa-Bay" ~ "Homa Bay",
      TRUE ~ county
    )
  ) 

temp_precip_clean <- temp_precip_clean %>%
  mutate(
    county = case_when(
      county == "Homabay" ~ "Homa Bay",
      county == "Elgeyo-Marakwet" ~ "Elgeyo Marakwet",
      county == "Tharaka-Nithi" ~ "Tharaka Nithi",
      county == "Trans-Nzoia" ~ "Trans Nzoia",
      county == "Uasin-Gishu" ~ "Uasin Gishu",
      county == "West-Pokot" ~ "West Pokot",
      county == "Taita-Taveta" ~ "Taita Taveta",
      county == "Tana-River" ~ "Tana River",
      county == "Homa-Bay" ~ "Homa Bay",
      TRUE ~ county
    )
  ) 

climate_data <- temp_precip_clean %>%
  left_join(ndvi_clean, by = c("county", "year", "month"))

# Filter to years matching entomological data (2016-2024)
climate_data <- climate_data %>%
  filter(year >= 2016 & year <= 2024)

# First, create lagged variables from the CLIMATE data (which has continuous months)
climate_with_lags <- climate_data %>%
  arrange(county, year, month) %>%
  group_by(county) %>%
  mutate(
    # Temperature lags
    mean_temp_lag1 = lag(mean_temp, 1),
    mean_temp_lag2 = lag(mean_temp, 2),
    mean_temp_lag3 = lag(mean_temp, 3),
    
    # Precipitation lags
    total_precip_lag1 = lag(total_precip, 1),
    total_precip_lag2 = lag(total_precip, 2),
    total_precip_lag3 = lag(total_precip, 3),
    
    # NDVI lags
    mean_ndvi_lag1 = lag(ndvi, 1),
    mean_ndvi_lag2 = lag(ndvi, 2),
    mean_ndvi_lag3 = lag(ndvi, 3)
  ) %>%
  ungroup()


# Aggregate anopheles data by county, year, month
monthly_vectors <- anopheles_final %>%
  filter(!is.na(year), !is.na(month_num), !is.na(county)) %>%
  group_by(county, year, month = month_num) %>%
  summarise(
    # Abundance
    total_abundance = n(),
    n_gambiae = sum(species == "An. gambiae s.l.", na.rm = TRUE),
    n_funestus = sum(species == "An. funestus s.l.", na.rm = TRUE),
    
    # Subspecies counts
    n_arabiensis = sum(subspecies == "An. arabiensis", na.rm = TRUE),
    n_gambiae_ss = sum(subspecies == "An. gambiae s.s.", na.rm = TRUE),
    n_funestus_ss = sum(subspecies == "An. funestus s.s.", na.rm = TRUE),
    
    # Sporozoite
    n_tested_spz = sum(!is.na(sporozoite_positive)),
    n_positive_spz = sum(sporozoite_positive, na.rm = TRUE),
    spz_rate = ifelse(n_tested_spz > 0, n_positive_spz / n_tested_spz * 100, NA),
    
    # Blood meal / HBI
    n_bloodmeal = sum(!is.na(human_bloodmeal)),
    n_human_blood = sum(human_bloodmeal, na.rm = TRUE),
    hbi = ifelse(n_bloodmeal > 0, n_human_blood / n_bloodmeal * 100, NA),
    
    # Malaria burden (take first value since it's constant per county)
    malaria_burden = first(malaria_burden),
    
    .groups = "drop"
  )


# Now merge the monthly vectors with climate data that has lags
vector_climate <- monthly_vectors %>%
  left_join(climate_with_lags, by = c("county", "year", "month"))

# Add season
vector_climate <- vector_climate %>%
  mutate(
    season = case_when(
      month %in% c(1, 2) ~ "Short Dry (Jan-Feb)",
      month %in% c(3, 4, 5) ~ "Long Rains (Mar-May)",
      month %in% c(6, 7, 8, 9) ~ "Long Dry (Jun-Sep)",
      month %in% c(10, 11, 12) ~ "Short Rains (Oct-Dec)"
    ),
    season = factor(season, levels = c(
      "Short Dry (Jan-Feb)",
      "Long Rains (Mar-May)",
      "Long Dry (Jun-Sep)",
      "Short Rains (Oct-Dec)"
    )),
    month_name = month.abb[month]
  )



# -----------------------------------------------------------------------------
# 4.1 CREATE RISK SCORES BY COUNTY
# -----------------------------------------------------------------------------

county_risk <- vector_climate %>%
  group_by(county, malaria_burden) %>%
  dplyr::summarise(
    mean_abundance = mean(total_abundance),
    total_abundance = sum(total_abundance),
    mean_spz_rate = mean(spz_rate, na.rm = TRUE),
    mean_hbi = mean(hbi, na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Normalize scores
    abundance_score = scales::rescale(mean_abundance, to = c(0, 100)),
    spz_score = scales::rescale(replace_na(mean_spz_rate, 0), to = c(0, 100)),
    hbi_score = scales::rescale(replace_na(mean_hbi, 0), to = c(0, 100)),
    
    # Composite risk score
    risk_score = 0.4 * abundance_score + 0.4 * spz_score + 0.2 * hbi_score,
    
    # Risk category
    risk_category = case_when(
      risk_score >= 75 ~ "Very High",
      risk_score >= 50 ~ "High",
      risk_score >= 25 ~ "Moderate",
      TRUE ~ "Low"
    ),
    risk_category = factor(risk_category, 
                           levels = c("Low", "Moderate", "High", "Very High"))
  ) %>%
  arrange(desc(risk_score))

cat("\n=== COUNTY RISK SCORES ===\n")
print(county_risk, n = 20)

# Risk score visualization
p_risk <- county_risk %>%
  mutate(county = fct_reorder(county, risk_score)) %>%
  ggplot(aes(x = risk_score, y = county, fill = risk_category)) +
  geom_col() +
  geom_text(aes(label = round(risk_score, 1)), hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("Low" = "#1a9850", "Moderate" = "#fee08b",
                               "High" = "#fc8d59", "Very High" = "#d73027")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Malaria Vector Risk Score by County",
    x = "Risk Score (0-100)",
    y = NULL,
    fill = "Risk Category"
  )

print(p_risk)


# -----------------------------------------------------------------------------
# 4.2 SEASONAL RISK BY MALARIA BURDEN ZONE
# -----------------------------------------------------------------------------
seasonal_zone_risk <- vector_climate %>%
  group_by(malaria_burden, season) %>%
  dplyr::summarise(
    mean_abundance = mean(total_abundance),
    mean_spz_rate = mean(spz_rate, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_zone_season <- seasonal_zone_risk %>%
  ggplot(aes(x = season, y = mean_abundance, fill = malaria_burden)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = zone_colors) +
  labs(
    title = "Seasonal Vector Abundance by Malaria Burden Zone",
    x = NULL,
    y = "Mean Abundance",
    fill = "Malaria Burden"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_zone_season)



# -----------------------------------------------------------------------------
# 4.3 CONTROL RECOMMENDATIONS
# -----------------------------------------------------------------------------
control_recommendations <- county_risk %>%
  mutate(
    # IRS priority
    irs_priority = case_when(
      risk_category == "Very High" ~ "Immediate - Full coverage",
      risk_category == "High" ~ "High - Targeted coverage",
      risk_category == "Moderate" ~ "Medium - Focal spraying",
      TRUE ~ "Low - Reactive only"
    ),
    
    # LLIN priority
    llin_priority = case_when(
      risk_category %in% c("Very High", "High") ~ "Universal coverage",
      risk_category == "Moderate" ~ "Targeted distribution",
      TRUE ~ "Routine replacement"
    ),
    
    # Best timing based on seasonal patterns
    timing = case_when(
      malaria_burden %in% c("Very High Burden", "High Burden") ~ 
        "Year-round with intensification before long rains (Feb-Mar)",
      malaria_burden == "Moderate-High Burden" ~ 
        "Seasonal - Pre-rainy season (Feb & Sep)",
      TRUE ~ "Reactive - Based on surveillance"
    )
  ) %>%
  dplyr::select(county, malaria_burden, risk_category, risk_score,
                irs_priority, llin_priority, timing)

cat("\n=== CONTROL RECOMMENDATIONS ===\n")
print(control_recommendations, n = 20)

# Define colors for IRS priority
irs_colors <- c(
  "Immediate - Full coverage" = "#d73027",
  "High - Targeted coverage" = "#fc8d59",
  "Medium - Focal spraying" = "#fee08b",
  "Low - Reactive only" = "#1a9850"
)

p_irs <- control_recommendations %>%
  mutate(
    county = fct_reorder(county, risk_score),
    irs_priority = factor(irs_priority, levels = c(
      "Immediate - Full coverage",
      "High - Targeted coverage",
      "Medium - Focal spraying",
      "Low - Reactive only"
    ))
  ) %>%
  ggplot(aes(x = risk_score, y = county, fill = irs_priority)) +
  geom_col() +
  geom_text(aes(label = round(risk_score, 1)), hjust = -0.2, size = 3) +
  scale_fill_manual(values = irs_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Indoor Residual Spraying (IRS) Priority by County",
    subtitle = "Based on composite risk score",
    x = "Risk Score",
    y = NULL,
    fill = "IRS Priority"
  ) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 11, hjust = 0.5))

print(p_irs)

# Define colors for LLIN priority
llin_colors <- c(
  "Universal coverage" = "#d73027",
  "Targeted distribution" = "#fee08b",
  "Routine replacement" = "#1a9850"
)

p_llin <- control_recommendations %>%
  mutate(
    county = fct_reorder(county, risk_score),
    llin_priority = factor(llin_priority, levels = c(
      "Universal coverage",
      "Targeted distribution",
      "Routine replacement"
    ))
  ) %>%
  ggplot(aes(x = risk_score, y = county, fill = llin_priority)) +
  geom_col() +
  geom_text(aes(label = round(risk_score, 1)), hjust = -0.2, size = 3) +
  scale_fill_manual(values = llin_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Long-Lasting Insecticidal Net (LLIN) Priority by County",
    subtitle = "Based on composite risk score",
    x = "Risk Score",
    y = NULL,
    fill = "LLIN Priority"
  ) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))

print(p_llin)

# Summary by priority levels
priority_summary <- control_recommendations %>%
  dplyr::count(irs_priority) %>%
  mutate(
    irs_priority = factor(irs_priority, levels = c(
      "Immediate - Full coverage",
      "High - Targeted coverage",
      "Medium - Focal spraying",
      "Low - Reactive only"
    ))
  )

p_priority_summary <- priority_summary %>%
  ggplot(aes(x = irs_priority, y = n, fill = irs_priority)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = irs_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Number of Counties by IRS Priority Level",
    x = NULL,
    y = "Number of Counties"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_priority_summary)


# Risk category colors
risk_colors <- c(
  "Very High" = "#d73027",
  "High" = "#fc8d59",
  "Moderate" = "#fee08b",
  "Low" = "#1a9850"
)

p_risk_dist <- control_recommendations %>%
  mutate(risk_category = factor(risk_category, levels = c("Low", "Moderate", "High", "Very High"))) %>%
  dplyr::count(risk_category) %>%
  ggplot(aes(x = risk_category, y = n, fill = risk_category)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = risk_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Distribution of Counties by Risk Category",
    x = "Risk Category",
    y = "Number of Counties"
  )

print(p_risk_dist)


library(patchwork)
# Combine all plots
(p_irs + p_llin) /
  (p_priority_summary + p_risk_dist) +
  plot_annotation(
    title = "Vector Control Recommendations Dashboard",
    subtitle = "Based on integrated risk assessment (2016-2024)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

# Timing Recommendations
# Timing colors
timing_colors <- c(
  "Year-round with intensification before long rains (Feb-Mar)" = "#fc8d59",
  "Seasonal - Pre-rainy season (Feb & Sep)" = "#d73027",
  "Reactive - Based on surveillance" = "#1a9850"
)

p_timing <- control_recommendations %>%
  dplyr::count(timing) %>%
  mutate(timing = str_wrap(timing, width = 30)) %>%
  ggplot(aes(x = reorder(timing, n), y = n, fill = timing)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("#1a9850", "#d73027", "#fc8d59")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Intervention Timing Recommendations",
    subtitle = "Number of counties per timing strategy",
    x = NULL,
    y = "Number of Counties"
  )+
  theme(
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_timing)





##########################################
# library(caret)
# 
df_analysis <- ml_data %>%
  mutate(
    # Convert categorical variables to factors
    county = as.factor(county),
    season = as.factor(season),
    malaria_burden = as.factor(malaria_burden),
    # Ensure month is a factor for seasonality checks
    # month_f = as.factor(month)
  )
df_analysis <- vector_climate %>%
  mutate(
    # Convert categorical variables to factors
    county = as.factor(county),
    season = as.factor(season),
    malaria_burden = as.factor(malaria_burden),
    # Ensure month is a factor for seasonality checks
    month_f = as.factor(month)
  )


# Sporozoite Rate & HBI Summary ---
infectivity_summary <- df_analysis %>%
  group_by(county) %>%
  summarise(
    Total_Mosquitoes = sum(total_abundance),
    Avg_Spz_Rate = mean(spz_rate, na.rm = TRUE),
    Avg_HBI = mean(hbi, na.rm = TRUE),
    Total_Positive = sum(n_positive_spz)
  )

print(infectivity_summary)


# Summarize Data by County
map_data <- vector_climate %>%
  group_by(county) %>%
  dplyr::summarise(
    Total_Abundance = sum(total_abundance, na.rm = TRUE),
    Mean_Spz_Rate = mean(spz_rate, na.rm = TRUE)
  )

# Merge Data with Shapefile
kenya_map_merged <- kenya_counties %>%
  left_join(map_data, by = c("NAME_1" = "county"))

# Map 1: Vector Abundance
map1<-ggplot(data = kenya_map_merged) +
  geom_sf(aes(fill = Total_Abundance)) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       na.value = "grey90", name = "Abundance") +
  theme_minimal() +
  labs(title = "Vector Distribution: Total Abundance (2016-2024)",
       subtitle = "Grey regions indicate no sampling data") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"))
map1

# Map 2: Sporozoite Rate (Infectivity)
map2<-ggplot(data = kenya_map_merged) +
  geom_sf(aes(fill = Mean_Spz_Rate)) +
  scale_fill_viridis_c(option = "inferno",
                       na.value = "grey90",direction = -1, name = "Spz Rate") +
  theme_minimal() +
  labs(title = "Malaria Risk: Mean Sporozoite Rate",
       subtitle = "Grey regions indicate no sampling data") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"))
map2




