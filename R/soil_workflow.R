r_base_directory <- file.path(getwd(), "R")
source(file.path(r_base_directory, "soil_helpers.R"))
source(file.path(r_base_directory, "soil_fetch.R"))
source(file.path(r_base_directory, "soil_processing.R"))
source(file.path(r_base_directory, "soil_export.R"))

library(dplyr)


# Main workflow
lonlat_tibble <- read.csv("coord_to_SSURGO_import.csv")

# Initialize lists to store results
all_soil_data <- list()
all_soil_profiles <- list()
failed_soil_fetch <- c()

# Step 2: Fetch soil profiles
for (i in 1:nrow(lonlat_tibble)) {
  print(paste("Processing row:", i))
  lonlat <- c(lonlat_tibble$longitude[i], lonlat_tibble$latitude[i])
  environment_id <- lonlat_tibble$environment_id[i]
  
  # Try to fetch soil profile
  ssurgo_profiles <- tryCatch({
    get_ssurgo_soil_profile(environment_id, lonlat, nsoil = 1)
  }, error = function(e) {
    message(paste("Failed to fetch soil profile for environment", environment_id, ":", conditionMessage(e)))
    failed_soil_fetch <<- c(failed_soil_fetch, environment_id)
    return(NULL)
  })
  
  if (is.null(ssurgo_profiles)) next
  
  # Process profiles for each soil
  for (soil_index in 1:length(ssurgo_profiles)) {
    soil_names <- names(attributes(ssurgo_profiles[[soil_index]]))
    soil <- attributes(ssurgo_profiles[[soil_index]])[!soil_names %in% names(attributes(data.frame()))]
    soil <- as.data.frame(lapply(soil, function(x) if (length(x) == 0) NA else x))
    
    # Add environment ID to soil data
    soil$environment_id <- environment_id
    
    profile <- ssurgo_profiles[[soil_index]]
    profile$environment_id <- environment_id
    
    # Store results
    all_soil_data[[length(all_soil_data) + 1]] <- soil
    all_soil_profiles[[length(all_soil_profiles) + 1]] <- profile
  }
}

# Step 3: Combine results
combined_soil_data <- bind_rows(all_soil_data)
combined_soil_profiles <- bind_rows(all_soil_profiles)

# Step 4: Save results
# Save combined soil data and profiles to an Excel file
write_soil_data(combined_soil_data, combined_soil_profiles)

# Write missing environment IDs to a text file if any fetches failed
if (length(failed_soil_fetch) > 0) {
  write.table(failed_soil_fetch, file = "failed_soil_fetch.txt", row.names = FALSE, col.names = FALSE)
  message("Failed soil fetch environment IDs written to 'failed_soil_fetch.txt'")
} else {
  message("No failed soil fetch environment IDs to write.")
}

# Save combined data as RDS for further analysis
saveRDS(list(soil_data = combined_soil_data, soil_profiles = combined_soil_profiles), file = "combined_soils.RDS")

# Done!
message("Workflow completed successfully!")
