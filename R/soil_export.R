library(writexl)
library(dplyr)


write_soil_data <- function(combined_soil_data, combined_soil_profiles) {
  writexl::write_xlsx(list(soil = combined_soil_data, soil_profile = combined_soil_profiles),
                      path = "./combined_soil_data.xlsx")
  saveRDS(list(soil = combined_soil_data, soil_profile = combined_soil_profiles), "soils.RDS")
}
