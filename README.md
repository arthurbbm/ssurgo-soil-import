# Soil Data Workflow

This repository provides a modular workflow to fetch, process, and analyze soil data using the SSURGO database. It accepts a CSV file containing environment information, geographical coordinates, and metadata, processes the soil profiles, and outputs the results in multiple formats.

------------------------------------------------------------------------

## **Repository Structure**

-   **`soil_helpers.R`**: Contains utility functions used throughout the workflow.
-   **`soil_fetch.R`**: Handles data fetching from the SSURGO database.
-   **`soil_processing.R`**: Processes and transforms raw soil data.
-   **`soil_export.R`**: Saves the processed data to files.
-   **`soil_workflow.R`**: Main workflow script orchestrating the entire process.

------------------------------------------------------------------------

## **Input File**

The input file must be a CSV file with the following columns:

| year | environment_id | latitude | longitude | plantingDate | harvestDate | state | source |
|---------|---------|---------|---------|---------|---------|---------|---------|
| 2024 | Env1 | 40.123456 | -94.123456 | 2024-05-15 | 2024-09-20 | IA | Dataset A |
| 2024 | Env2 | 40.234567 | -94.234567 | 2024-05-20 | 2024-09-25 | MO | Dataset B |

### **Columns Description**

1.  **`year`**: Year of the data collection or observation (e.g., `2024`).
2.  **`environment_id`**: A unique identifier for the environment or location (e.g., `Env1`).
3.  **`latitude`**: Latitude coordinate in decimal degrees (e.g., `40.123456`).
4.  **`longitude`**: Longitude coordinate in decimal degrees (e.g., `-94.123456`).
5.  **`plantingDate`**: The planting date for the crop in `YYYY-MM-DD` format.
6.  **`harvestDate`**: The harvest date for the crop in `YYYY-MM-DD` format.
7.  **`state`**: The U.S. state where the environment is located (e.g., `IA` for Iowa).
8.  **`source`**: The dataset or source of the data (e.g., `Dataset A`).

**Example Input File**: `coord_to_SSURGO_import.csv`

Ensure the file has the exact column names and structure as shown above.

------------------------------------------------------------------------

## **Setup Instructions**

### **1. Clone the Repository**

Clone this repository to your local machine:

``` bash
git clone <repository-url>
cd <repository-folder>
```

### **2. Install Required R Packages**

Install the required R packages using the following script:

``` r
required_packages <- c(
  "stringr", "soiltexture", "assertthat", "soilDB", "aqp",
  "spData", "sf", "writexl", "dplyr", "httr", "xml2", "apsimx"
)

to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install)
}
```

### **3. Prepare Input Data**

Ensure your input file is named `coord_to_SSURGO_import.csv` and placed in the same directory as the R scripts.

------------------------------------------------------------------------

## **Running the Workflow**

### **1. Execute the Workflow**

Run the `soil_workflow.R` script in R:

``` r
source("soil_workflow.R")
```

### **2. Outputs**

The workflow generates the following files in the working directory:

1.  **`combined_soil_data.xlsx`**: An Excel file with combined soil data and soil profiles.
2.  **`combined_soils.RDS`**: A serialized R object for further analysis in R.
3.  **`failed_soil_fetch.txt`**: A text file listing `environment_id`s where soil profiles could not be fetched.

------------------------------------------------------------------------

## **Troubleshooting**

1.  **Dependency Errors**: Ensure all required packages are installed. If errors persist, rerun the package installation script above.

2.  **Input File Issues**: Verify that the CSV file has the correct structure and valid coordinates.

3.  **Failed Fetches**: Check `failed_soil_fetch.txt` for a list of `environment_id`s where soil profiles could not be retrieved.

------------------------------------------------------------------------

## **Contact**

If you encounter any issues or have questions, feel free to contact the repository maintainer.
