README

Description of the Dataset
--------------------------
This folder contains meteorological temperature data for the months June to September 2019 from the weather staion Aarhus South.
The data are stored as CSV files, with one file per month.

File Naming
----------------------
The files follow the naming scheme:

    YYYY_MM_temp.csv

- YYYY = year (here: 2019)
- MM = month (06 = June, 07 = July, 08 = August, 09 = September)
- temp = temperature

File Overview
----------------------
| File             | Size (approx.) | Parameter   | Month          |
|------------------|----------------|-------------|----------------|
| 2019_06_temp.csv | 26 KB          | Temperature | June 2019      |
| 2019_07_temp.csv | 27 KB          | Temperature | July 2019      |
| 2019_08_temp.csv | 27 KB          | Temperature | August 2019    |
| 2019_09_temp.csv | 26 KB          | Temperature | September 2019 |

Structure of the CSV Files
--------------------------
Each CSV file contains a time series of meteorological observations. The column names are in Danish (as provided by DMI).

Examples:
---------------------
| DateTime           | Middeltemperatur | Maksimumtemperatur | Minimumtemperatur |
|--------------------|------------------|--------------------|-------------------|
| 01.06.2019 00:00   | 13.1             | 13.2               | 12.9              |
| 01.06.2019 01:00   | 13.1             | 13.2               | 12.9              |

- "Middeltemperatur": mean temperature (°C)
- "Maksimumtemperatur": maximum temperature (°C)
- "Minimumtemperatur": minimum temperature (°C)

Data Source
-----------
The data were obtained from the Danish Meteorological Institute (DMI) via the official free data portal:

    https://www.dmi.dk/friedata/observationer/

- Station: Aarhus South
- Download date: 26 August 2023

License and Citation
--------------------
The dataset is distributed under the open data policy of the Danish Meteorological Institute (DMI).
See: https://www.dmi.dk/friedata/

If these data are used in a scientific publication, the source should be acknowledged as:

    Danish Meteorological Institute (DMI), Aarhus South station, data retrieved from the DMI free data portal (https://www.dmi.dk/friedata/observationer/), downloaded on 26 August 2023.
