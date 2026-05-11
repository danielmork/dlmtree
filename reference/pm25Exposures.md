# PM2.5 Exposure data

Data.frame containing a sample of weekly average PM2.5 exposures across
a range of states/counties. The PM2.5 data was downloaded from US EPA
(https://aqs.epa.gov/aqsweb/airdata/download_files.html) daily data
summaries and averaged by week. Forty-week ranges were assess for
non-missingness and grouped for this dataset.

## Usage

``` r
data(pm25Exposures)
```

## Format

data.frame; columns: S = state, C = city, 1-40 = weekly exposure data

## Source

<https://aqs.epa.gov/aqsweb/airdata/download_files.html>

## References

<https://www.epa.gov/outdoor-air-quality-data>
