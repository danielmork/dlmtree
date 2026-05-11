# Download simulated data for dlmtree articles

Download simulated data for dlmtree articles

## Usage

``` r
get_sbd_dlmtree()
```

## Value

A data frame with 10000 rows (observations) and 202 variables. All data
is simulated. The variables are:

- bwgaz:

  Outcome to be used. Simulated birth weight for gestational age
  z-score.

- ChildSex:

  Binary sex of child.

- MomAge:

  Continuous age in years.

- GestAge:

  Continuous estimated gestational age at birth in weeks.

- MomHeightIn:

  Continuous maternal height in inches.

- MomPriorWeightLbs:

  Continuous mothers pre-pregnancy weight in pounds.

- MomPriorBMI:

  Continuous mothers pre-pregnancy BMI.

- race:

  Categorical race.

- Hispanic:

  Binary indicator of Hispanic.

- MomEdu:

  Categorical maternal highest educational attainment.

- SmkAny:

  Binary indicator of any smoking during pregnancy.

- Marital:

  Categorical maternal marital status.

- Income:

  Categorical income.

- EstDateConcept:

  Estimated date of conception.

- EstMonthConcept:

  Estimated month of conception.

- EstYearConcept:

  Estimated year of conception.

- pm25_1 - pm25_37:

  Weekly average exposure to PM2.5 for weeks 1 to 37. The columns are
  already scaled by the exposure IQR of 0.35.

- no2_1 - no2_37:

  Weekly average exposure to NO2 for weeks 1 to 37. The columns are
  already scaled by the exposure IQR of 9.13.

- so2_1 - so2_37:

  Weekly average exposure to SO2 for weeks 1 to 37. The columns are
  already scaled by the exposure IQR of 0.96.

- co2_1 - co2_37:

  Weekly average exposure to CO for weeks 1 to 37. The columns are
  already scaled by the exposure IQR of 0.15.

- temp_1 - temp_37:

  Weekly average exposure to temperature for weeks 1 to 37. The columns
  are already scaled by the exposure IQR of 27.93

- source:

  Variable indicating that the data came from the bdlim package.

## Examples

``` r
sbd_dlmtree <- get_sbd_dlmtree()

```
