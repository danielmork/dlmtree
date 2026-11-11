# Calculates posterior inclusion probabilities (PIPs) for modifiers in HDLM & HDLMM

Method for calculating posterior inclusion probabilities (PIPs) for
modifiers in HDLM & HDLMM

## Usage

``` r
pip(object, type = 1)
```

## Arguments

- object:

  An object of class dlmtree.

- type:

  Type=1 indicates single modifier PIPs. Type=2 indicates joint modifier
  PIPs for two modifiers.

## Value

numeric vector of PIPs named with modifiers (type=1) or data.frame of
PIPs with the following columns (type=2):

- var1:

  first modifier of joint modifiers

- var2:

  second modifier of joint modifiers

- pip:

  joint PIPs for the two modifiers

## Details

pip

## Examples

``` r
# \donttest{
# Posterior inclusion probability with HDLM 
D <- sim.hdlmm(sim = "B", n = 1000)
fit <- dlmtree(y ~ ., 
               data = D$dat,
               exposure.data = D$exposures,
               dlm.type = "linear",
               family = "gaussian",
               het = TRUE)
#> Preparing data...
#> 
#> Running shared HDLM:
#> Burn-in % complete 
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> MCMC iterations (est time: 4 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
pip(fit)
#>   mod_num   mod_bin mod_scale        c1        c2        c3        c4        c5 
#>     1.000     0.250     1.000     0.300     0.310     0.295     0.210     0.280 
#>        b1        b2        b3        b4        b5 
#>     0.270     0.240     0.230     0.295     0.395 
pip(fit, type = 2)
#>          var1      var2   pip
#> 3     mod_num mod_scale 1.000
#> 27  mod_scale   mod_num 1.000
#> 1     mod_num   mod_num 0.690
#> 29  mod_scale mod_scale 0.630
#> 39  mod_scale        b5 0.160
#> 159        b5 mod_scale 0.160
#> 4     mod_num        c1 0.155
#> 40         c1   mod_num 0.155
#> 13    mod_num        b5 0.150
#> 157        b5   mod_num 0.150
#> 30  mod_scale        c1 0.120
#> 42         c1 mod_scale 0.120
#> 38  mod_scale        b4 0.115
#> 146        b4 mod_scale 0.115
#> 36  mod_scale        b2 0.110
#> 120        b2 mod_scale 0.110
#> 6     mod_num        c3 0.105
#> 12    mod_num        b4 0.105
#> 34  mod_scale        c5 0.105
#> 66         c3   mod_num 0.105
#> 94         c5 mod_scale 0.105
#> 144        b4   mod_num 0.105
#> 31  mod_scale        c2 0.095
#> 55         c2 mod_scale 0.095
#> 35  mod_scale        b1 0.090
#> 107        b1 mod_scale 0.090
#> 32  mod_scale        c3 0.085
#> 68         c3 mod_scale 0.085
#> 5     mod_num        c2 0.080
#> 16    mod_bin mod_scale 0.080
#> 28  mod_scale   mod_bin 0.080
#> 53         c2   mod_num 0.080
#> 2     mod_num   mod_bin 0.070
#> 7     mod_num        c4 0.070
#> 9     mod_num        b1 0.070
#> 10    mod_num        b2 0.070
#> 14    mod_bin   mod_num 0.070
#> 79         c4   mod_num 0.070
#> 105        b1   mod_num 0.070
#> 118        b2   mod_num 0.070
#> 11    mod_num        b3 0.060
#> 33  mod_scale        c4 0.060
#> 81         c4 mod_scale 0.060
#> 131        b3   mod_num 0.060
#> 18    mod_bin        c2 0.055
#> 54         c2   mod_bin 0.055
#> 8     mod_num        c5 0.050
#> 92         c5   mod_num 0.050
#> 37  mod_scale        b3 0.040
#> 75         c3        b2 0.040
#> 78         c3        b5 0.040
#> 104        c5        b5 0.040
#> 123        b2        c3 0.040
#> 130        b2        b5 0.040
#> 133        b3 mod_scale 0.040
#> 162        b5        c3 0.040
#> 164        b5        c5 0.040
#> 166        b5        b2 0.040
#> 48         c1        b1 0.035
#> 61         c2        b1 0.035
#> 62         c2        b2 0.035
#> 88         c4        b2 0.035
#> 108        b1        c1 0.035
#> 109        b1        c2 0.035
#> 122        b2        c2 0.035
#> 124        b2        c4 0.035
#> 44         c1        c2 0.030
#> 52         c1        b5 0.030
#> 56         c2        c1 0.030
#> 57         c2        c2 0.030
#> 65         c2        b5 0.030
#> 77         c3        b4 0.030
#> 86         c4        c5 0.030
#> 98         c5        c4 0.030
#> 115        b1        b3 0.030
#> 117        b1        b5 0.030
#> 139        b3        b1 0.030
#> 149        b4        c3 0.030
#> 160        b5        c1 0.030
#> 161        b5        c2 0.030
#> 165        b5        b1 0.030
#> 20    mod_bin        c4 0.025
#> 23    mod_bin        b2 0.025
#> 43         c1        c1 0.025
#> 49         c1        b2 0.025
#> 50         c1        b3 0.025
#> 59         c2        c4 0.025
#> 63         c2        b3 0.025
#> 71         c3        c3 0.025
#> 74         c3        b1 0.025
#> 80         c4   mod_bin 0.025
#> 83         c4        c2 0.025
#> 89         c4        b3 0.025
#> 91         c4        b5 0.025
#> 100        c5        b1 0.025
#> 101        c5        b2 0.025
#> 103        c5        b4 0.025
#> 110        b1        c3 0.025
#> 112        b1        c5 0.025
#> 116        b1        b4 0.025
#> 119        b2   mod_bin 0.025
#> 121        b2        c1 0.025
#> 125        b2        c5 0.025
#> 134        b3        c1 0.025
#> 135        b3        c2 0.025
#> 137        b3        c4 0.025
#> 151        b4        c5 0.025
#> 152        b4        b1 0.025
#> 163        b5        c4 0.025
#> 19    mod_bin        c3 0.020
#> 46         c1        c4 0.020
#> 51         c1        b4 0.020
#> 58         c2        c3 0.020
#> 60         c2        c5 0.020
#> 67         c3   mod_bin 0.020
#> 70         c3        c2 0.020
#> 82         c4        c1 0.020
#> 85         c4        c4 0.020
#> 96         c5        c2 0.020
#> 129        b2        b4 0.020
#> 142        b3        b4 0.020
#> 143        b3        b5 0.020
#> 147        b4        c1 0.020
#> 153        b4        b2 0.020
#> 154        b4        b3 0.020
#> 156        b4        b5 0.020
#> 167        b5        b3 0.020
#> 168        b5        b4 0.020
#> 17    mod_bin        c1 0.015
#> 21    mod_bin        c5 0.015
#> 24    mod_bin        b3 0.015
#> 25    mod_bin        b4 0.015
#> 26    mod_bin        b5 0.015
#> 41         c1   mod_bin 0.015
#> 45         c1        c3 0.015
#> 47         c1        c5 0.015
#> 64         c2        b4 0.015
#> 69         c3        c1 0.015
#> 72         c3        c4 0.015
#> 76         c3        b3 0.015
#> 84         c4        c3 0.015
#> 90         c4        b4 0.015
#> 93         c5   mod_bin 0.015
#> 95         c5        c1 0.015
#> 102        c5        b3 0.015
#> 114        b1        b2 0.015
#> 126        b2        b1 0.015
#> 132        b3   mod_bin 0.015
#> 136        b3        c3 0.015
#> 138        b3        c5 0.015
#> 145        b4   mod_bin 0.015
#> 148        b4        c2 0.015
#> 150        b4        c4 0.015
#> 158        b5   mod_bin 0.015
#> 73         c3        c5 0.010
#> 87         c4        b1 0.010
#> 97         c5        c3 0.010
#> 99         c5        c5 0.010
#> 111        b1        c4 0.010
#> 22    mod_bin        b1 0.005
#> 106        b1   mod_bin 0.005
#> 128        b2        b3 0.005
#> 140        b3        b2 0.005
#> 15    mod_bin   mod_bin 0.000
#> 113        b1        b1 0.000
#> 127        b2        b2 0.000
#> 141        b3        b3 0.000
#> 155        b4        b4 0.000
#> 169        b5        b5 0.000
# }
```
