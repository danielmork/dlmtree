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
#> MCMC iterations (est time: 6 seconds)
#> [0--------25--------50--------75--------100]
#>  ''''''''''''''''''''''''''''''''''''''''''
#> Compiling results...
pip(fit)
#>   mod_num   mod_bin mod_scale        c1        c2        c3        c4        c5 
#>     1.000     0.295     1.000     0.360     0.350     0.375     0.435     0.395 
#>        b1        b2        b3        b4        b5 
#>     0.320     0.300     0.300     0.330     0.315 
pip(fit, type = 2)
#>          var1      var2   pip
#> 3     mod_num mod_scale 1.000
#> 27  mod_scale   mod_num 1.000
#> 29  mod_scale mod_scale 0.635
#> 1     mod_num   mod_num 0.430
#> 33  mod_scale        c4 0.225
#> 81         c4 mod_scale 0.225
#> 34  mod_scale        c5 0.175
#> 94         c5 mod_scale 0.175
#> 7     mod_num        c4 0.165
#> 79         c4   mod_num 0.165
#> 6     mod_num        c3 0.150
#> 66         c3   mod_num 0.150
#> 32  mod_scale        c3 0.145
#> 68         c3 mod_scale 0.145
#> 10    mod_num        b2 0.140
#> 31  mod_scale        c2 0.140
#> 55         c2 mod_scale 0.140
#> 118        b2   mod_num 0.140
#> 8     mod_num        c5 0.125
#> 36  mod_scale        b2 0.125
#> 92         c5   mod_num 0.125
#> 120        b2 mod_scale 0.125
#> 30  mod_scale        c1 0.120
#> 39  mod_scale        b5 0.120
#> 42         c1 mod_scale 0.120
#> 159        b5 mod_scale 0.120
#> 37  mod_scale        b3 0.115
#> 133        b3 mod_scale 0.115
#> 5     mod_num        c2 0.110
#> 11    mod_num        b3 0.110
#> 53         c2   mod_num 0.110
#> 131        b3   mod_num 0.110
#> 2     mod_num   mod_bin 0.100
#> 9     mod_num        b1 0.100
#> 14    mod_bin   mod_num 0.100
#> 35  mod_scale        b1 0.100
#> 105        b1   mod_num 0.100
#> 107        b1 mod_scale 0.100
#> 38  mod_scale        b4 0.085
#> 146        b4 mod_scale 0.085
#> 4     mod_num        c1 0.075
#> 16    mod_bin mod_scale 0.075
#> 28  mod_scale   mod_bin 0.075
#> 40         c1   mod_num 0.075
#> 12    mod_num        b4 0.070
#> 144        b4   mod_num 0.070
#> 85         c4        c4 0.065
#> 13    mod_num        b5 0.060
#> 59         c2        c4 0.060
#> 83         c4        c2 0.060
#> 157        b5   mod_num 0.060
#> 91         c4        b5 0.055
#> 163        b5        c4 0.055
#> 20    mod_bin        c4 0.045
#> 43         c1        c1 0.045
#> 46         c1        c4 0.045
#> 57         c2        c2 0.045
#> 58         c2        c3 0.045
#> 62         c2        b2 0.045
#> 64         c2        b4 0.045
#> 70         c3        c2 0.045
#> 71         c3        c3 0.045
#> 75         c3        b2 0.045
#> 80         c4   mod_bin 0.045
#> 82         c4        c1 0.045
#> 122        b2        c2 0.045
#> 123        b2        c3 0.045
#> 148        b4        c2 0.045
#> 156        b4        b5 0.045
#> 168        b5        b4 0.045
#> 17    mod_bin        c1 0.040
#> 18    mod_bin        c2 0.040
#> 19    mod_bin        c3 0.040
#> 26    mod_bin        b5 0.040
#> 41         c1   mod_bin 0.040
#> 54         c2   mod_bin 0.040
#> 67         c3   mod_bin 0.040
#> 77         c3        b4 0.040
#> 90         c4        b4 0.040
#> 99         c5        c5 0.040
#> 104        c5        b5 0.040
#> 149        b4        c3 0.040
#> 150        b4        c4 0.040
#> 158        b5   mod_bin 0.040
#> 164        b5        c5 0.040
#> 52         c1        b5 0.035
#> 63         c2        b3 0.035
#> 65         c2        b5 0.035
#> 76         c3        b3 0.035
#> 117        b1        b5 0.035
#> 129        b2        b4 0.035
#> 135        b3        c2 0.035
#> 136        b3        c3 0.035
#> 153        b4        b2 0.035
#> 160        b5        c1 0.035
#> 161        b5        c2 0.035
#> 165        b5        b1 0.035
#> 21    mod_bin        c5 0.030
#> 50         c1        b3 0.030
#> 72         c3        c4 0.030
#> 84         c4        c3 0.030
#> 86         c4        c5 0.030
#> 87         c4        b1 0.030
#> 88         c4        b2 0.030
#> 93         c5   mod_bin 0.030
#> 98         c5        c4 0.030
#> 111        b1        c4 0.030
#> 124        b2        c4 0.030
#> 130        b2        b5 0.030
#> 134        b3        c1 0.030
#> 166        b5        b2 0.030
#> 25    mod_bin        b4 0.025
#> 44         c1        c2 0.025
#> 48         c1        b1 0.025
#> 49         c1        b2 0.025
#> 56         c2        c1 0.025
#> 73         c3        c5 0.025
#> 97         c5        c3 0.025
#> 102        c5        b3 0.025
#> 103        c5        b4 0.025
#> 108        b1        c1 0.025
#> 121        b2        c1 0.025
#> 138        b3        c5 0.025
#> 145        b4   mod_bin 0.025
#> 151        b4        c5 0.025
#> 47         c1        c5 0.020
#> 51         c1        b4 0.020
#> 60         c2        c5 0.020
#> 74         c3        b1 0.020
#> 89         c4        b3 0.020
#> 95         c5        c1 0.020
#> 96         c5        c2 0.020
#> 100        c5        b1 0.020
#> 110        b1        c3 0.020
#> 112        b1        c5 0.020
#> 128        b2        b3 0.020
#> 137        b3        c4 0.020
#> 140        b3        b2 0.020
#> 143        b3        b5 0.020
#> 147        b4        c1 0.020
#> 167        b5        b3 0.020
#> 22    mod_bin        b1 0.015
#> 24    mod_bin        b3 0.015
#> 61         c2        b1 0.015
#> 78         c3        b5 0.015
#> 101        c5        b2 0.015
#> 106        b1   mod_bin 0.015
#> 109        b1        c2 0.015
#> 125        b2        c5 0.015
#> 132        b3   mod_bin 0.015
#> 142        b3        b4 0.015
#> 154        b4        b3 0.015
#> 162        b5        c3 0.015
#> 45         c1        c3 0.010
#> 69         c3        c1 0.010
#> 115        b1        b3 0.010
#> 139        b3        b1 0.010
#> 23    mod_bin        b2 0.005
#> 114        b1        b2 0.005
#> 116        b1        b4 0.005
#> 119        b2   mod_bin 0.005
#> 126        b2        b1 0.005
#> 152        b4        b1 0.005
#> 15    mod_bin   mod_bin 0.000
#> 113        b1        b1 0.000
#> 127        b2        b2 0.000
#> 141        b3        b3 0.000
#> 155        b4        b4 0.000
#> 169        b5        b5 0.000
# }
```
