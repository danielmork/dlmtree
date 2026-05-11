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
#>     1.000     0.625     1.000     0.625     0.725     0.530     0.665     0.710 
#>        b1        b2        b3        b4        b5 
#>     0.485     0.555     0.645     0.570     0.600 
pip(fit, type = 2)
#>          var1      var2   pip
#> 3     mod_num mod_scale 1.000
#> 27  mod_scale   mod_num 1.000
#> 29  mod_scale mod_scale 1.000
#> 31  mod_scale        c2 0.300
#> 55         c2 mod_scale 0.300
#> 34  mod_scale        c5 0.255
#> 94         c5 mod_scale 0.255
#> 37  mod_scale        b3 0.245
#> 133        b3 mod_scale 0.245
#> 1     mod_num   mod_num 0.225
#> 16    mod_bin mod_scale 0.215
#> 28  mod_scale   mod_bin 0.215
#> 30  mod_scale        c1 0.205
#> 42         c1 mod_scale 0.205
#> 32  mod_scale        c3 0.200
#> 68         c3 mod_scale 0.200
#> 33  mod_scale        c4 0.190
#> 81         c4 mod_scale 0.190
#> 7     mod_num        c4 0.180
#> 38  mod_scale        b4 0.180
#> 79         c4   mod_num 0.180
#> 146        b4 mod_scale 0.180
#> 2     mod_num   mod_bin 0.170
#> 14    mod_bin   mod_num 0.170
#> 39  mod_scale        b5 0.170
#> 159        b5 mod_scale 0.170
#> 6     mod_num        c3 0.155
#> 36  mod_scale        b2 0.155
#> 66         c3   mod_num 0.155
#> 120        b2 mod_scale 0.155
#> 5     mod_num        c2 0.145
#> 53         c2   mod_num 0.145
#> 8     mod_num        c5 0.130
#> 11    mod_num        b3 0.130
#> 60         c2        c5 0.130
#> 92         c5   mod_num 0.130
#> 96         c5        c2 0.130
#> 131        b3   mod_num 0.130
#> 4     mod_num        c1 0.125
#> 40         c1   mod_num 0.125
#> 12    mod_num        b4 0.120
#> 45         c1        c3 0.120
#> 69         c3        c1 0.120
#> 103        c5        b4 0.120
#> 144        b4   mod_num 0.120
#> 151        b4        c5 0.120
#> 20    mod_bin        c4 0.115
#> 59         c2        c4 0.115
#> 80         c4   mod_bin 0.115
#> 83         c4        c2 0.115
#> 18    mod_bin        c2 0.110
#> 35  mod_scale        b1 0.110
#> 54         c2   mod_bin 0.110
#> 58         c2        c3 0.110
#> 70         c3        c2 0.110
#> 86         c4        c5 0.110
#> 98         c5        c4 0.110
#> 107        b1 mod_scale 0.110
#> 10    mod_num        b2 0.105
#> 47         c1        c5 0.105
#> 52         c1        b5 0.105
#> 63         c2        b3 0.105
#> 75         c3        b2 0.105
#> 91         c4        b5 0.105
#> 95         c5        c1 0.105
#> 100        c5        b1 0.105
#> 112        b1        c5 0.105
#> 118        b2   mod_num 0.105
#> 123        b2        c3 0.105
#> 135        b3        c2 0.105
#> 160        b5        c1 0.105
#> 163        b5        c4 0.105
#> 23    mod_bin        b2 0.100
#> 62         c2        b2 0.100
#> 73         c3        c5 0.100
#> 76         c3        b3 0.100
#> 88         c4        b2 0.100
#> 97         c5        c3 0.100
#> 119        b2   mod_bin 0.100
#> 122        b2        c2 0.100
#> 124        b2        c4 0.100
#> 136        b3        c3 0.100
#> 25    mod_bin        b4 0.095
#> 44         c1        c2 0.095
#> 46         c1        c4 0.095
#> 50         c1        b3 0.095
#> 56         c2        c1 0.095
#> 57         c2        c2 0.095
#> 65         c2        b5 0.095
#> 74         c3        b1 0.095
#> 82         c4        c1 0.095
#> 101        c5        b2 0.095
#> 110        b1        c3 0.095
#> 115        b1        b3 0.095
#> 125        b2        c5 0.095
#> 134        b3        c1 0.095
#> 139        b3        b1 0.095
#> 145        b4   mod_bin 0.095
#> 161        b5        c2 0.095
#> 13    mod_num        b5 0.090
#> 49         c1        b2 0.090
#> 51         c1        b4 0.090
#> 72         c3        c4 0.090
#> 78         c3        b5 0.090
#> 84         c4        c3 0.090
#> 99         c5        c5 0.090
#> 121        b2        c1 0.090
#> 147        b4        c1 0.090
#> 157        b5   mod_num 0.090
#> 162        b5        c3 0.090
#> 9     mod_num        b1 0.085
#> 19    mod_bin        c3 0.085
#> 21    mod_bin        c5 0.085
#> 67         c3   mod_bin 0.085
#> 93         c5   mod_bin 0.085
#> 102        c5        b3 0.085
#> 104        c5        b5 0.085
#> 105        b1   mod_num 0.085
#> 138        b3        c5 0.085
#> 164        b5        c5 0.085
#> 17    mod_bin        c1 0.080
#> 26    mod_bin        b5 0.080
#> 41         c1   mod_bin 0.080
#> 77         c3        b4 0.080
#> 85         c4        c4 0.080
#> 90         c4        b4 0.080
#> 116        b1        b4 0.080
#> 128        b2        b3 0.080
#> 140        b3        b2 0.080
#> 143        b3        b5 0.080
#> 149        b4        c3 0.080
#> 150        b4        c4 0.080
#> 152        b4        b1 0.080
#> 156        b4        b5 0.080
#> 158        b5   mod_bin 0.080
#> 167        b5        b3 0.080
#> 168        b5        b4 0.080
#> 22    mod_bin        b1 0.075
#> 87         c4        b1 0.075
#> 89         c4        b3 0.075
#> 106        b1   mod_bin 0.075
#> 111        b1        c4 0.075
#> 129        b2        b4 0.075
#> 130        b2        b5 0.075
#> 137        b3        c4 0.075
#> 153        b4        b2 0.075
#> 166        b5        b2 0.075
#> 24    mod_bin        b3 0.070
#> 114        b1        b2 0.070
#> 117        b1        b5 0.070
#> 126        b2        b1 0.070
#> 132        b3   mod_bin 0.070
#> 165        b5        b1 0.070
#> 64         c2        b4 0.065
#> 148        b4        c2 0.065
#> 43         c1        c1 0.060
#> 142        b3        b4 0.060
#> 154        b4        b3 0.060
#> 61         c2        b1 0.055
#> 109        b1        c2 0.055
#> 48         c1        b1 0.050
#> 71         c3        c3 0.050
#> 108        b1        c1 0.050
#> 15    mod_bin   mod_bin 0.000
#> 113        b1        b1 0.000
#> 127        b2        b2 0.000
#> 141        b3        b3 0.000
#> 155        b4        b4 0.000
#> 169        b5        b5 0.000
# }
```
