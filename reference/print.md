# print

print generic function for S3method

## Usage

``` r
print(x, ...)

# S3 method for class 'tdlnm'
print(x, ...)

# S3 method for class 'tdlm'
print(x, ...)

# S3 method for class 'tdlmm'
print(x, ...)

# S3 method for class 'hdlm'
print(x, ...)

# S3 method for class 'hdlmm'
print(x, ...)

# S3 method for class 'monotone'
print(x, ...)

# S3 method for class 'summary.hdlm'
print(x, digits = 3, ...)

# S3 method for class 'summary.hdlmm'
print(x, digits = 3, ...)

# S3 method for class 'summary.monotone'
print(x, digits = 3, ...)

# S3 method for class 'summary.tdlm'
print(x, digits = 3, ...)

# S3 method for class 'summary.tdlmm'
print(x, digits = 3, ...)

# S3 method for class 'summary.tdlnm'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class 'tdlm', 'tdlmm', 'tdlnm', 'hdlm', 'hdlmm',
  'monotone', representing a fitted model using dlmtree(); or a summary
  object produced by applying summary() to one of these model objects.

- ...:

  additional parameters

- digits:

  number of decimal places to round the numeric values to

## Value

For a fitted model object, prints an assorted model output including
model formula call and available methods. For a summary object, prints a
summary output of a model fit in the R console.
