# Simulate genotype frequencies for 00, 01, 10, and 11 (maternal-paternal)
freq = runif(4)
freq = freq/sum(freq)
# freq = c(0.25, 0.25, 0.25, 0.25)

# Genotypic effects
a = 1
d = 0
i = 0.5

# Vector of genetic values
g = c(-a, d+i, d-i, a)
g = g - sum(freq*g) # Center

## Genotype dosage
x = c(0, 1, 1, 2)
x = x - sum(freq*x) # Center

## Breeding value regressor 
# Genetic contribution of "a"
x_a = c(-1, 0, 0, 1)
x_a = x_a - sum(freq*x_a) # Center

## Dominance deviation regressor
# Genetic contribution of "d"
x_d = c(0, 1, 1, 0)
x_d = x_d - sum(freq*x_d) # Center

# Regression coefficient from regressing x_d on x_a
m = sum(freq*x_a*x_d) / sum(freq*x_a^2)

# Construct orthogonal regressor using lack-of-fit from regression of x_d on x_a
x_d = x_d - x_a*m

# Check orthogonality (should be zero within numeric precision)
sum(freq*x_a*x_d)

## Imprinting regressor
# Genetic contribution of "i"
x_i = c(0, 1, -1, 0)
x_i = x_i - sum(freq*x_i) # Center

# Regression coefficient from regressing x_i on x_a
m_a = sum(freq*x_a*x_i) / sum(freq*x_a^2)

# Regression coefficient from regressing x_i on x_d
m_d = sum(freq*x_d*x_i) / sum(freq*x_d^2)

# Construct orthogonal regressor using lack-of-fit
x_i = x_i - x_a*m_a - x_d*m_d

# Check orthogonality (all should be zero within numeric precision)
sum(freq*x_i*x_a)
sum(freq*x_i*x_d)

## Calculate variances

# Additive genetic variance
alpha = sum(freq*x_a*g) / sum(freq*x_a^2)
bv = x_a*alpha
sum(freq*bv^2) 

# Dominance genetic variance
beta = sum(freq*x_d*g) / sum(freq*x_d^2)
dd = x_d*beta
sum(freq*dd^2)

# Imprinting genetic variance
gamma = sum(freq*x_i*g) / sum(freq*x_i^2)
id = x_i*gamma
sum(freq*id^2) 




