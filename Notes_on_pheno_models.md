# Notes on Phenotype Observation Models

## Date and status

- Drafted on April 14, 2026
- Updated on April 15, 2026
- Low-level expert helpers now reflect part of this discussion in the main repo:
  - `asLogNormal(x, meanLogShift = NULL)`
  - `asPoisson(x, meanLogShift = NULL)`
  - `asCategorical(x, p = ..., threshold = ...)`, with `p` taking precedence when both are supplied
- The higher-level stored observation-model API proposed below remains a design note and is not implemented

## Purpose

This note captures a proposed redesign of the current `asX()` approach for non-Gaussian observed phenotypes. The aim is to keep the statistical model coherent while making the API easier to use and more consistent with how AlphaSimR already treats founder-based trait scaling and founder-based residual calibration.

The immediate motivation was the addition of an `asPoisson()` transform, but the proposal applies more broadly to:

- Gaussian observed traits
- log-normal observed traits
- categorical observed traits
- Poisson observed traits
- TODO: what-other-distribution-we-should-consider observed traits
- TODO: can we add user defined traits - how would this be done?

## Problem statement

The current `asX()` functions are useful expert tools, but they are awkward as the main user interface because they expose latent-scale parameters that most users do not think in.

Current low-level pattern:

- `asLogNormal()` uses `meanLogShift` as an additional shift on the latent log scale
- `asCategorical()` can use observed-scale category probabilities `p` or latent `threshold`; if both are supplied, `p` currently takes precedence
- `asPoisson()` uses `meanLogShift` as an additional shift on the latent log scale

This creates two problems:

1. The API is inconsistent across distributions
2. The user must reason about latent-scale calibration even when they mostly care about observed-scale outcomes

There is also an execution-hook design question - we can do transforms via `finalizePop()` or `finalizePheno()`, but these are advanced tools, so ideally we would have something slicker/simpler.

## Core design principle

Use a two-layer phenotype model API:

- low-level expert transforms remain available as `asLogNormal()`, `asCategorical()`, `asPoisson()`
- high-level user-facing phenotype observation models are stored in `SimParam` and applied automatically whenever phenotypes are generated
- this stored-model layer should integrate with the existing `finalizePheno()` execution point
- `finalizePop()` remains the broader population hook, while `finalizePheno()` is the correct hook for phenotype-value recoding

The high-level API should let the user specify targets on the observed scale and should compile those into fixed latent-scale parameters.

## Temporal semantics

This is the most important design rule.

User-facing phenotype model parameters should be interpreted as setup-time targets, usually relative to a reference population, and then should remain fixed as simulation proceeds.

That means:

- users set the phenotype observation model once, usually at the start
- the model is calibrated on a reference population, usually `founderPop`
- this calibration produces fixed latent-scale parameters
- those fixed parameters are reused over time
- current observed means, variances, and category frequencies are then allowed to evolve naturally as populations evolve genetically

This matches current AlphaSimR behavior for:

- `setVarE()`, which converts `h2` / `H2` into fixed absolute error variances using founder-scale variances
- `rescaleTraits()`, which explicitly targets founder-population means and variances

### What should stay fixed

- distribution family
- link function
- structural choice such as whether latent Gaussian residual is included
- compiled latent-scale parameters such as:
  - log-link family `meanLogShift`
  - categorical thresholds

### What should evolve over time

- current population mean genetic value
- current population genetic variance
- current observed mean
- current observed variance
- current category frequencies

### What should not happen by default

The system should not silently recalibrate the observation model on every generation or every call to `setPheno()`. That would erase evolutionary change on the observed scale.

If recalibration is ever supported, it should be explicit and advanced.

## Statistical framing

### General latent-to-observed pipeline

For a given trait, phenotype generation should be thought of as:

1. compute latent genetic value
2. optionally add latent Gaussian residual
3. apply observation model to get observed phenotype

This separates:

- latent biological signal
- latent Gaussian noise
- observed-scale sampling / recoding

### Poisson case

For Poisson, the natural model is:

```text
y_i | x_i ~ Poisson(lambda_i)
log(lambda_i) = meanLogShift + x_i
```

where `x_i` is the latent value actually used by the model.

Important distinction:

- conditional on `x_i`, the response is Poisson
- marginally, if `x_i` is Gaussian, the response is Poisson-lognormal

Therefore Poisson sampling does not replace latent environmental deviation. It adds observation noise conditional on the rate. A Gaussian residual in the latent predictor creates additional rate heterogeneity and overdispersion on the observed scale.

This implies a user-facing switch is needed:

- pure Poisson sampling around a latent predictor
- Poisson-lognormal with latent Gaussian residual plus Poisson sampling

Recommended default for Poisson:

- `includeVarE = FALSE`

Reason:

- Poisson already contributes sampling variance
- this gives the cleanest default interpretation
- extra heterogeneity can still be requested explicitly

## Public API proposal

### Stored model in `SimParam`

Proposed primary user API:

```r
SP$setPhenoModel(
  trait = 1,
  model = poissonModel(mean = 10, reference = "founder", includeVarE = FALSE)
)

SP$setPhenoModel(
  trait = 2,
  model = categoricalModel(
    p = c(0.2, 0.5, 0.3),
    reference = "founder",
    includeVarE = TRUE
  )
)

SP$setPhenoModel(
  trait = 3,
  model = lognormalModel(mean = 1, reference = "founder", includeVarE = TRUE)
)
```

Supporting methods:

```r
SP$setPhenoModel(trait, model)
SP$getPhenoModel(trait = NULL)
SP$clearPhenoModel(trait = NULL)
SP$recompilePhenoModel(trait = NULL)
```

### Model constructor functions

```r
gaussianModel()

lognormalModel(
  mean = NULL,
  meanLogShift = NULL,
  reference = "founder",
  includeVarE = TRUE
)

categoricalModel(
  p = NULL,
  threshold = NULL,
  reference = "founder",
  includeVarE = TRUE
)

poissonModel(
  mean = NULL,
  meanLogShift = NULL,
  reference = "founder",
  includeVarE = FALSE
)
```

### API rules

- for log-link families, exactly one user-scale or expert-scale parameterization must be supplied
- user-scale parameters:
  - `mean`
  - `p`
- expert-scale parameters:
  - `meanLogShift`
  - `threshold`
- for categorical models, a user-friendly option is to keep the current low-level behavior where `p` takes precedence over `threshold`
- `reference` defaults to `"founder"`
- `includeVarE` specifies whether latent Gaussian residual variance participates in the latent predictor used by the observation model

### One-off override in `setPheno()`

Optional extension:

```r
pop = setPheno(
  pop,
  h2 = 0.3,
  obsModel = poissonModel(mean = 8, includeVarE = FALSE),
  simParam = SP
)
```

This should be treated as an override for the current call, not as a replacement for the model stored in `SimParam`.

## Meaning of user-facing parameters

### Reference-based interpretation

User-facing parameters should be interpreted relative to the reference population at compile time.

Examples:

- `poissonModel(mean = 10)` means:
  - compile a fixed `meanLogShift` so the reference population has average observed count about 10
- `categoricalModel(p = c(0.1, 0.8, 0.1))` means:
  - compile thresholds so the reference population has those category frequencies
- `lognormalModel(mean = 1)` means:
  - compile a fixed `meanLogShift` so the reference population has observed mean about 1

The compiled parameters then remain fixed while later populations evolve.

## Internal design proposal

### Where model state should live

Observation-model state should live in `SimParam`, and should be executed through the phenotype-finalization path.

Reason:

- `finalizePheno()` now already provides the correct execution point for phenotype recoding across new and existing populations
- this removes the earlier need to invent a second phenotype-finalization hook
- but the observation-model state itself should still be structured data in `SimParam`, not hidden inside ad hoc user code
- therefore the design should combine:
  - structured model state in `SimParam`
  - package-managed model application at phenotype finalization time
  - optional user customization through `finalizePheno()`

### Proposed internal state per trait

Each trait should have a phenotype-model record with at least:

```r
list(
  family = "poisson",
  spec = list(
    mean = 10,
    meanLogShift = NULL,
    reference = "founder",
    includeVarE = FALSE
  ),
  compiled = list(
    meanLogShift = 1.234
  ),
  dirty = FALSE
)
```

Or, more generally:

```r
list(
  family = <character>,
  spec = <user specification>,
  compiled = <latent-scale parameters>,
  reference = <reference definition>,
  dirty = <logical>
)
```

### Proposed internal helper functions

```r
compilePhenoModel(model, trait, simParam, varE = NULL, reps = 1, ...)
applyPhenoModel(x, model, rng = TRUE, ...)
getLatentForModel(pop, trait, includeVarE, ...)
```

Possible R-level naming if package-private:

```r
.compilePhenoModel(...)
.applyPhenoModel(...)
.getLatentPhenoInput(...)
```

## Compile-once, apply-many-times contract

This should be the core implementation contract.

### Compile step

The compile step:

- evaluates the model against a reference population or reference latent distribution
- translates user-facing targets into fixed latent-scale parameters
- stores the compiled parameters in `SimParam`

### Apply step

The apply step:

- uses the current population's latent values
- uses the already compiled fixed latent parameters
- produces observed phenotypes

No recalibration occurs during normal phenotype generation.

## When compiled models become dirty

Compiled phenotype models should be invalidated when latent-scale trait architecture changes in ways that affect calibration.

At minimum:

- `setVarE()`
- `rescaleTraits()`
- `switchTrait()`
- trait addition
- trait removal
- change of `founderPop`

Recommended behavior:

- preserve user `spec`
- mark compiled model `dirty = TRUE`
- lazily recompile on next phenotype generation, or explicitly through `SP$recompilePhenoModel()`

## Calibration formulas

### Empirical calibration preferred

For user-scale targets, calibration should preferably use the empirical reference population rather than rely only on Normal approximations.

### Log-normal target mean

If latent values used by the model for the reference population are `x_ref`, then:

```r
meanLogShift = log(target_mean / mean(exp(x_ref)))
```

### Poisson target mean

If latent values used by the model for the reference population are `x_ref`, then:

```r
meanLogShift = log(target_mean / mean(exp(x_ref)))
```

Then at application time:

```r
lambda = exp(meanLogShift + x_current)
y = rpois(n = length(lambda), lambda = lambda)
```

### Categorical target probabilities

For a target probability vector `p`, compile thresholds from the empirical or modeled reference latent distribution so that the reference population approximately realizes those category frequencies.

Two options:

- empirical quantiles from `x_ref`
- Normal-theory quantiles using reference mean and variance

Recommended default:

- empirical quantiles, because that aligns better with the actual latent reference population

## Semantics of `includeVarE`

This argument must be very clear in both docs and implementation.

### `includeVarE = TRUE`

The latent predictor for the observation model includes Gaussian residual variance.

Examples:

- log-normal from latent phenotype
- ordered categorical from latent phenotype
- Poisson-lognormal from latent phenotype

### `includeVarE = FALSE`

The latent predictor for the observation model excludes Gaussian residual variance.

Examples:

- Poisson counts around the genetic predictor only
- direct transform from genetic value without latent Gaussian environmental noise

## Integration into phenotype pipeline

### Desired pipeline around `finalizePheno()`

Conceptually:

```r
setPheno(...)
  -> compute latent trait values
  -> add Gaussian residual if requested
  -> apply registered phenotype observation models if present
  -> apply user `finalizePheno()` customizations
  -> return/store observed phenotype
```

The exact internal implementation may choose to compute both:

- latent value without Gaussian residual
- latent value with Gaussian residual

and then select between them according to `includeVarE`.

### Recommended composition with the current hook

Recommended order:

```r
pheno = calcPheno(...)
pheno = .applyRegisteredPhenoModels(pheno, pop = pop, simParam = simParam, ...)
pheno = simParam$finalizePheno(pheno, pop = pop, simParam = simParam, ...)
```

This keeps:

- built-in observation models as package behavior
- `finalizePheno()` as an advanced post-processing hook for users
- compatibility with the new main-repo phenotype finalization design

### Why not only `finalizePop()`

`finalizePop()` is useful as a general hook, but not sufficient for a persistent phenotype observation system because it only applies on newly created populations and is not the central phenotype-generation pipeline.

### Why not rely only on ad hoc `finalizePheno()` code

`finalizePheno()` is now the right execution point, but a user should not have to hand-code:

- founder/reference calibration
- compile-once versus apply-many-times logic
- dirty-state invalidation
- composition across multiple traits and model families

Those concerns should remain package-managed, with `finalizePheno()` available for extra customization.

## Suggested defaults by family

```r
gaussianModel()
```

- identity transform

```r
lognormalModel(..., includeVarE = TRUE)
```

- default should use the latent phenotype

```r
categoricalModel(..., includeVarE = TRUE)
```

- default should use the latent phenotype

```r
poissonModel(..., includeVarE = FALSE)
```

- default should be pure Poisson sampling around the latent predictor

## Expert layer retained

The existing `asX()` functions should remain available as low-level tools.

Proposed role:

- `asLogNormal()`: low-level latent-scale transform
- `asCategorical()`: low-level transform using thresholds or category probabilities
- `asPoisson()`: low-level transform using a latent log-scale shift

But they should no longer be the primary user-facing phenotype-model API.

## Minimal spec for `asPoisson()`

If `asPoisson()` exists as an expert tool, it should remain simple:

```r
asPoisson(x, meanLogShift = NULL)
```

where:

- `x` is the latent predictor on the log scale
- `meanLogShift` is a fixed additive offset on the log scale
- output is sampled from `Poisson(exp(meanLogShift + x))`

This function should not try to do full founder/reference calibration by itself. That belongs in the stored-model layer.

## Naming note for low-level expert helpers

The low-level expert helpers should prefer names that make two things explicit:

- this is an additional shift applied during recoding
- the scale of that shift

This is why plain names such as `mean` or `meanlog` can be misleading:

- `mean` sounds like an observed-scale mean
- `meanlog` sounds like a full latent mean, while in practice it is usually an additional shift applied to `x`

For log-link families, a clearer naming pattern is:

- `meanLogShift`

This works well for:

- `asLogNormal()`
- `asPoisson()`
- any future log-link count or positive-valued helper

Among the naming options discussed so far:

- `shiftMeanLog` is awkward and harder to read
- `meanLogShift` is clearer and keeps the important word `shift`

This is now the direction used by the current low-level expert helpers:

```r
asLogNormal(x, meanLogShift = NULL)
asPoisson(x, meanLogShift = NULL)
```

Semantics:

- `x` already contains the latent trait values from AlphaSimR
- `meanLogShift` is an extra additive shift to the log-scale location during recoding

Even with such an argument available, the primary place to set trait mean in normal AlphaSimR workflows remains:

- `SP$addTrait(..., mean = ...)`

The low-level helper shift should usually be left at default unless the user intentionally wants an extra recoding-specific shift.

For ordered categorical traits, a different pattern remains preferable:

- expert parameterization by `threshold`
- convenience parameterization by `p`
- current low-level helper behavior gives `p` precedence over `threshold` for user convenience

so there is less need for a shared shift-style name there.

## Other distributions to consider

Beyond Gaussian, log-normal, categorical, and Poisson, the most useful future candidates are likely:

- negative binomial
- gamma
- beta
- Bernoulli / binomial
- zero-inflated count models

### Negative binomial

Use case:

- overdispersed counts
- a natural extension beyond Poisson when Poisson sampling variance is too restrictive

Parameters to be mindful of:

- log-scale location shift, ideally using the same naming convention as Poisson
  - e.g. `meanLogShift`
- one dispersion parameter
  - likely `size` or `theta`

Candidate low-level shape:

```r
asNegBinom(x, meanLogShift = NULL, size)
```

### Gamma

Use case:

- strictly positive continuous traits with right skew

Parameters to be mindful of:

- log-scale location shift
  - e.g. `meanLogShift`
- one shape-like parameter
  - likely `shape`
  - alternatively a coefficient-of-variation style parameter

Candidate low-level shape:

```r
asGamma(x, meanLogShift = NULL, shape)
```

### Beta

Use case:

- proportions constrained to `(0, 1)`

Parameters to be mindful of:

- logit-scale location shift
  - e.g. `meanLogitShift`
- one precision parameter
  - typically `phi`

Candidate low-level shape:

```r
asBeta(x, meanLogitShift = NULL, phi)
```

### Bernoulli / binomial

Use case:

- binary traits
- bounded count outcomes

Parameters to be mindful of:

- logit-scale location shift
  - e.g. `meanLogitShift`
- binomial trial count for the bounded-count case
  - e.g. `size`

Note:

- binary threshold traits are already partly covered by `asCategorical()`
- an explicit Bernoulli / binomial helper could still be useful if sampling under a GLM interpretation is desired

Candidate low-level shapes:

```r
asBernoulli(x, meanLogitShift = NULL)
asBinomial(x, meanLogitShift = NULL, size)
```

### Zero-inflated count models

Use case:

- count traits with excess zeros beyond standard Poisson or negative binomial behavior

Parameters to be mindful of:

- count-model location shift
  - e.g. `meanLogShift`
- zero-inflation probability or logit-scale shift
  - e.g. `zeroProb` or `zeroLogitShift`
- count dispersion if using a negative-binomial count component

Candidate low-level shapes:

```r
asZIPoisson(x, meanLogShift = NULL, zeroProb)
asZINegBinom(x, meanLogShift = NULL, size, zeroProb)
```

### General parameterization guideline

For future low-level expert helpers, aim for:

- one location shift on the link scale
- one family-specific dispersion / shape / precision parameter when needed
- names that make the link scale explicit

This suggests a useful naming pattern across families:

- `meanLogShift` for log-link families
- `meanLogitShift` for logit-link families
- `threshold` for ordered categorical traits

## Proposed code format

Below is a suggested code shape only. It is not yet implemented.

### User-facing constructors

```r
poissonModel <- function(
  mean = NULL,
  meanLogShift = NULL,
  reference = "founder",
  includeVarE = FALSE
) {
  stopifnot(xor(is.null(mean), is.null(meanLogShift)))

  structure(
    list(
      family = "poisson",
      spec = list(
        mean = mean,
        meanLogShift = meanLogShift,
        reference = reference,
        includeVarE = includeVarE
      )
    ),
    class = "AlphaSimR_pheno_model"
  )
}
```

```r
lognormalModel <- function(
  mean = NULL,
  meanLogShift = NULL,
  reference = "founder",
  includeVarE = TRUE
) {
  stopifnot(xor(is.null(mean), is.null(meanLogShift)))

  structure(
    list(
      family = "lognormal",
      spec = list(
        mean = mean,
        meanLogShift = meanLogShift,
        reference = reference,
        includeVarE = includeVarE
      )
    ),
    class = "AlphaSimR_pheno_model"
  )
}
```

```r
categoricalModel <- function(
  p = NULL,
  threshold = NULL,
  reference = "founder",
  includeVarE = TRUE
) {
  if (!is.null(p)) {
    threshold = NULL
  }

  structure(
    list(
      family = "categorical",
      spec = list(
        p = p,
        threshold = threshold,
        reference = reference,
        includeVarE = includeVarE
      )
    ),
    class = "AlphaSimR_pheno_model"
  )
}
```

### SimParam storage API

```r
setPhenoModel = function(trait, model) {
  stopifnot(inherits(model, "AlphaSimR_pheno_model"))
  private$.phenoModels[[trait]] = list(
    family = model$family,
    spec = model$spec,
    compiled = NULL,
    dirty = TRUE
  )
  invisible(self)
}
```

```r
getPhenoModel = function(trait = NULL) {
  if (is.null(trait)) {
    return(private$.phenoModels)
  }
  private$.phenoModels[trait]
}
```

```r
clearPhenoModel = function(trait = NULL) {
  if (is.null(trait)) {
    private$.phenoModels = vector("list", self$nTraits)
  } else {
    private$.phenoModels[trait] = list(NULL)
  }
  invisible(self)
}
```

### Package-side application before user `finalizePheno()`

```r
.applyRegisteredPhenoModels <- function(pheno, pop, simParam, ...) {
  for (trt in seq_len(ncol(pheno))) {
    model = simParam$getPhenoModel(trt)[[1]]

    if (is.null(model)) {
      next
    }

    if (isTRUE(model$dirty) || is.null(model$compiled)) {
      xRef = .getReferenceLatentInput(
        trait = trt,
        model = model,
        simParam = simParam
      )
      model$compiled = .compilePhenoModel(model, xRef = xRef)
      model$dirty = FALSE
      simParam$setPhenoModel(trt, model)
    }

    x = .getLatentInputForCurrentPheno(
      pheno = pheno,
      pop = pop,
      trait = trt,
      model = model,
      simParam = simParam
    )

    pheno[, trt] = .applyPhenoModel(x = x, model = model)
  }

  pheno
}
```

### Compile helper skeleton

```r
.compilePhenoModel <- function(model, xRef) {
  if (model$family == "poisson") {
    if (!is.null(model$spec$meanLogShift)) {
      compiled = list(meanLogShift = model$spec$meanLogShift)
    } else {
      compiled = list(
        meanLogShift = log(model$spec$mean / mean(exp(xRef)))
      )
    }
  } else if (model$family == "lognormal") {
    if (!is.null(model$spec$meanLogShift)) {
      compiled = list(meanLogShift = model$spec$meanLogShift)
    } else {
      compiled = list(
        meanLogShift = log(model$spec$mean / mean(exp(xRef)))
      )
    }
  } else if (model$family == "categorical") {
    if (!is.null(model$spec$threshold)) {
      compiled = list(threshold = model$spec$threshold)
    } else {
      compiled = list(
        threshold = as.numeric(stats::quantile(
          xRef,
          probs = cumsum(model$spec$p)
        ))
      )
    }
  } else if (model$family == "gaussian") {
    compiled = list()
  } else {
    stop("Unknown phenotype model family")
  }

  compiled
}
```

### Apply helper skeleton

```r
.applyPhenoModel <- function(x, model) {
  if (model$family == "gaussian") {
    return(x)
  }

  if (model$family == "lognormal") {
    return(exp(model$compiled$meanLogShift + x))
  }

  if (model$family == "categorical") {
    return(as.numeric(cut(
      x = x,
      breaks = model$compiled$threshold,
      include.lowest = TRUE,
      right = FALSE
    )))
  }

  if (model$family == "poisson") {
    lambda = exp(model$compiled$meanLogShift + x)
    return(stats::rpois(n = length(lambda), lambda = lambda))
  }

  stop("Unknown phenotype model family")
}
```

### Phenotype generation skeleton

```r
.finalizeGeneratedPheno <- function(pheno, pop, simParam, ...) {
  pheno = .applyRegisteredPhenoModels(
    pheno = pheno,
    pop = pop,
    simParam = simParam,
    ...
  )

  pheno = simParam$finalizePheno(
    pheno,
    pop = pop,
    simParam = simParam,
    ...
  )

  pheno
}
```

## Open design questions

### 1. Should reference always be founder-based?

Recommended default:

- yes, `"founder"`

Optional advanced support:

- explicit population object as reference

### 2. Should category calibration use empirical quantiles or Normal theory?

Recommended default:

- empirical quantiles from the reference latent values

### 3. Should Poisson support `reps` in a special way?

Current recommendation:

- do not add special count-specific `reps` semantics in the first implementation
- keep current phenotype replication behavior separate from the observation-model redesign

### 4. Should observation-model summaries be exposed?

Likely yes:

```r
summary(SP$getPhenoModel(1)[[1]])
```

but not required for the first implementation.

### 5. Should built-in models run before or after user `finalizePheno()`?

Recommended answer:

- before

Reason:

- package-managed observation models are the structured default behavior
- user `finalizePheno()` should remain an advanced post-processing hook
- this avoids forcing users to manually reimplement the built-in model logic

## Recommended implementation order

1. Add internal storage for per-trait phenotype observation models in `SimParam`
2. Add user-facing model constructor functions
3. Add compile and apply helpers
4. Add package-managed phenotype-model application in the same phenotype-finalization path that now feeds `SimParam$finalizePheno`
5. Preserve `finalizePheno()` as a user post-processing hook
6. Keep `finalizePop()` examples as optional advanced customization, not the primary pattern
7. Add tests for:
   - founder-calibrated Poisson mean
   - founder-calibrated log-normal mean
   - founder-calibrated categorical frequencies
   - persistence of fixed compiled parameters across evolving populations
   - invalidation after `rescaleTraits()` / `setVarE()`
   - composition of built-in phenotype models with user-defined `finalizePheno()`

## Bottom line

The proposed system should be thought of as a persistent, founder-calibrated phenotype observation layer:

- users specify observed-scale goals once
- AlphaSimR compiles those goals into fixed latent-scale parameters
- later populations evolve naturally
- the observation model is stored in `SimParam`
- the built-in model logic is applied through the phenotype-finalization path
- `finalizePheno()` remains available as the user-facing advanced hook on top of that path

This preserves statistical clarity, improves usability, and fits the existing founder-based design of AlphaSimR better than using standalone `asX()` transforms as the primary user interface.
