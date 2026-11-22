# dendroTools 2.0.0

## New features

* `daily_response()` can now detrend the aggregated climate series of each
  moving window with the smoothers used by `dplR::detrend()`. Besides the
  existing `"SLD"`, `dc_method` accepts `"Spline"`, `"Mean"`, `"Ar"` and
  `"Friedman"`, controlled by the new arguments `dc_nyrs`, `dc_f`,
  `dc_difference`, `dc_span`, `dc_bass` and `dc_scale`.

* For strictly positive series the results are identical to `dplR::detrend()`
  in dplR 1.7.8 and later, which scales a prewhitened (`"Ar"`) series by the
  mean of the original series.
  The ARSTAN safeguards that dplR applies when a fitted curve is not strictly
  positive are not reproduced, because they would silently replace the selected
  method for climate variables that cross zero; such windows return `NA`
  instead.

* `dc_difference = FALSE` (the dplR default) performs ratio detrending and is
  rejected for variables with negative values, such as temperature in degrees
  Celsius. Use `dc_difference = TRUE` for these.

* The returned object gains a `$detrending` element recording the detrending
  method and settings that were used.

* Results for `dc_method = "SLD"` are unchanged.

* `monthly_response()`, `daily_response_seascorr()` and
  `monthly_response_seascorr()` gain the same climate detrending, with the
  same arguments, defaults and checks as `daily_response()`, and their results
  gain `$detrending`. Results without detrending are unchanged, and for
  `dc_method = "SLD"` the (partial) correlation surfaces and
  `$optimized_return` are unchanged.

* In the two seascorr functions the primary and the control windows are
  detrended with the same settings by default. The new argument
  `dc_env_data_control` replaces any of them for the control windows. For
  example, `dc_method = "Spline", dc_difference = TRUE,
  dc_env_data_control = list(dc_difference = FALSE)` detrends temperature
  (primary) by differences and precipitation (control) by ratios. Their
  `$detrending` has the elements `primary` and `control`.

* In all four functions the dplR methods are fitted in chronological order
  of the years in the row names, so that results do not depend on how the
  rows are sorted. This matters most for `"Ar"`, because the seascorr
  functions sort the years in descending order before the analysis, and
  `$optimized_return_all` uses the climate data in the order they were
  supplied. `"SLD"` keeps the row order, as before; its residuals do not
  depend on the direction of time.

* `dc_method = "Ar"` changes `cor_na_use` or `pcor_na_use` to
  `"pairwise.complete.obs"`, with a warning, unless it already excludes
  missing values. `daily_response()` previously switched to `"complete.obs"`,
  which stopped with "no complete element pairs" as soon as ratio detrending
  left a window empty (for example a short precipitation sum that is zero in
  one year) and made the cluster filter skip. For a single window both give
  the same correlation.

* All four functions share one check of the `dc_*` arguments, so the
  arguments behave identically everywhere.

* New vignette "Detrending_climate_data" describes the climate detrending
  methods and gives examples for monthly and daily data, ratio and
  difference detrending, prewhitening and the seascorr functions.

* New vignette "Partial_correlations_seascorr" introduces
  `monthly_response_seascorr()` and `daily_response_seascorr()`, including
  the partial correlations of the control variable (`pcor_control`),
  `$variance_partitioning`, significance filters and climate detrending.

* `daily_response_seascorr()` and `monthly_response_seascorr()` gain
  `pcor_control`. With `pcor_control = TRUE` they also calculate, for every
  window, the partial correlation of the response with the control variable,
  controlling for the primary variable. This is what a run with the primary
  and control data, and their settings, swapped returns. These partial
  correlations are filtered with the same significance method and settings as
  the primary ones, with the permutation methods against the same
  permutations, and returned as `$calculations_control` with
  `$boot_lower_control`, `$boot_upper_control`,
  `$significance_filter_control` and a heatmap, `plot(x, type = 3)`.
  `summary()` adds the strongest of them, its window and the value in the
  optimal window of the primary variable. The primary results do not change.

* The seascorr functions return `$variance_partitioning` for the optimal
  window: the three simple correlations, both partial correlations, the
  variance of the response explained by both climate variables together, and
  how it splits into parts unique to each variable and a shared part.

* The significance filters of the two seascorr functions now share one
  implementation. Their results are unchanged.

* `aggregate_function = "SPEI"` is a new option in `daily_response()` and
  `monthly_response()`, and for the primary and the control variable of
  `daily_response_seascorr()` and `monthly_response_seascorr()`. It calculates
  the Standardized Precipitation-Evapotranspiration Index (SPEI) of every
  window from the climatic water balance in the climate data, i.e.
  precipitation minus potential evapotranspiration. The water balance is
  summed over the window and standardised with a log-logistic distribution
  fitted by unbiased probability weighted moments, the default of
  `SPEI::spei()`. Each window is standardised against all years of the
  climate data, before they are matched with the response. A monthly window
  of K months ending in month m therefore has the value of
  `SPEI::spei(scale = K)` in month m, the SPEI of a year does not depend on
  the analysed years, and a window in a previous year is standardised as the
  same window in the current year. Previous years, climate detrending (by
  differences only), the significance filters, the bootstrap and
  `pcor_control` work with SPEI as with the other aggregations. Unlike them,
  SPEI does not skip missing values: a window that includes a missing day or
  month is `NA` in that year, with a warning, and `cor_na_use` or
  `pcor_na_use` is changed to `"pairwise.complete.obs"` unless it already
  excludes missing values. DOY 366 of non-leap years is skipped. The SPEI
  package is not needed.

* New function `water_balance()` calculates the daily or monthly climatic
  water balance from precipitation and minimum and maximum temperatures, with
  potential evapotranspiration from the Hargreaves equation. Daily data use
  the extraterrestrial radiation of each day of year; for monthly data the
  values are those of `SPEI::hargreaves()`.

* New datasets from the black pine (*Pinus nigra*) site Sežana in Slovenia
  (Izmir et al. 2024, Dendrochronologia 88, 126236): `SEZ_rwl`, the raw ring
  widths of 41 cores from 21 trees (1882–2021), and `SEZ_daily_tmin`,
  `SEZ_daily_tavg`, `SEZ_daily_tmax` and `SEZ_daily_psum`, the daily minimum,
  mean and maximum temperatures and precipitation sums from E-OBS
  (1950–2020).

* New vignette "SPEI_drought_index" shows, with the data from Sežana, how to
  calculate the climatic water balance with `water_balance()` and how to use
  `aggregate_function = "SPEI"` with monthly and daily data, subsets of
  years, climate detrending and partial correlations.

## Changes in behaviour

* With `boot = TRUE`, all four functions now report the coefficient of the
  observed data (correlation, partial correlation, R squared or adjusted R
  squared), as with `boot = FALSE`. This value enters `$calculations`, the
  plots and the choice of the optimal window. Earlier versions reported the
  mean of the bootstrap replicates, which differs from the coefficient of the
  observed data and depends on the random resamples. The confidence intervals
  in `$boot_lower` and `$boot_upper` are unchanged. Results with
  `boot = TRUE` therefore change, and the optimal window can change with
  them.

## Bug fixes

* `dc_method = "SLD"` no longer recycles fitted values when the climate series
  contains missing values. This now also applies to `monthly_response()`,
  `daily_response_seascorr()` and `monthly_response_seascorr()`.

* `$optimized_return_all` keeps its year row names when detrending is applied.
  This now also applies to `monthly_response()`, `daily_response_seascorr()`
  and `monthly_response_seascorr()`. In the seascorr functions,
  `dc_method = "SLD"` previously also matched the primary and control series of
  `$optimized_return_all` by row position instead of by year and sorted the
  rows as text ("1", "10", "100", ...), so values from different years were
  paired. They are now merged by year.

* With `dc_method = "SLD"`, `$aggregated_climate` (and
  `$aggregated_climate_primary` and `$aggregated_climate_control`) is now a
  matrix, as documented and as without detrending. It used to be a data frame.
  The values are unchanged.

* Warnings about missing values, in all four functions, are no longer
  triggered by scattered gaps that the window aggregation removes, such as
  DOY 366 in non-leap years, or by gaps in years that are not analysed, such
  as an incomplete recent year after the end of the tree-ring record. They
  now report only the analysed years containing a gap at least as long as
  the shortest window width, and name those years. The seascorr functions
  report them separately for `env_data_primary` and `env_data_control`.

* With `aggregate_function = "sum"`, a window without any value in a year is
  now `NA`, in all four functions. It used to be 0, i.e. a dry window that was
  never observed. In `daily_response()` and `daily_response_seascorr()`,
  `"min"` and `"max"` returned `Inf` and `-Inf` for such windows, with a
  warning, and now also return `NA`. Windows with some values are aggregated
  over the available days or months, as before. Results for data without
  such gaps are unchanged.

* Without detrending, results are again exactly those of 1.2.17: integer
  climate data aggregated with `"sum"`, `"min"`, `"max"` or `"median"` were
  returned as doubles in `$optimized_return` and `$aggregated_climate` of
  `daily_response()`. The values were unchanged.

* `remove_insignificant_method = "sidak"` (all four functions) estimated the
  effective number of tests only from windows without a missing year. After
  `dc_method = "Ar"`, which removes the first years of nearly every window,
  that left a handful of windows (3 to 41 in the examples tested), so the
  number of tests was underestimated and the filter was too permissive. A
  year missing in every window removed every window and switched the filter
  off. Incomplete years are now left out instead of incomplete windows, and
  windows are dropped only if fewer than 10 complete years would remain.
  Results for complete data are unchanged.

* `$cross_validation` of `daily_response()` and `monthly_response()` returned
  `NA` for every fold that contained a year without a climate value, for
  example a missing year or the first years after `dc_method = "Ar"`. Such
  years are now left out before the folds are formed. Results for complete
  data are unchanged.

* The p values in `$temporal_stability` of `daily_response_seascorr()` and
  `monthly_response_seascorr()` ignored the degree of freedom spent on the
  control variable (`nrow(dataset_temp - 2)` evaluated to the number of
  rows), so they were slightly too small. They now use the number of complete
  years minus one, as `psych::corr.p()` asks for partial correlations. This
  changes these p values for every analysis.

* The p values in `$temporal_stability` of `monthly_response()` ignored
  `exact = FALSE`, which was misspelled in the calls to `cor.test()`, so for
  `cor_method = "kendall"` and `"spearman"` they were exact p values whenever
  the data had no ties. They are now the asymptotic p values, as in
  `daily_response()`. Pearson correlations are not affected.

* The documentation of `pcor_na_use` in the seascorr functions gave
  `"pairwise.complete.obs"` as the default. The default is `"everything"`, as
  in earlier versions.

* `boot_ci_type = "bca"` with `boot_n` smaller than the number of analysed
  years gave only `NA` intervals without explanation, because
  `boot::boot.ci()` cannot estimate bca intervals from so few replicates. All
  four functions now give a warning that names the `boot_n` needed. The
  documentation of `boot_ci_type` now lists the available interval types;
  `"stud"` was listed but is not available.

# dendroTools 1.2.17

## Correlation-specific critical values (changes Kendall results)

`critical_r()` derived its threshold from the Pearson sampling distribution
regardless of which coefficient was actually computed, and the criteria that
rely on it — `"alpha"`, `"fdr_BH"`, `"fdr_BY"` and `"sidak"` — inherited that
assumption. Kendall's tau has a different null distribution, so those criteria
were far too strict on tau surfaces and discarded genuine climate responses.

`critical_r()` gains a `cor_method` argument, accepting `"pearson"` (the
default, unchanged), `"spearman"` (which uses the same transformation applied
to ranks) and `"kendall"`. For Kendall's tau the critical value now follows the
normal approximation with variance 2(2n + 5) / 9n(n − 1), the same null
distribution `stats::cor.test()` uses when `exact = FALSE`. The four response
functions pass their `cor_method` or `pcor_method` through, and the false
discovery rate criteria build their p values from the matching distribution.

**Results change for anyone using `cor_method = "kendall"` or
`pcor_method = "kendall"` together with `remove_insignificant = TRUE`.** More
windows are now retained, correctly so. The magnitude depends on sample size:
the old threshold was 41 per cent too strict at n = 24 and 47 per cent too
strict at n = 73. On one test surface the `"sidak"` threshold fell from 0.443,
which retained no windows at all, to 0.276, against 0.297 from the permutation
test that had always been correct.

`"maxstat"` and `"cluster"` were never affected, because they build their null
by permuting the statistic actually being computed and so make no
distributional assumption. They remain the most reliable choice for rank
correlations, and for partial correlations with Kendall's tau in particular,
where the effective-sample-size adjustment used by the analytic criteria is
itself an approximation.

## Significance filtering: new interface

The logical argument `remove_insignificant_boot` has been removed and replaced
by `remove_insignificant_method`, which selects the criterion applied when
`remove_insignificant = TRUE`. This keeps one switch for whether to filter and
one for how, rather than adding an argument per method.

`remove_insignificant_method` accepts:

* `"alpha"` — per-window significance threshold. This is the default and
  reproduces the behaviour of earlier versions.
* `"boot"` — the bootstrap confidence interval excludes zero. Requires
  `boot = TRUE`.
* `"fdr_BY"` — false discovery rate, Benjamini–Yekutieli. Valid under
  arbitrary dependence between windows.
* `"fdr_BH"` — false discovery rate, Benjamini–Hochberg. Valid under positive
  dependence and less conservative than BY.
* `"maxstat"` — max-statistic permutation test.
* `"sidak"` — Šidák correction on the effective number of tests.
* `"cluster"` — cluster-extent inference.

The `"alpha"` and `"boot"` criteria are applied one window at a time and
therefore do not account for how many windows were screened. The remaining
criteria do.

`"maxstat"` builds a null distribution of the largest absolute correlation
obtained anywhere in the search by repeatedly shuffling the response, and
retains only windows exceeding the resulting threshold. The number of
permutations is set by the new argument `perm_n` (default 999). It makes no
distributional assumptions, but requires `perm_n` re-evaluations of the whole
search. It is available for `method = "cor"` and for `method = "lm"` with a
single predictor, but not for `method = "brnn"`, where the permutation would
require refitting the network at every window.

`"sidak"` targets the same family-wise error rate as `"maxstat"` but estimates
it analytically. Because neighbouring windows overlap, the number of windows
screened greatly overstates the number of independent tests. The effective
number is obtained from the eigenvalues of the correlation matrix among the
screened windows, following Li and Ji (2005), and a Šidák correction is applied
to that number rather than to the nominal count. The correlation matrix has
rank at most the number of years, so its non-zero eigenvalues are computed from
a small years-by-years matrix; no permutations are needed and the criterion is
effectively instantaneous. It is available for every method, including
`method = "brnn"`, where `"maxstat"` cannot run.

`"sidak"` is an approximation to `"maxstat"`. Across 54 combinations of site and
climate variable the two thresholds agreed to within 0.016 on average, with a
maximum discrepancy of 0.038, for ordinary correlations; for partial
correlations the corresponding figures were 0.020 and 0.046. The approximation
is close to exact for strongly autocorrelated variables such as temperature and
tends to be slightly permissive for weakly autocorrelated variables such as
daily precipitation, because the eigenvalue approach does not represent
temporal autocorrelation. `"maxstat"` is therefore preferable for final results,
while `"sidak"` suits exploratory work and `method = "brnn"`.

`"cluster"` is cluster-extent inference, adapted from neuroimaging. The surface
is thresholded at the new argument `cluster_threshold`, connected regions of
supra-threshold windows are identified, and each region's mass, the sum of
absolute correlations inside it, is compared against the largest mass arising
in permuted surfaces. It controls the family-wise error rate over *clusters*
rather than over windows, so it establishes that a region contains signal
without identifying which windows within it carry that signal. It is available
for `method = "cor"`, because the permuted statistic has to be laid out on the
same surface as the observed one.

`"cluster"` has little discriminating power on tree-ring response surfaces, and
is provided mainly for completeness and comparison. A real climate signal forms
one large connected region whose mass greatly exceeds anything the permuted
surfaces produce, so on daily surfaces every cluster that forms is retained and
the permutation step never rejects; the number of windows kept is then set by
`cluster_threshold` rather than by the data. On a daily search of about 3,400
windows, raising the threshold from the nominal critical correlation to 0.45
reduced the retained area from 95 to 40 per cent while the max-statistic
threshold on the same data stayed near 0.36. On the much smaller monthly
surfaces the permutation occasionally does reject a cluster, but not stably
between runs. `"maxstat"` or `"sidak"` are preferable when the question is
which windows carry signal.

For the seascorr functions the degrees of freedom account for the control
variable throughout, and the effective number of tests is estimated from the
primary windows.

Applies to `daily_response()`, `monthly_response()`,
`daily_response_seascorr()` and `monthly_response_seascorr()`.

The returned list element `$boot_ci_filter` has been renamed
`$significance_filter`. It records which method was applied, the criterion
used, and how many matrix cells were retained or removed. `"alpha"`, `"sidak"`,
`"maxstat"` and `"cluster"` also report the threshold applied, as
`critical_r`. For `"maxstat"` it additionally reports the observed maximum
(`observed_max_r`), the global p value (`p_global`), the number of
permutations, and the number of windows searched. For `"sidak"` it reports the
effective number of tests (`effective_tests`), the per-test significance level
(`alpha_per_test`), and the number of windows searched. For `"cluster"` it
reports the forming threshold used (`forming_threshold`), the critical cluster
mass (`cluster_mass_critical`), how many clusters formed and how many were
retained (`n_clusters_formed`, `n_clusters_retained`), the largest observed
mass, and the max-statistic threshold from the same permutations
(`maxstat_critical_r`), so the two criteria can be compared directly.

When `"maxstat"` removes every window, the message now reports the global p
value and states that no correlation exceeded what the search alone would
produce, rather than suggesting that the threshold be relaxed.

## Missing values: consistent behaviour across functions

Missing values are no longer handled silently. All four functions now default
to `cor_na_use = "everything"` or `pcor_na_use = "everything"`, so a window
overlapping a missing value returns `NA` and the user decides explicitly how to
proceed. Previously the seascorr functions defaulted to
`pcor_na_use = "pairwise.complete.obs"` and quietly computed over the available
years, which was inconsistent with `daily_response()` and `monthly_response()`.

All four functions now warn when the climate data contain missing values and
the chosen rule does not handle them, naming the argument to set. Previously
only `daily_response()` did so.

For `method = "cor"` and for the seascorr partial correlations, a window whose
bootstrap produces any missing replicate is now `NA` in all four functions,
rather than being computed from the remaining replicates. For `method = "lm"`
and `method = "brnn"` the opposite rule applies, again in all four functions: a
missing replicate there indicates a model fit that did not converge rather than
missing data, so it is dropped and the window is computed from the replicates
that succeeded.

Windows are still aggregated with `na.rm = TRUE`, so a window with a small
number of missing days is aggregated from the days that are present.

## Bootstrap performance

The bootstrap for `method = "cor"`, and for the partial correlations in the
seascorr functions, is now computed by matrix algebra rather than by
resampling loops. Each replicate correlation is obtained from weighted sums
over the original observations, so a whole window's replicates come from a
single matrix product.

Approximate speedups on a search of about 6500 windows with 1000 replicates:

| correlation | before | after |
|---|---|---|
| Pearson | ~3 min | ~4 s |
| Spearman | ~9 min | ~45 s |
| Kendall | ~19 min | ~55 s |

Confidence intervals of type `"norm"`, `"perc"` and `"basic"` are now computed
directly from the replicates and reproduce `boot::boot.ci()` exactly. The
`"bca"` type, and the `"lm"` and `"brnn"` methods, are unchanged and continue
to use `boot::boot()`.

**Results will differ slightly from earlier versions even with the same seed.**
Resample indices are now drawn once and shared by every window, instead of
being redrawn for each window. Differences are of the order of Monte Carlo
error (typically below 0.02 in correlation units). The shared resamples mean
all windows are evaluated on the same resampled years, so bootstrap surfaces
are internally comparable rather than each window carrying independent noise.

`cor_na_use` and `pcor_na_use` are honoured by the new code path, including
`"everything"`, under which any replicate that draws a missing value returns
`NA`, matching `stats::cor()`.

## Bug fixes

* `boot_ci_type = "stud"` is now rejected, with a message explaining that
  studentized intervals require a variance estimate for each bootstrap
  replicate, which is not available for correlation statistics. Previously
  this combination was accepted, ran the full search, and then failed part
  way through with the uninformative error `number of items to replace is not
  a multiple of replacement length`. `boot_ci_type` is validated once at the
  start of the call, so an invalid value now fails immediately rather than
  after minutes of computation.

* When `remove_insignificant = TRUE`, the confidence-interval matrices
  `$boot_lower` and `$boot_upper` are now masked with the same criterion as
  `$calculations`. Previously they retained values for windows the filter had
  removed, so plotting intervals over a filtered surface could suggest
  significance for windows that had in fact been rejected.

## References

Kendall, M. G. (1938). A new measure of rank correlation. *Biometrika* 30,
81–93.

Li, J. and Ji, L. (2005). Adjusting multiple testing in multilocus analyses
using the eigenvalues of a correlation matrix. *Heredity* 95, 221–227.

Nichols, T. E. and Holmes, A. P. (2002). Nonparametric permutation tests for
functional neuroimaging: a primer with examples. *Human Brain Mapping* 15,
1–25.

# dendroTools 1.2.16
* Added new quantile-probability arguments:
  - `quantile_prob` for `daily_response()` and `monthly_response()`
  - `quantile_prob_env_data_primary` and `quantile_prob_env_data_control` for `daily_response_seascorr()` and `monthly_response_seascorr()`
* This allows analyses based on lower or upper climate extremes, for example the 5th, 90th, or 95th percentile of temperature, precipitation, or other environmental variables within the selected time window.
* Extended previous-year analyses in `daily_response()`, `daily_response_seascorr()`, `monthly_response()`, and `monthly_response_seascorr()`.
* Added the new argument `number_previous_years`, which allows users to include more than one previous year in response-function analyses.
* Users can now analyse climate windows spanning up to five previous years plus the current year.
* The argument `previous_year` now acts as the main switch for previous-year analyses. If `previous_year = FALSE`, previous-year settings are ignored and the analysis is performed for the current year only.
* If negative values are supplied in `day_interval` or `month_interval` while `previous_year = FALSE`, the interval is ignored with a warning and the analysis is performed for the current year only.
* Updated plotting methods to support multi-year climate windows. Heatmaps and extreme-value plots now show relative-year labels such as `Y-3`, `Y-2`, `Y-1`, and `Y`, with year separators where appropriate.
* Updated x-axis labels in response-function plots to indicate whether the reference position is the start, end, or middle of the analysed window.
* Updated summary methods so that optimal windows from multi-year analyses are reported with the correct relative year.
* Improved robustness of the `brnn` method. Occasional unstable BRNN model fits are now handled more safely, returning `NA` for failed windows instead of interrupting the full analysis.
* Improved BRNN cross-validation handling so that failed folds return correctly structured missing-value outputs.
* Updated documentation and examples to demonstrate the new quantile aggregation option and the extended previous-year functionality.

# dendroTools 1.2.15
* there was a minor bug in sd(), which returned NA in case of missing values in climate data. This is now corrected and users are informed about the missing values with specific warning

# dendroTools 1.2.14
* there was a minor bug in aggregating climate-data when creatng the dataset for optimized_return
* Authors@R field was added to DESCRIPTION 

# dendroTools 1.2.13
* all examples have been streamlined to shorten the time required for executing R package checks
* added detailed explanations for the necessity of using donttest{}

# dendroTools 1.2.12
* PCA option is removed from the daily_response and other similar functions
* plot specific options is also removed
* New climate detrending option is implemented, 'SLD', which stands for 'Simple Linear Detrending' and should be preferred
* Two new arguments are implemented, 'skip_window_position' and 'skip_window_length' which allow for controlling the granularity of the analysis and can help in reducing computation time by focusing on a subset of the data
* Spline detrending is now removed since it caused instability of the dendroTools package

# dendroTools 1.2.11
* response data frame can now also have a missing value, but only if "cor_na_use" argument is one of "complete.obs", "na.or.complete" or "pairwise.complete.obs"
* an undocumented argument in internal function is now removed

# dendroTools 1.2.10
* A bug in some special cases for the reference window "end" is now removed
* In case of no significant correlations, daily_response() now returns a warning
rather than error

# dendroTools 1.2.9
* new options for the aggregate function - minimum and maximum
* additional feature on reference_window in the monthly_response() and monthly_response_seascorr()

# dendroTools 1.2.8
* daily_response() default remove_insignificant argument is set to FALSE
* updated description
* default method in daily_response() and monthly_response() is set to 'cor'
* corrected all spelling errors

# dendroTools 1.2.7
* comments are added to the code, which is now more readable
* unused part of the scripts are removed
* new arguments are implemented to control the behavior of calculations of correlation coefficients in case of missing values
* the calculation of significance of partial correlation coefficients is improved
* daily_transform() is also improved, so it gives the correct behavior in case of missing values

# dendroTools 1.2.6
* detrending climate data is now available in all response type functions

# dendroTools 1.2.4
* updated vignettes
* new arguments are added which can be used to detrend climate
* climate detrending is based on the dplR R package
* gridExtra is removed from dependences

# dendroTools 1.2.1
* Minor correction for ggplot2 object -> guide = FALSE --> guide = 'none'

# dendroTools 1.2.0
* There is a major upgrade related to new arguments that now allow you to specify the time of interest to be used for calculating correlation coefficients or other statistical metrics.
For example, if you are using daily_response() or daily_response_seascorr() and want to calculate correlations only from the beginning of last October to the current end of September, use the argument daily_interval =c(-274, 273). Note that a negative sign indicates the previous year's doy, while a positive sign indicates the current year doy. In normal (non-leap ) years, October 1 represents doy 274, and September 30 represents doy 273.
* If you use monthly_response() or monthly_response_seascorr() and want to calculate correlations only from the beginning of the previous October to the current end of September, for example, use the argument monthly_interval =c(-10, 9). Note that a negative sign indicates the month of the previous year, while a positive sign indicates the month of the current year.
Other functions have been changed only slightly, so they will work even if new arguments are used.

# dendroTools 1.1.4
* new error is implemented in case of wrong subset_years
* minor corrections to partial correlation functions

# dendroTools 1.1.3
* deleted example in vignette Examples_with_daily_climate_data which caused problems with CRAN check

# dendroTools 1.1.2
* a bug corrected in the data_transform()
* minor bug is removed which happened in a case of non significant calculations or only one significant correlation

# dendroTools 1.1.1
* dendroTools version is changed to 1.1.1
* Corrections related to testthat changes.
* imports is updated: oce(>= 1.2-0)

# dendroTools 1.0.8
* Minor corrections after new publication
* Sys.setlocale("LC_TIME", "C") is now added, so optimal dates will always show in English language
* The vignette titles are changed

# dendroTools 1.0.7.
* Bootstrapping is now available for all metrics, monthly and daily functions
* Seacorr is now available also for monthly data, including bootstrapping
* swit272 daily precipitation dataset is added
* swit272 monthly temperature dataset is removed
* KNMI_daily_transform() is now called data_transform(). The entire function is updated and now also enables transformation of daily data into monthly.

# dendroTools 1.0.6.
* Bootstrapping of correlation coefficients is introduced. To do so, use the argument boot in daily_response() and monthly_response(). Bootstrapping is currently available only for correlation coefficients. In future version, it will also be available for model fitting. 
* Package version is changed to 1.0.6

# dendroTools 1.0.5.
* new function: monthly_response_seascorr()
* titles of plot_extreme, plot_specific and plot_heatmap are updated
* for daily_response() and monthly_response(), it is now possible to define method for correlation coefficient: "pearson", "kendall", "spearman". To do so, use , cor_method argument
* daily_response output list is class "dmrs"
* Using newly defined class "dmrs", new summary function is defined
* Package version is changed to 1.0.5

# dendroTools 1.0.4.
* New data is available for examples demonstration: swit272_monthly_temperatures
* New function is available: daily_response_seascorr(). This function implements partial correlations in the analysis of daily response functions
* Package version is changed to 1.0.4

# dendroTools 1.0.3.
* Package version is changed to 1.0.3
* There is new function introduced: KNMI_daily_transform() which transforms daily data obtained from KNMI Climate explorer into data frame suitable for daily_response(). 
* New TRW chronology is included, swit272
* New function: monthly_response()
* New data: LJ_monthly_temperatures
* New data: LJ_monthly_precipitation

# dendroTools 1.0.2.
* Package version is changed to 1.0.2
* There have been many issues with RWeka package, which depends on rJava. Therefore, functions from RWeka are replaced with other functions. For random forest model we now use randomForest package
* There are six new output elements, all of them are residual diagnostic plots for calibration, holdout and edge data.
* Due to many problems related to rJava in R, MT and BMT methods from RWeka are now joined and replaced by cubist function from Cubist r package. Cubist fits a model tree or ensemble of model trees, if committees argument is more than 1.
