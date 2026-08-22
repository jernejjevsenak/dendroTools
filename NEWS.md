# dendroTools 1.2.18

## New features

* `daily_response()` can now detrend the aggregated climate series of each
  moving window with the smoothers used by `dplR::detrend()`. Besides the
  existing `"SLD"`, `dc_method` accepts `"Spline"`, `"Mean"`, `"Ar"` and
  `"Friedman"`, controlled by the new arguments `dc_nyrs`, `dc_f`,
  `dc_difference`, `dc_span`, `dc_bass` and `dc_scale`.

* For strictly positive series the results are identical to `dplR::detrend()`.
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

## Bug fixes

* `dc_method = "SLD"` no longer recycles fitted values when the climate series
  contains missing values.

* `$optimized_return_all` keeps its year row names when detrending is applied.

* Warnings about missing values are no longer triggered by scattered gaps, such
  as DOY 366 in non-leap years, which the window aggregation removes. They now
  report only the years containing a gap at least as long as the shortest
  window width, and name those years.

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
