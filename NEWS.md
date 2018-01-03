# dendroTools 0.0.5.

* new feature is added to the daily_response(): the possibility of using rowMedians instead of rowMeans
* some functions were made internal: smooth_matrix(), plot_extreme(), plot_heatmap(), plot_specific(), round_df(), count_ones()
* new data is added: LJ_daily_precipitation
* the support for the tidy data for the env_data data frame is implemented
* new output elements are given: temporal_stability, corss_validation, results of PC regression (PC_output) and transfer_function
* plot_extreme(), plot_specifi() and plot_heatmap() are now integrated in daiily_reponse() and given as output elements.
* The word measure is replaced with word "metric". This way is more correct. 
* The compare_methods() and daily_response() have a new feature: the possibility of using PCA transformation. 
* New regression methods are added to the compare_methods(): Ridge and Lasso regression from the glmnet package and polynomial regression
* Titles of plot_extreme() and plot_specific() have new title attribute: Optimal Selection description
* The compare_methods() has a new feature: blocked_CV
* The daily_response() has a new feature: the possibility of using PCA on response data. New examples are added. 
* Package version is changed to 0.0.5
