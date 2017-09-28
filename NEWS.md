# dendroTools 0.0.3.

* Package was previously called dendroExtra, but authors decided to rename it and submit it as a new package

* Changes in comparison to the package dendroExtra 0.0.2
* new dataset is added that is used to test new function compare_methods()
* two new functions are added: compare_methods() and calculate_measures()
* required packages were added to the Description: dcv (>= 0.1.1), MLmetrics (>= 1.1.1), RWeka (>= 0.4-34),
        dplyr (>= 0.7.0), reshape (>= 0.8.6)
* the brnn option from the daily_response() was removed. We recived some comments that this function does not converge in some examples,
so we decided to remove it. Consequently, also smooth_matrix() was removed since it is not needed anymore. 
