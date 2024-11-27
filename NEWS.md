# MMVmalaria 1.1.0

* Longitudinal chemoprevention modelling update.

* New functions
-getTimeKRaboveGR : Function to calculate the time that a defined kill rate is above a defined growth rate.
-generate_bite_time_poisson : Function to generate a string of bite times based on a poisson distribution.
-add_bites_to_trial_files : Function to generate a series of trial files from some defined "baseline" (containing PKPD parameters) and add infectious bites events based on defined parameters such as amount of parasites innoculated and a user-defined function describing the distribution of bite times.
-summarizeTrial_ChemoSurvival_MonoPD : Summary function for Simulate_VirtualTrials. Performs the summary as defined by summarizeTrial_ChemoSurvival, and additionally calculates time above MIC and time above MPC90.
-summarizeTrial_ChemoSurvival_ComboPD : Summary function for Simulate_VirtualTrials. Performs the summary as defined by summarizeTrial_ChemoSurvival, and additionally calculates the time that kill rate is above growth rate. 

* Fixes
-summarizeTrial_ChemoSurvival : Fixed a critical bug where outputNames were not being passed into the function correctly. 

* Updated remotes
-MMVverse/MMVbase added as remote in DESCRIPTION file to allow MMVbase to be imported. 
 