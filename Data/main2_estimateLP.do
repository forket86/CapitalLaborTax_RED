/*******************************************************************************
 This program estimates local projections on fiscal variables. 
 The code is adapted from the code of Ramey and Zubairy (2018). 
 The data are yearly fiscal and economic variables, loaded from data_out.csv
 which is created by the Matlab file main1_dataprep.m
 The results are saved in "lp_results_sample`sample'_lags`p'.csv" which can 
 then be plotted using the Matlab file main3_plotLP.m
 
 !! IMPORTANT: User must manually change loadpath to the appropriate path for
               where they have stored the replication package !!
*******************************************************************************/



set more off
set varabbrev off
pause on


// TO BE CHANGED BY USER: Set this path to be the folder where the Data folder
// of the replication package is located. For example, if a user on a Mac
// computer called alex saved the replication package to their Downloads folder,
// the path could be:
// local loadpath "/Users/alex/Downloads/Replication Package/Data"
// On a Windows machine, the path could be:
// local loadpath "C:/Users/alex/Downloads/Replication Package/Data"
// Insert the correct path here:
local loadpath "REPLACE_WITH_YOUR_OWN_PATH/Replication Package/Data"
// set working directory to loadpath
cd "`loadpath'"

set matsize 800



/*******************************************************************************
  SET PARAMETERS THAT GOVERN SPECIFICATION
  Baseline specification in the text: sample = 1, p = 2, hmax = 10
*******************************************************************************/

local sample = 1  /*1 = full sample, 2 and 3 are defined below  */

local p = 2 /*number of lags of control variables*/

local hmax = 10 // how far forward to run projections


*******************************************************************************
** RAW DATA IMPORTATION AND DATA SETUP
*******************************************************************************


// import
import delimited data_out.csv, varnames(1) clear

// set date
tsset date, y

//check sorted by date
sort date

//*****************************
// SAMPLE SELECTION

// Original FVV moment sample (1971-2013)
if `sample'==2 {
	//drop pre-WWII
	drop if date<1971
	drop if date>2013
}

// Original Burnside IRF sample (1947-1995)
if `sample'==3 {
	//drop pre-WWII
	drop if date<1947
	drop if date>1995
}




//*****************************


// time t runs from 1 onwards
gen t = _n
// powers of t
gen t2 = t^2
gen t3 = t^3
gen t4 = t^4

//h is the horizon of the IRF
gen h = t - 1

// derived variables
gen ly = ln(y)
gen lg = ln(g)
gen lb = ln(rdebt_pub)


//list of main Y variables (also in control list)
global varlist ly lg taul tauk lb

//list of controls
global Zlist L(1/`p').newsy L(1/`p').ly L(1/`p').lg L(1/`p').taul L(1/`p').tauk L(1/`p').lb t t2 t3 t4

//list of other Y variables (not in main control list)
global varlist_extra tauc taup



*******************************************************************************
** RUN REGRESSIONS
*******************************************************************************


// store coef and se for each horizon here (this is the final output)
foreach var in $varlist $varlist_extra {
	gen b_`var' = .
	gen se_`var' = .
}

gen b_temp = .
gen se_temp = .

// regress for each horizon i and variable var
forvalues i = 0/`hmax' { 

	foreach var in $varlist {
		
		//regression and save coefficients in temp variables
		ivreg2 F`i'.`var' newsy $Zlist , robust bw(auto)
		replace b_temp = _b[newsy]
		replace se_temp = _se[newsy]
		
		//pull temp variables into the final output
		replace b_`var' = b_temp if h==`i'
		replace se_`var' = se_temp if h==`i'
  
	}
	
	
	foreach var in $varlist_extra {
		
		//regression and save coefficients in temp variables
		ivreg2 F`i'.`var' newsy $Zlist L(1/`p').`var', robust bw(auto)
		replace b_temp = _b[newsy]
		replace se_temp = _se[newsy]
		
		//pull temp variables into the final output
		replace b_`var' = b_temp if h==`i'
		replace se_`var' = se_temp if h==`i'
  
	}

}  

// drop temp variables
drop b_temp se_temp



*******************************************************************************
** EXPORT RESULTS
*******************************************************************************


// save temporary dataset so can reload later
save "temp", replace

// keep only h, coefs and standard errors
keep h b_* se_*
drop if h > `hmax'

//export to csv with sample number appended to file name
format * %15.10f // IMPORTANT: outsheet exports in display format, so need this to avoid rounding
format h %15.0f
outsheet using "lp_results_sample`sample'_lags`p'.csv" , replace comma nolabel

//reload full workspace
use temp.dta, clear
erase temp.dta













 
