********************************************************************************
*                                                                              *
*                              TI UTE Assignment                               *
*                       Part I -  descriptives and estimation                  *
********************************************************************************

clear all
clear
set more off
eststo clear

* Set working directory

	// your code here
	cd "C:\Users\Regi Kusumaatmadja\Dropbox\TI year 2\Block 03\Urban code"
	global RESULTS "C:\Users\Regi Kusumaatmadja\Dropbox\TI year 2\Block 03\Urban code\results"

* import data
use "TI_UTE_London_data.dta", clear


*turn off automatic abbreviation by stata
set varabbrev off, permanent


* define group IDs for location FEs
egen msoai_n = group(msoai)
egen msoaj_n = group(msoaj)

* rename/define variables to allign with ASRW

gen Qi = price_per_m2
gen Ki = landha

* define share open space
gen Oi = openspaceha/areaha




********************************************************************************
*                     1 - Exploring data 
********************************************************************************

* 1.1 Inspect data  
	
	// your code here
	tabstat ncomm commtime Qi Ki Oi, stat(mean sd p1 p25 p50 p75 p99) col(stat)
	
* 1.2 make descriptives table  

	// your code here
	tabstat ncomm commtime landha Qi Ki Oi, stat(mean sd p1 p25 p50 p75 p99) col(stat)
	tabstat ncomm commtime landha price_per_m2, stat(mean sd p1 p25 p50 p75 p99) col(stat)
// 	bys laname: tabstat ncomm commtime landha price_per_m2, stat(mean sd p25 p50 p75) col(stat)
	

* 1.3 make graphs   

	// your code here
// 	gen alt_lncomm = ln(1 + ncomm)
// 	gen llandha = ln(landha)
// 	gen llprice_per_m2 = ln(price_per_m2)
	
	** one histogram
	*** Qi
	twoway (histogram Qi), by(laname, title("Histogram of Residential Price per Square Meter (Qi)"))
	
	graph export "$RESULTS/Q13c_01.png", as(png) replace
		
	*** Oi
	twoway (histogram Oi), by(laname, title("Histogram of Share of Open Space (Qi)"))
	
	graph export "$RESULTS/Q13c_02.png", as(png) replace
	
	*** Commtime
	twoway (histogram commtime), by(laname, title("Histogram of Commuting Time (commtime)"))
	
	graph export "$RESULTS/Q13c_03.png", as(png) replace

	** one scatter plot
	g lncomm_al = log(ncomm)
	reg lncomm_al commtime msoai_n msoaj_n, r nocons
	predict hat_lncomm_al
	label variable hat_lncomm_al "Predicted Log of Number of Commuters"
	
	twoway(scatter hat_lncomm_al commtime), by(laname, title("Predicted Log of Commuters v. Commuting Time: OLS (nocons)"))
	
	graph export "$RESULTS/Q13c_04.png", as(png) replace
	
	poi2hdfe ncomm commtime, id1(msoai_n) id2(msoaj_n)
	predict hat_ncomm
	
	twoway (scatter hat_ncomm commtime), by(laname, title("Number of Commuters v. Commuting Time: Poisson"))
	
	graph export "$RESULTS/Q13c_05.png", as(png) replace
	

********************************************************************************
*                     2 - Estimation 
********************************************************************************


********************************************************************************
*                     2.0 - Set parameters                                    
********************************************************************************
scalar ?? = 0.80    // Share of labour in firm costs (from Ahlfeldt et al)
scalar ?? = 0.75    // Share of expenses on other goods (from Ahlfeldt et al)
scalar ?? = 0.75    // Share of capital in construction costs (from Ahlfeldt et al)

scalar H_total = 4864000 // Total working population in London

scalar thresholdwages = 0.001 // to export to mata later



********************************************************************************
*                      2.1 Obtain commuting elasticity using poi2hdfe          *
********************************************************************************

	//your code here
	gen log_ncomm = log(ncomm)
	gen log_commtime = log(commtime)
	label variable log_ncomm "log of flows"
	eststo clear
	eststo ols01: reg log_ncomm commtime msoai_n msoaj_n, r nocons
	eststo poi01: poi2hdfe ncomm commtime, id1(msoai_n) id2(msoaj_n)
// 	eststo poi02: poi2hdfe log_ncomm commtime, id1(msoai_n) id2(msoaj_n)
	
	esttab ols01 poi01 using "$RESULTS/Q21.tex", se star(* 0.1 ** 0.05 *** 0.01) ///
		keep(commtime) ///
		title(OLS and Poisson Regressions \label{Q21Tab1}) stats(N r2_a) ///
		label booktabs replace
	eststo clear
	
	drop log_ncomm log_commtime
	
	// define elasticity for use below
	scalar ?? = -_b[commtime] //your code here
	
	
********************************************************************************
*    2.2 Obtain number of people working (H???) and living (H???) in each area     *
********************************************************************************

* Residential population
qui total ncomm
qui replace ncomm = ncomm*( H_total /_b[ncomm])   //Total number of jobs in Greater London in 2011

* calculate residents for each location
qui bysort msoai: egen H???i = sum(ncomm) // people living in i
qui label variable H???i "Number of workers living in location i"

qui gen  tempj_resident = H???i if msoai==msoaj // copy H???i when i == j
qui bysort msoaj: egen H???j = mean(tempj_resident) // take mean to copy value to all j in group
qui label variable H???j "Number of workers living in location j"

* jobs
qui bysort msoaj: egen H???j = sum(ncomm) // people working in j
qui label variable H???j "Number of workers working in location j"

qui gen  tempi_worker = H???j if msoai==msoaj // copy H???j when i == j
qui bysort msoai: egen H???i = mean(tempi_worker) //take mean to copy value to all i in group
qui label variable H???i "Number of workers working in location i"

* remove helper columns
qui drop tempi_worker tempj_resident


*make graphs with population data
	
	//your code here
	
// 	** standardize residential population and working population
// 	for any H???i H???j H???i H???j: egen zX = std(X)
//	
// 	** descriptive statistics
// 	tabstat H???i H???j, stat(mean sd median) col(stat)
// 	tabstat H???i H???j, stat(mean sd median) col(stat)
//	
// 	** residential (i vs j)
// 	twoway (kdensity zH???i) (kdensity zH???j)
//		
// 	** worker (i vs j)
// 	twoway (kdensity zH???i) (kdensity zH???j)
//	
// 	** origin v. destination
// 	twoway (kdensity zH???i) (kdensity zH???j)
// 	twoway (kdensity zH???j) (kdensity zH???i)
	
	** population/area 
	gen densityH???i = H???i/landha
	label variable densityH???i "Density of workers living in location i"
	gen densityH???i = H???i/landha
	label variable densityH???i "Density of workers working in location i"
	gen densityH???j = H???j/landha
	label variable densityH???j "Density of workers living in location j"
	gen densityH???j = H???j/landha 
	label variable densityH???j "Density of workers working in location j"
	tabstat H???i H???i H???j H???j densityH???i densityH???i densityH???j densityH???j, stat(mean median) col(stat)
	
	*** try to make histogram/kdensity of population area
	preserve 
// 	twoway (histogram densityH???i) (histogram densityH???i), by(laname) ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist_no.png", as(png) replace
	
	twoway  (histogram densityH???j, width(4) color(green%30)) ///
			(histogram densityH???j, width(4) color(red%30)), ///
			by(laname, title("Histogram of Density H???j and H???j")) ///
			legend( order(1 "Day/densityH???j/green" 2 "Night/densityH???j/red") )
	graph export "$RESULTS/q22b_hist_no3.png", as(png) replace
	
// 	twoway  (histogram H???j, width(4) color(green%30)) ///
// 			(histogram H???j, width(4) color(red%30)), ///
// 			by(laname, title("Histogram of Density H???j and H???j")) ///
// 			legend( order(1 "Day/densityH???j/green" 2 "Night/densityH???j/red") )
// 	graph export "$RESULTS/q22b_hist_no4.png", as(png) replace
	restore
	
// 	preserve 
// 	collapse H???i H???i H???j H???j densityH???i densityH???i densityH???j densityH???j, by(laname)
// 	twoway (kdensity H???i) (kdensity H???i), ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b.png", as(png) replace 
//	
// 	twoway (histogram densityH???i) (histogram densityH???i), ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist.png", as(png) replace
//	
// 	twoway (histogram densityH???j) (histogram densityH???j), by(laname) ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist_no1.png", as(png) replace
//	
// 	restore
	
	*** day is H???i and night is H???i
	preserve
	collapse (sum) H???i (sum) H???i (sum) H???i (sum) H???i, by(laname)
	gen dayNightRatio = H???i/H???i
	save "dayNight.dta", replace
	export excel using "dayNight.xls", sheetreplace firstrow(variables)
	restore

********************************************************************************
*   2.3 . Solve for unobserved transformed wages ??j=wj^?? and obtain wj       *
********************************************************************************

* get data into mata and reshape into matrices and vectors
mata iNfrom = max(st_data(.,"msoai_n"))
mata iNto = max(st_data(.,"msoaj_n"))
mata v = st_numscalar("??")

* get commtime with rows j and cols i
sort msoaj_n msoai_n // make sure to sort here, essential!
mata mCommtime = rowshape(st_data(.,"commtime"), iNfrom)
mata vHri = rowshape(st_data(.,"H???i"), iNfrom)[1,.]'
mata vHmj = rowshape(st_data(.,"H???j"), iNfrom)[.,1]


* set starting values for wage vectors (vW0 = 1 and vW1 = 0 ) 
mata vW0 = J(iNto, 1, 1)
mata vW1 = J(iNto, 1, 0)

* get/set parameters for while loop
mata error = sum(abs(vW0-vW1))
mata thresholdwages = st_numscalar("thresholdwages")
mata j = 1

mata while(error > thresholdwages & j < 100){

	// take weighted sum for each i
	vWSums_i = colsum(J(1,iNto, vW0):/exp(v * mCommtime))

	// take weighted sum for each j,  to get commuting probabilities
	mCommProbs = (J(1,iNto, vW0) :/ exp(v * mCommtime)):/vWSums_i

	// matrix multiply with vHri to get nr of inbound commutes for each j
	vInbound_j =mCommProbs * vHri

	// update the wages
	vW1 = (vHmj:/vInbound_j) :* vW0

	// if missing replace with random
	 if ( missing(vW1) >0 ){
		vW0 = 0.95+(1.05-0.95)*runiform(1,iNto)
		vW1 =0.95+(1.05-0.95)*runiform(1,iNto)
	  }
		
	// update the wage vector
	vW0 = vW1
			
	// normalise by geometric mean (i.e. exp(mean(log(vW0)))
	vW0=vW0/exp(mean(log(vW0)))
		
	error = sum(abs(vW0-vW1))
	j = j + 1
}
	
	* put mata output in stata
	mata vWj = vec(J(1,iNto , vW0)') //repeat for each j
	mata vWi = J(iNto, 1, vW0) //repeat for each i

	mata st_addvar("float","??i")
	mata st_addvar("float","??j")
	mata st_store(.,"??i",vWi)
	mata st_store(.,"??j",vWj)

	* get final update of mCommProbs and export as vector to stata
	mata mCommProbs = (J(1,iNto, vW0) :/ exp(v * mCommtime)):/vWSums_i
	mata vCommProbs = vec(mCommProbs')
	mata st_addvar("float","??i")
	mata st_store(.,"??i",vCommProbs)
	
	// my code 
// 	bys laname: tabstat ??i ??j, stat(mean sd p1 p99) col(stat)
// 	for any ??i ??j: bys laname: tabstat X, stat(mean sd p1 p99) col(stat)
	preserve
	collapse (mean) mean_??i = ??i (mean) mean_??j = ??j (max) max_??i = ??i (max) max_??j = ??j, by(laname)
	egen max??i = max(max_??i)
	egen max??j = max(max_??j)
	gen maxIndicator??i = 1 if max_??i == max??i
	gen maxIndicator??j = 1 if max_??j == max??j
	save "transformedWages.dta", replace
	export excel using "transformedWages.xls", sheetreplace firstrow(variables)
	restore
	
	preserve
	twoway 	(histogram ??j, width(4) color(green%30)) ///
			(histogram ??i, width(4) color(red%30)) ///
			, by(laname, title("Histogram $\omega_j$ and $\omega_i$")) ///
			legend( order(1 "$\omega_j$/Green" 2 "$\omega_i$/Red") )
	graph export "$RESULTS/q23a.png", as(png) replace		
	restore
	
	
********************************************************************************
*    2.4. calulate commuting heterogeneity and ?? *
********************************************************************************
	
	* calculate standard deviations within each local authority 
	qui bysort laname: egen ??log?? = sd(ln(??j)) if msoai==msoaj
	qui bysort laname: egen ??logmedinc = sd(ln(medincj)) if msoai==msoaj
	
	* turn these SDs into variances (stata does not want to do is in one step)
	qui gen ??2log?? = ??log??^2
	qui gen ??2logmedinc = ??logmedinc^2
	drop ??log?? ??logmedinc

	* find (1/??)^2 that minimizes squared difference
	
		// your code here
		eststo clear
		eststo ols01: reg ??2logmedinc ??2log?? , r nocons
		esttab ols01 using "$RESULTS/Q24.tex", se star(* 0.1 ** 0.05 *** 0.01) ///
		title(OLS for Step 3 \label{Q24Tab1}) stats(N r2_a) ///
		label booktabs replace
		
		eststo clear
		
	* calculate ?? 
	
		// your code here
		scalar varepsilon = (_b[??2log??])^(-1/2)
	
	* calculate ?? 
	
		//your code here
		scalar kappa = ??/varepsilon
		

********************************************************************************
*    2.5. Recover amenities and production endowments (structural residuals)     *
********************************************************************************

*normalize prices
qui ameans Qi if msoai==msoaj
qui scalar Qi_gm = r(mean_g)         // Geometric means, see equation (S.47)

*normalize adjusted wages
qui ameans ??i if msoai==msoaj
qui scalar ??i_gm = r(mean_g)         // Geometric means, see equation (S.47)

*normalize residential population
qui ameans H???i if msoai==msoaj
qui scalar H???i_gm = r(mean_g)        // Geometric means, see equation (S.47)

*normalize wages
qui bysort msoai: egen Wi = total(??j/exp(??*commtime))
qui ameans Wi if msoai==msoaj
qui scalar Wi_gm = r(mean_g)         // Geometric means, see equation (S.47)



*calculate productiviy levels (Ai) and amenity levels (Bi)

	//your code here
	*** Use Eq 27 and 28 of ARSW
	gen termAi1 = (1-??)*ln(Qi/Qi_gm)
	gen termAi2 = (??/varepsilon)*ln(??i/??i_gm)
	gen Ai = exp( termAi1 + termAi2 )
	gen termBi1 = (1/varepsilon)*ln(H???i/H???i_gm)
	gen termBi2 = (1-??)*ln(Qi/Qi_gm)
	gen termBi3 = (1/varepsilon)*ln(Wi/Wi_gm)
	gen Bi = exp( termBi1 + termBi2 - termBi3 )
	
	drop termAi1 termAi2 termBi1 termBi2 termBi3
	
	** Table?
	preserve 
	collapse (mean) meanAi = Ai (median) medAi = Ai (max) maxAi = Ai (mean) meanBi = Bi (median) medBi = Bi (max) maxBi = Bi, by(laname)
	
	egen maximumAi = max(maxAi)
	gen indicatorAi = 1 if maxAi == maximumAi
	egen maximumBi = max(maxBi)
	gen indicatorBi = 1 if maxBi == maximumBi
	drop maximumAi maximumBi
	
	egen maxMedAi = max(medAi)
	egen maxMedBi = max(medBi)
	
	gen normMedAi = (medAi/maxMedAi)*100
	gen normMedBi = (medBi/maxMedBi)*100
	
	drop maxMedAi maxMedBi
	
	save "productivityAmenities.dta", replace
	export excel using "productivityAmenities.xls", sheetreplace firstrow(variables)
	restore

	
********************************************************************************
*    2.6 Estimate effect of open space on Amenity level Bi                     *
********************************************************************************

	//your code here
	gen logBi = log(Bi)
	eststo clear
	eststo robust01: reg logBi Oi, r
	eststo cluster01: reg logBi Oi, vce(cluster msoai_n)
	
	esttab robust01 cluster01 using "$RESULTS/Q26.tex", se star(* 0.1 ** 0.05 *** 0.01) ///
		title(OLS for Step 4: Amenities \label{Q26Tab1}) stats(N r2_a) ///
		label booktabs replace
		
	eststo clear



********************************************************************************
*    Closing: export data under a new name for use in Part II
********************************************************************************

	// your code here
	save "TI_UTE_London_data_outP1.dta", replace
