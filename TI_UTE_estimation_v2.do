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
scalar α = 0.80    // Share of labour in firm costs (from Ahlfeldt et al)
scalar β = 0.75    // Share of expenses on other goods (from Ahlfeldt et al)
scalar μ = 0.75    // Share of capital in construction costs (from Ahlfeldt et al)

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
	scalar ν = -_b[commtime] //your code here
	
	
********************************************************************************
*    2.2 Obtain number of people working (Hᴹ) and living (Hᴿ) in each area     *
********************************************************************************

* Residential population
qui total ncomm
qui replace ncomm = ncomm*( H_total /_b[ncomm])   //Total number of jobs in Greater London in 2011

* calculate residents for each location
qui bysort msoai: egen Hᴿi = sum(ncomm) // people living in i
qui label variable Hᴿi "Number of workers living in location i"

qui gen  tempj_resident = Hᴿi if msoai==msoaj // copy Hᴿi when i == j
qui bysort msoaj: egen Hᴿj = mean(tempj_resident) // take mean to copy value to all j in group
qui label variable Hᴿj "Number of workers living in location j"

* jobs
qui bysort msoaj: egen Hᴹj = sum(ncomm) // people working in j
qui label variable Hᴹj "Number of workers working in location j"

qui gen  tempi_worker = Hᴹj if msoai==msoaj // copy Hᴹj when i == j
qui bysort msoai: egen Hᴹi = mean(tempi_worker) //take mean to copy value to all i in group
qui label variable Hᴹi "Number of workers working in location i"

* remove helper columns
qui drop tempi_worker tempj_resident


*make graphs with population data
	
	//your code here
	
// 	** standardize residential population and working population
// 	for any Hᴿi Hᴿj Hᴹi Hᴹj: egen zX = std(X)
//	
// 	** descriptive statistics
// 	tabstat Hᴿi Hᴿj, stat(mean sd median) col(stat)
// 	tabstat Hᴹi Hᴹj, stat(mean sd median) col(stat)
//	
// 	** residential (i vs j)
// 	twoway (kdensity zHᴿi) (kdensity zHᴿj)
//		
// 	** worker (i vs j)
// 	twoway (kdensity zHᴹi) (kdensity zHᴹj)
//	
// 	** origin v. destination
// 	twoway (kdensity zHᴿi) (kdensity zHᴹj)
// 	twoway (kdensity zHᴿj) (kdensity zHᴹi)
	
	** population/area 
	gen densityHᴿi = Hᴿi/landha
	label variable densityHᴿi "Density of workers living in location i"
	gen densityHᴹi = Hᴹi/landha
	label variable densityHᴹi "Density of workers working in location i"
	gen densityHᴿj = Hᴿj/landha
	label variable densityHᴿj "Density of workers living in location j"
	gen densityHᴹj = Hᴹj/landha 
	label variable densityHᴹj "Density of workers working in location j"
	tabstat Hᴿi Hᴹi Hᴿj Hᴹj densityHᴿi densityHᴹi densityHᴿj densityHᴹj, stat(mean median) col(stat)
	
	*** try to make histogram/kdensity of population area
	preserve 
// 	twoway (histogram densityHᴹi) (histogram densityHᴿi), by(laname) ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist_no.png", as(png) replace
	
	twoway  (histogram densityHᴹj, width(4) color(green%30)) ///
			(histogram densityHᴿj, width(4) color(red%30)), ///
			by(laname, title("Histogram of Density Hᴹj and Hᴿj")) ///
			legend( order(1 "Day/densityHᴹj/green" 2 "Night/densityHᴿj/red") )
	graph export "$RESULTS/q22b_hist_no3.png", as(png) replace
	
// 	twoway  (histogram Hᴹj, width(4) color(green%30)) ///
// 			(histogram Hᴿj, width(4) color(red%30)), ///
// 			by(laname, title("Histogram of Density Hᴹj and Hᴿj")) ///
// 			legend( order(1 "Day/densityHᴹj/green" 2 "Night/densityHᴿj/red") )
// 	graph export "$RESULTS/q22b_hist_no4.png", as(png) replace
	restore
	
// 	preserve 
// 	collapse Hᴹi Hᴿi Hᴹj Hᴿj densityHᴹi densityHᴿi densityHᴹj densityHᴿj, by(laname)
// 	twoway (kdensity Hᴹi) (kdensity Hᴿi), ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b.png", as(png) replace 
//	
// 	twoway (histogram densityHᴹi) (histogram densityHᴿi), ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist.png", as(png) replace
//	
// 	twoway (histogram densityHᴹj) (histogram densityHᴿj), by(laname) ///
// 		legend( order(1 "Day" 2 "Night") )
// 	graph export "$RESULTS/q22b_hist_no1.png", as(png) replace
//	
// 	restore
	
	*** day is Hᴹi and night is Hᴿi
	preserve
	collapse (sum) Hᴹi (sum) Hᴿi (sum) Hᴹi (sum) Hᴿi, by(laname)
	gen dayNightRatio = Hᴹi/Hᴿi
	save "dayNight.dta", replace
	export excel using "dayNight.xls", sheetreplace firstrow(variables)
	restore

********************************************************************************
*   2.3 . Solve for unobserved transformed wages ωj=wj^ε and obtain wj       *
********************************************************************************

* get data into mata and reshape into matrices and vectors
mata iNfrom = max(st_data(.,"msoai_n"))
mata iNto = max(st_data(.,"msoaj_n"))
mata v = st_numscalar("ν")

* get commtime with rows j and cols i
sort msoaj_n msoai_n // make sure to sort here, essential!
mata mCommtime = rowshape(st_data(.,"commtime"), iNfrom)
mata vHri = rowshape(st_data(.,"Hᴿi"), iNfrom)[1,.]'
mata vHmj = rowshape(st_data(.,"Hᴹj"), iNfrom)[.,1]


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

	mata st_addvar("float","ωi")
	mata st_addvar("float","ωj")
	mata st_store(.,"ωi",vWi)
	mata st_store(.,"ωj",vWj)

	* get final update of mCommProbs and export as vector to stata
	mata mCommProbs = (J(1,iNto, vW0) :/ exp(v * mCommtime)):/vWSums_i
	mata vCommProbs = vec(mCommProbs')
	mata st_addvar("float","πi")
	mata st_store(.,"πi",vCommProbs)
	
	// my code 
// 	bys laname: tabstat ωi ωj, stat(mean sd p1 p99) col(stat)
// 	for any ωi ωj: bys laname: tabstat X, stat(mean sd p1 p99) col(stat)
	preserve
	collapse (mean) mean_ωi = ωi (mean) mean_ωj = ωj (max) max_ωi = ωi (max) max_ωj = ωj, by(laname)
	egen maxωi = max(max_ωi)
	egen maxωj = max(max_ωj)
	gen maxIndicatorωi = 1 if max_ωi == maxωi
	gen maxIndicatorωj = 1 if max_ωj == maxωj
	save "transformedWages.dta", replace
	export excel using "transformedWages.xls", sheetreplace firstrow(variables)
	restore
	
	preserve
	twoway 	(histogram ωj, width(4) color(green%30)) ///
			(histogram ωi, width(4) color(red%30)) ///
			, by(laname, title("Histogram $\omega_j$ and $\omega_i$")) ///
			legend( order(1 "$\omega_j$/Green" 2 "$\omega_i$/Red") )
	graph export "$RESULTS/q23a.png", as(png) replace		
	restore
	
	
********************************************************************************
*    2.4. calulate commuting heterogeneity and κ *
********************************************************************************
	
	* calculate standard deviations within each local authority 
	qui bysort laname: egen σlogω = sd(ln(ωj)) if msoai==msoaj
	qui bysort laname: egen σlogmedinc = sd(ln(medincj)) if msoai==msoaj
	
	* turn these SDs into variances (stata does not want to do is in one step)
	qui gen σ2logω = σlogω^2
	qui gen σ2logmedinc = σlogmedinc^2
	drop σlogω σlogmedinc

	* find (1/ε)^2 that minimizes squared difference
	
		// your code here
		eststo clear
		eststo ols01: reg σ2logmedinc σ2logω , r nocons
		esttab ols01 using "$RESULTS/Q24.tex", se star(* 0.1 ** 0.05 *** 0.01) ///
		title(OLS for Step 3 \label{Q24Tab1}) stats(N r2_a) ///
		label booktabs replace
		
		eststo clear
		
	* calculate ε 
	
		// your code here
		scalar varepsilon = (_b[σ2logω])^(-1/2)
	
	* calculate κ 
	
		//your code here
		scalar kappa = ν/varepsilon
		

********************************************************************************
*    2.5. Recover amenities and production endowments (structural residuals)     *
********************************************************************************

*normalize prices
qui ameans Qi if msoai==msoaj
qui scalar Qi_gm = r(mean_g)         // Geometric means, see equation (S.47)

*normalize adjusted wages
qui ameans ωi if msoai==msoaj
qui scalar ωi_gm = r(mean_g)         // Geometric means, see equation (S.47)

*normalize residential population
qui ameans Hᴿi if msoai==msoaj
qui scalar Hᴿi_gm = r(mean_g)        // Geometric means, see equation (S.47)

*normalize wages
qui bysort msoai: egen Wi = total(ωj/exp(ν*commtime))
qui ameans Wi if msoai==msoaj
qui scalar Wi_gm = r(mean_g)         // Geometric means, see equation (S.47)



*calculate productiviy levels (Ai) and amenity levels (Bi)

	//your code here
	*** Use Eq 27 and 28 of ARSW
	gen termAi1 = (1-α)*ln(Qi/Qi_gm)
	gen termAi2 = (α/varepsilon)*ln(ωi/ωi_gm)
	gen Ai = exp( termAi1 + termAi2 )
	gen termBi1 = (1/varepsilon)*ln(Hᴿi/Hᴿi_gm)
	gen termBi2 = (1-β)*ln(Qi/Qi_gm)
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
