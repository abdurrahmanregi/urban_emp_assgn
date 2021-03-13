********************************************************************************
*                                                                              *
*                                  TI UTE Assignment                           *
*                              Part II - Counterfactuals                       *
********************************************************************************

* Set working directory
	
	//your code here
// 	cd "~/Dropbox/1educ/0_TI_Urban_Econ/TI_UTE_emp_assign"
// 	use "TI_UTE_London_data_results.dta", replace
	
	cd "C:\Users\Regi Kusumaatmadja\Dropbox\TI year 2\Block 03\Urban code"
	use "TI_UTE_London_data_outP1.dta", clear


********************************************************************************
*                      3. Set  parameters                            *
********************************************************************************
scalar α = 0.80    // Share of labour in firm costs (from Ahlfeldt et al)
scalar β = 0.75    // Share of expenses on other goods (from Ahlfeldt et al)
scalar μ = 0.75    // Share of capital in construction costs (from Ahlfeldt et al)
scalar H_total = 4864000 // Total working population in London
scalar ς = 0.25 // weighted average for updating, see (S.83)


* parameters obtained in part I (replace with your own)
// scalar ε = 3
// scalar κ = 0.03
// scalar ν = 0.1

scalar ε = 3.3263865
scalar κ = 0.03020185
scalar ν = 0.10046303


*define estimates for counterfactual, obtained in question 2.6 

	
	//your code here
	reg logBi Oi, vce(cluster msoai_n)
	scalar alphaQ26 = _b[_cons]
	scalar betaQ26 = _b[Oi]

********************************************************************************
*                 3.1. Create counterfactual Bi and floorsplace use           *
********************************************************************************

* create counterfactual open space: Oi_c

	
	// your code here
	*** 50% reduction of share of open area
	gen Oi_c = 0.5*Oi 

* create counterfactual amenety levels: Bi_c 
	
	// your code here
	*** London does not have open space. That is, Oi = 0
// 	gen Bi_c = exp(alphaQ26 + betaQ26*Oi_c) if Oi > 0
// 	replace Bi_c = Bi if Oi == 0
	
	gen Bi_c = Bi*exp(betaQ26*(Oi_c-Oi))
	
* create counterfactual land area: Ki_c
	
	// your code here
// 	gen Ki_c = landha + (Oi_c - Oi)*openspaceha
	gen Ki_c = landha + (Oi_c)*areaha
	

* scatter plots to compare initial and counterfactual cases

	// your code here
	*** For amenities
	twoway (scatter Bi_c Bi) (lfit Bi_c Bi), title("Amenities: Counterfactuals v. Levels: All")
	graph export "$RESULTS/amenitiesCounterfactual.png", as(png) replace
	
	twoway (scatter Bi_c Bi) (lfit Bi_c Bi), by(laname, title("Amenities: Counterfactuals v. Levels: By District"))
	graph export "$RESULTS/amenitiesCount_laname.png", as(png) replace	
		
	*** For land
	twoway (scatter Ki_c Ki) (lfit Ki_c Ki), title("Land Area: Counterfactuals v. Levels: All")
	graph export "$RESULTS/landCounterfactual.png", as(png) replace
	
	twoway (scatter Ki_c Ki) (lfit Ki_c Ki), by(laname, title("Land Area: Counterfactuals v. Levels: By District"))
	graph export "$RESULTS/landCounterfactual_laname.png", as(png) replace
	

* !! NB: the algorithm below assumes counterfactuals Bi_c and Ki_c exist !!

********************************************************************************
*                        3.2 Set/update starting values                         *
********************************************************************************

*calculate land use
qui g Lᴹi = (((ωi^(1/ε))/(α*Ai))^(1/(1-α)))*Hᴹi
qui bysort msoai: egen sumi = total(πi*(ωj^(1/ε))/exp(κ*commtime))
qui g Lᴿi = (1-β)*sumi*(Hᴿi/Qi)
qui g θi = Lᴹi/(Lᴹi+Lᴿi)
qui g φi = (Lᴹi + Lᴿi)/Ki_c^(1-μ)	

* parameters for while loop
local threshold = 1
gen error = .
local error = 100000
local r = 1

qui while `error' > `threshold' & `r'<25 {
	
	*for first run get starting values from Part I
	qui if `r'  == 1 {
			g wi_c = ωi^(1/ε)
			g wj_c = ωj^(1/ε)
			g Qi_c = Qi
			g θi_c = θi
			g Lᴹi_c = Lᴹi
			g Lᴿi_c = Lᴿi
			g Hᴿi_c = Hᴿi
			g Hᴹi_c = Hᴹi
			g Hᴹj_c = Hᴹj
			g Hᴿj_c = Hᴿj
	}

	*for subsequent runs get updated data
	qui if `r'  > 1 {
			* update variables
			replace wi_c = ς*wi_c + (1-ς)*wi_c1
			replace Qi_c = ς*Qi_c + (1-ς)*Qi_c1
			replace θi_c = ς*θi_c + (1-ς)*θi_c1
			
			* transfer wages from i to j
			drop wj_c
			g temp = wi_c if msoai==msoaj
			bysort msoaj: egen wj_c = mean(temp)
			drop temp
			
			* drop vars that will be updates
			drop wi_c1 wj_c1 Qiᴹ_c1 Qiᴿ_c1 Qi_c1 Yi_c1 θi_c1
			drop Hᴿi_c1 Hᴹj_c1 Hᴿj_c1 Hᴹi_c1
			drop dij πij_c1 πiji_c1 vi_c1   
	}

	
********************************************************************************
*                       Update quilibrium equations                            *
********************************************************************************
	* NB: here we use the equations as in ASRW supplement page 56-57

	* equation (S.71)
	gen dij =  exp(κ * commtime)
	gen prob = ((dij * Qi_c^(1-β))^-ε) * (Bi_c * wj_c)^ε
	egen sumij = total(prob)
	g πij_c1 = prob/sumij
	drop prob sum*
	
	* equation (S.72)
	qui bysort msoai: egen sumi = total((wj_c/dij)^ε)
	qui g πiji_c1 = (wj_c / dij)^ε / sumi
	drop sumi
	
	
	* equation (S.73)
	bysort msoai: egen Hᴿi_c1 = total(πij_c1 * H_total)
	g temp = Hᴿi_c if msoai==msoaj
	bysort msoaj: egen Hᴿj_c1 = mean(temp)
	drop temp
		
		
	* equation (S.74)
	bysort msoaj: egen Hᴹj_c1 = total(πij_c1 * H_total)
	g temp = Hᴹj_c1 if msoaj==msoai
	bysort msoai: egen Hᴹi_c1 = mean(temp)
	drop temp

	* equation (S.75)
	gen Yi_c1 = Ai * ((Hᴹi_c1)^α) * (θi_c * φi * (Ki_c^(1-μ)))^(1-α)
	
	* equation (S.76)
	gen wi_c1 = (α * Yi_c1)/(Hᴹi_c1)
		
	* also apply to j
	g temp = wi_c1 if msoai==msoaj
	bysort msoaj: egen wj_c1 = mean(temp)
	drop temp
	
	* equation (S.77)
	bysort msoai: egen vi_c1 = total(πiji_c1 * (wj_c) / dij)
	
	* equation (S.78) - commercial floorprice 
	gen Qiᴹ_c1 = ((1 - α)*Yi_c1)/(θi_c * φi * Ki_c^(1-μ))
	
	* equation (S.79)  - residential floorprice
	gen Qiᴿ_c1 = ((1 - β) * vi_c1 * Hᴿi_c1) / ((1-θi_c) * φi * (Ki_c^(1 - μ)))
	 
	* combine into one floorprice
	gen Qi_c1 = max(Qiᴹ_c1, Qiᴿ_c1)
 	
		
	* equation (S.80-82) 
	gen θi_c1 = .
	replace θi_c1 = 1  if Hᴿi_c1 == 0  // only commercial
	replace θi_c1 = 0  if Hᴹi_c1 == 0  // only residential
	replace θi_c1 = ((1 - α)*Yi_c1)/(Qi_c1 * φi * Ki_c^(1 - μ))
	
	
	* re-normalise wages with geometric means
	ameans wi_c1 if msoai==msoaj
	replace wi_c1 = wi_c1/r(mean_g)	

	* update and print error 
	replace error = abs(Qi_c1 - Qi_c) + abs(wi_c1 - wi_c) + abs(θi_c1 - θi_c) 
	noisily total error if msoai==msoaj
	qui local error = _b[error]


	* update loop counter
	qui local r = `r'+1
	
}



* produce results

	gen Oi_diff = Oi_c - Oi
	gen Qi_diff = Qi_c - Qi
	scatter Qi_diff Oi_diff if msoai==msoaj
	
	
	// your code here
	*** what should we do here?
	*** Part (a)
	twoway (scatter Qi_diff Oi_diff if msoai==msoaj) (lfit Qi_diff Oi_diff), title("Change in Open Space i on the House Prices i: All") 
	graph export "$RESULTS/Q32a.png", as(png) replace
	 
	*** Part (a) by district
	twoway (scatter Qi_diff Oi_diff if msoai==msoaj) (lfit Qi_diff Oi_diff), by(laname, title("Change in Open Space i on the House Prices i: By District"))
	graph export "$RESULTS/Q32a_laname.png", as(png) replace
	
	*** Part (b): Production and Open Space
	*** First, generate production at levels (initial)
	gen Yi = Ai * ((Hᴹi)^α) * (θi * φi * (Ki^(1-μ)))^(1-α)
// 	gen Yi = Ai * ((Hᴹi)^α) * ( Lᴹi )^(1-α)
	gen Yi_diff = Yi_c1 - Yi 
	
	*** Second, create scatter 
	twoway (scatter Yi_diff Oi_diff if msoai==msoaj) (lfit Yi_diff Oi_diff), title("Change in Open Space i on the Production i: All")
	graph export "$RESULTS/Q32b.png", as(png) replace
	
	twoway (scatter Yi_diff Oi_diff if msoai==msoaj) (lfit Yi_diff Oi_diff), by(laname, title("Change in Open Space i on the Production i: By District"))
	graph export "$RESULTS/Q32b_laname.png", as(png) replace
	
	*** Part (c)
// 	NOT YET
	
	*** Part (d)
	*** Define change in workplace
	gen Hᴹi_diff = Hᴹi_c1 - Hᴹi
	gen Hᴹj_diff = Hᴹj_c1 - Hᴹj
	
	twoway (scatter Hᴹi_diff Oi_diff if msoai==msoaj) (lfit Hᴹi_diff Oi_diff), title("Change in Open Space i on the Worker's Population i: All")
	graph export "$RESULTS/Q32d_01.png", as(png) replace
	
	twoway (scatter Hᴹj_diff Oi_diff if msoai==msoaj) (lfit Hᴹj_diff Oi_diff), title("Change in Open Space i on the Worker's Population j: All")
	graph export "$RESULTS/Q32d_02.png", as(png) replace
	
	twoway (scatter Hᴹi_diff Oi_diff if msoai==msoaj) (lfit Hᴹi_diff Oi_diff), by(laname, title("Change in Open Space i on the Worker's Population i: By District"))
	graph export "$RESULTS/Q32d_01_laname.png", as(png) replace
	
	twoway (scatter Hᴹj_diff Oi_diff if msoai==msoaj) (lfit Hᴹj_diff Oi_diff), by(laname, title("Change in Open Space i on the Worker's Population j: By District"))
	graph export "$RESULTS/Q32d_02_laname.png", as(png) replace
	
	*** Part (e)
	gen Hᴿi_diff = Hᴿi_c1 - Hᴿi
	gen Hᴿj_diff = Hᴿj_c1 - Hᴿj

	twoway (scatter Hᴿi_diff Oi_diff if msoai==msoaj) (lfit Hᴿi_diff Oi_diff), title("Change in Open Space i on the Residential's Population i: All")
	graph export "$RESULTS/Q32e_01.png", as(png) replace
	
	twoway (scatter Hᴿj_diff Oi_diff if msoai==msoaj) (lfit Hᴿj_diff Oi_diff), title("Change in Open Space i on the Residential's Population j: All")
	graph export "$RESULTS/Q32e_02.png", as(png) replace
	
	twoway (scatter Hᴿi_diff Oi_diff if msoai==msoaj) (lfit Hᴿi_diff Oi_diff), by(laname, title("Change in Open Space i on the Residential's Population i: By District"))
	graph export "$RESULTS/Q32e_01_laname.png", as(png) replace
	
	twoway (scatter Hᴿj_diff Oi_diff if msoai==msoaj) (lfit Hᴿj_diff Oi_diff), by(laname, title("Change in Open Space i on the Residential's Population j: By District"))
	graph export "$RESULTS/Q32e_02_laname.png", as(png) replace
	
	*** Part (f)
	gen Bi_diff = Bi_c - Bi
	gen Ki_diff = Ki_c - Ki 
	gen θi_diff = θi_c1 - θi
	
	****** Amenities v Open space
	twoway (scatter Bi_diff Oi_diff if msoai==msoaj) (lfit Bi_diff Oi_diff), title("Change in Open Space i on Amenities i: All")
	graph export "$RESULTS/Q32f_01.png", as(png) replace
	
	****** Land area v Open space
	twoway (scatter Ki_diff Oi_diff if msoai==msoaj) (lfit Ki_diff Oi_diff), title("Change in Open Space i on Land Area i: All")
	graph export "$RESULTS/Q32f_02.png", as(png) replace
	
	****** Share of commercial, thetai, v Open space
	twoway (scatter θi_diff Oi_diff if msoai==msoaj) (lfit θi_diff Oi_diff), title("Change in Open Space i on Share of Commercial i: All")
	graph export "$RESULTS/Q32f_03.png", as(png) replace
	
********************************************************************************
*                       3.3 Calculate utility levels
********************************************************************************

* utility levels in observed equilibrium
	
	// your code here
	gen uij = Bi*ωj*exp(-kappa*commtime)*Qi^(β-1)
	
* utility levels in counterfactual equilibrium
	
	// your code here
	gen uij_c = Bi_c*wj_c*exp(-kappa*commtime)*Qi_c^(β-1)
	
* produce output for report
	
	// your code here	
	
	***TBA
	
	save "TI_UTE_London_data_outP2.dta", replace