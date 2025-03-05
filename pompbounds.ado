/*This is the 'pobounds' subcommand, in which the non-MIV bounds on mean potential outcomes are computed,
as well as the CLR-type bounding functions for the bounds with MIV*/

program pompbounds, eclass properties(svyb)
	version 10.0
	#delimit ;
	syntax varlist (min=2 max=2 fv) [if] [in] [pw iw],  YMIN(string) YMAX(string) [ATT NMTS NMTR MIV(varname) bins(integer 5) DISCRETEMIV NMIV] UL(string) CONTRT(string) CONTRS(string) asmptn(string) ;
	marksample touse;
	#delimit cr
	
	/*Gather variables/options for use in the program*/

	gettoken y T : varlist
	
	tokenize "`ymin'"
	tokenize "`ymax'"
	tokenize "`bins'"
	tokenize "`ul'"
	tokenize "`contrt'"
	tokenize "`asmptn'"
	
	/*Preserve the sample in memory so nothing is changed for the user*/
	
	preserve
	
	/*Incorporate any restrictions from the user in addition to the enforcement
	of bounded support via ymin and ymax (there was an issue in testing so we
	subtract from/add to the ymin and ymax to ensure nothing dropped in error*/
	
	qui keep if `touse'
	qui keep if inrange(`y',(`ymin'-0.0000001),(`ymax'+0.0000001))

	/*Store levels of the treatment and then the min/max*/
	
	qui levelsof `T' if `touse', local(Tlevels)
	
	qui sum `T'
	local minT = r(min)
	local maxT = r(max)
	
	/*PRELIMINARY COMPUTATIONS FOR THE BOUNDS*/

	/*Create indicators for each treatment level*/
	tempvar i_treat itreat
	qui ta `T', gen(`i_treat')
	
	local i=1
	foreach t of local Tlevels {
		gen `itreat'`t' = `i_treat'`i'
		local ++i
	}
	
	tempname eb

	/*Estimate Pr[T=t] using the indicators above*/
	foreach t of local Tlevels  {
		qui reg `itreat'`t' [`weight' `exp']
		qui mat `eb'=e(b)
		global pr_`t' = `eb'[1,1]
	}
	
	/*Estimate Pr[T>t], Pr[T<t], Pr[T>=t], and Pr[T<=t]*/
	foreach t of local Tlevels {
	    tempvar Tgt`t' Tlt`t' Tgeq`t' Tleq`t'
		qui gen `Tgt`t'' = (`T'>`t')
		qui gen `Tlt`t'' = (`T'<`t')
		qui gen `Tgeq`t'' = (`T'>=`t')
		qui gen `Tleq`t'' = (`T'<=`t')
	}
	
	foreach t of local Tlevels  {
		qui reg `Tgt`t'' [`weight' `exp']
		qui mat `eb'=e(b)
		global prTgt`t' = `eb'[1,1]
		
		qui reg `Tlt`t'' [`weight' `exp']
		qui mat `eb'=e(b)
		global prTlt`t' = `eb'[1,1]
		
		qui reg `Tgeq`t'' [`weight' `exp']
		qui mat `eb'=e(b)
		global prTgeq`t' = `eb'[1,1]
		
		qui reg `Tleq`t'' [`weight' `exp']
		qui mat `eb'=e(b)
		global prTleq`t' = `eb'[1,1]
	}

	/*Estimate conditional means 
	E[y | T=t], E[y | T<t], E[y | T>t], E[y | T<=t], and E[y | T>=t]*/
	
		foreach t of local Tlevels  {
			qui reg `y' if `T'==`t' [`weight' `exp']
			qui mat `eb'=e(b)
			global m_y`t' = `eb'[1,1]
			
			/*Need to use capture for when T<minT ad T>maxT... replace with zero
			when given the `no observations' error*/

			cap qui reg `y' if `T' < `t' [`weight' `exp']
				if _rc==2000 {
					global m_yL`t' = 0
				}
				else {
					qui mat `eb'=e(b)
					global m_yL`t' = `eb'[1,1]
				}

			cap qui reg `y' if `T' > `t' [`weight' `exp']
				if _rc==2000 {
					global m_yG`t' = 0
				}
				else {
					qui mat `eb'=e(b)
					global m_yG`t' = `eb'[1,1]
				}
		
			qui reg `y' if `T' <= `t' [`weight' `exp']
			qui mat `eb'=e(b)
			global m_yLeq`t' = `eb'[1,1]

			qui reg `y' if `T' >= `t' [`weight' `exp']
			qui mat `eb'=e(b)
			global m_yGeq`t' = `eb'[1,1]
		}
		
		/*Ensure the means for T less than min treatment and T greater than max treatment
		are set to zero (useful below)*/
		
		global m_yL`minT' = 0
		global m_yG`maxT' = 0


	****************************************************************************
	****************************************************************************
	/*BOUNDING TREATMENT EFFECTS*/
	
	/*BOUNDING ATT(t,s), IF REQUESTED*/
	if "`att'" != "" {
		local subminT = `minT' + 1
		local submaxT = `maxT' - 1
		/*Bounds on the ATT, no MIV*/

		/*Non-positive MTS + non-negative MTR*/
		if ("`nmts'" != "" & "`nmtr'" == "") {
			forval r=`subminT'/`maxT' {
				local w=`r'-1
				forval s=`minT'/`w' {
			    if "`asmptn'" == "wc" {
					qui scalar LB_`s'_`r'_wc = `ymin'
					qui scalar UB_`s'_`r'_wc = `ymax'
					qui scalar LB_`r'_`r'_wc = ${m_y`r'}
					qui scalar UB_`r'_`r'_wc = ${m_y`r'}
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`s'_`r'_mts = `ymin'
					qui scalar UB_`s'_`r'_mts = ${m_y`s'}
					qui scalar LB_`r'_`r'_mts = ${m_y`r'}
					qui scalar UB_`r'_`r'_mts = ${m_y`r'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`s'_`r'_mtr = `ymin'
					qui scalar UB_`s'_`r'_mtr = ${m_y`r'}
					qui scalar LB_`r'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtr = ${m_y`r'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`s'_`r'_mtsr = `ymin'
					qui scalar UB_`s'_`r'_mtsr = min(${m_y`r'},${m_y`s'})
					qui scalar LB_`r'_`r'_mtsr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtsr = ${m_y`r'}
				}
				}
			}
		}
		
		/*Non-negative MTS + non-positive MTR*/
		
		else if ("`nmts'" == "" & "`nmtr'" != "") {
			forval r=`subminT'/`maxT' {
				local w=`r'-1
				forval s=`minT'/`w' {
			    if "`asmptn'" == "wc" {
					qui scalar LB_`s'_`r'_wc = `ymin'
					qui scalar UB_`s'_`r'_wc = `ymax'
					qui scalar LB_`r'_`r'_wc = ${m_y`r'}
					qui scalar UB_`r'_`r'_wc = ${m_y`r'}
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`s'_`r'_mts = ${m_y`s'}
					qui scalar UB_`s'_`r'_mts = `ymax'
					qui scalar LB_`r'_`r'_mts = ${m_y`r'}
					qui scalar UB_`r'_`r'_mts = ${m_y`r'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`s'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`s'_`r'_mtr = `ymax'
					qui scalar LB_`r'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtr = ${m_y`r'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`s'_`r'_mtsr = max(${m_y`r'},${m_y`s'})
					qui scalar UB_`s'_`r'_mtsr = `ymax'
					qui scalar LB_`r'_`r'_mtsr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtsr = ${m_y`r'}
				}
				}
			}
		}
		
		/*Non-positive MTS + non-positive MTR*/
		
		else if ("`nmts'" != "" & "`nmtr'" != "") {
			forval r=`subminT'/`maxT' {
				local w=`r'-1
				forval s=`minT'/`w' {		
			    if "`asmptn'" == "wc" {
					qui scalar LB_`s'_`r'_wc = `ymin'
					qui scalar UB_`s'_`r'_wc = `ymax'
					qui scalar LB_`r'_`r'_wc = ${m_y`r'}
					qui scalar UB_`r'_`r'_wc = ${m_y`r'}
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`s'_`r'_mts = `ymin'
					qui scalar UB_`s'_`r'_mts = ${m_y`s'}
					qui scalar LB_`r'_`r'_mts = ${m_y`r'}
					qui scalar UB_`r'_`r'_mts = ${m_y`r'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`s'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`s'_`r'_mtr = `ymax'
					qui scalar LB_`r'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtr = ${m_y`r'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`s'_`r'_mtsr = ${m_y`r'}
					qui scalar UB_`s'_`r'_mtsr = ${m_y`s'}
					qui scalar LB_`r'_`r'_mtsr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtsr = ${m_y`r'}
				}
				}
			}
		}
		
		/*Non-negative MTS + non-negative MTR*/
		
		else {
			forval r=`subminT'/`maxT' {
				local w=`r'-1
				forval s=`minT'/`w' {	
			    if "`asmptn'" == "wc" {
					qui scalar LB_`s'_`r'_wc = `ymin'
					qui scalar UB_`s'_`r'_wc = `ymax'
					qui scalar LB_`r'_`r'_wc = ${m_y`r'}
					qui scalar UB_`r'_`r'_wc = ${m_y`r'}
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`s'_`r'_mts = ${m_y`s'}
					qui scalar UB_`s'_`r'_mts = `ymax'
					qui scalar LB_`r'_`r'_mts = ${m_y`r'}
					qui scalar UB_`r'_`r'_mts = ${m_y`r'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`s'_`r'_mtr = `ymin'
					qui scalar UB_`s'_`r'_mtr = ${m_y`r'}
					qui scalar LB_`r'_`r'_mtr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtr = ${m_y`r'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`s'_`r'_mtsr = ${m_y`s'}
					qui scalar UB_`s'_`r'_mtsr = ${m_y`r'}
					qui scalar LB_`r'_`r'_mtsr = ${m_y`r'}
					qui scalar UB_`r'_`r'_mtsr = ${m_y`r'}
				}
				}
			}
		}
		
		/*Bounds on the mean potential outcomes (ATT version), no MIV*/
		
		forval t=`subminT'/`maxT' {
		    forval k=`minT'/`t' {
			if "`asmptn'" == "wc" {
				qui matrix LB_`k'_`t'_wc = LB_`k'_`t'_wc
				qui matrix UB_`k'_`t'_wc = UB_`k'_`t'_wc
			}
			else if "`asmptn'" == "mts" {
				qui matrix LB_`k'_`t'_mts = LB_`k'_`t'_mts
				qui matrix UB_`k'_`t'_mts = UB_`k'_`t'_mts
			}
			else if "`asmptn'" == "mtr" {
				qui matrix LB_`k'_`t'_mtr = LB_`k'_`t'_mtr
				qui matrix UB_`k'_`t'_mtr = UB_`k'_`t'_mtr
			}
			else if "`asmptn'" == "mtsr" {
				qui matrix LB_`k'_`t'_mtsr = LB_`k'_`t'_mtsr
				qui matrix UB_`k'_`t'_mtsr = UB_`k'_`t'_mtsr
			}
			}
			
		}
	restore 
		
	
	}
	
	
	
	
	
	/*BOUNDING mean potential outcomes  (ATE version)*/
	
	local subminT = `minT' + 1
	
	else {
	/*First, estimate all bounds on mean potential outcomes, not involving the MIV*/

		/*Non-positive MTS + non-negative MTR*/
		
		if ("`nmts'" != "" & "`nmtr'" == "") {
			foreach t of local Tlevels  {
			    if "`asmptn'" == "wc" {
					qui scalar LB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymin'*(1-${pr_`t'})
					qui scalar UB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymax'*(1-${pr_`t'})
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`t'_mts = `ymin'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
					qui scalar UB_`t'_mts = `ymax'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`t'_mtr = ${m_yLeq`t'}*${prTleq`t'} + `ymin'*${prTgt`t'}
					qui scalar UB_`t'_mtr = ${m_yGeq`t'}*${prTgeq`t'} + `ymax'*${prTlt`t'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`t'_mtsr = `ymin'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
					qui scalar UB_`t'_mtsr = `ymax'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
				}
			}
		}
		
		/*Non-negative MTS + non-positive MTR*/
		
		else if ("`nmts'" == "" & "`nmtr'" != "") {
			foreach t of local Tlevels  {
			    if "`asmptn'" == "wc" {
					qui scalar LB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymin'*(1-${pr_`t'})
					qui scalar UB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymax'*(1-${pr_`t'})
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`t'_mts = `ymin'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
					qui scalar UB_`t'_mts = `ymax'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`t'_mtr = ${m_yGeq`t'}*${prTgeq`t'} + `ymin'*${prTlt`t'}
					qui scalar UB_`t'_mtr = ${m_yLeq`t'}*${prTleq`t'} + `ymax'*${prTgt`t'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`t'_mtsr = `ymin'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
					qui scalar UB_`t'_mtsr = `ymax'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
				}
			}
		}
		
		/*Non-positive MTS + non-positive MTR*/
		
		else if ("`nmts'" != "" & "`nmtr'" != "") {
			foreach t of local Tlevels  {		
			    if "`asmptn'" == "wc" {
					qui scalar LB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymin'*(1-${pr_`t'})
					qui scalar UB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymax'*(1-${pr_`t'})
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`t'_mts = `ymin'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
					qui scalar UB_`t'_mts = `ymax'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`t'_mtr = ${m_yGeq`t'}*${prTgeq`t'} + `ymin'*${prTlt`t'}
					qui scalar UB_`t'_mtr = ${m_yLeq`t'}*${prTleq`t'} + `ymax'*${prTgt`t'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`t'_mtsr = ${m_yG`t'}*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
					qui scalar UB_`t'_mtsr = ${m_yL`t'}*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
				}
			}
		}
		
		/*Non-negative MTS + non-negative MTR*/
		
		else {
			foreach t of local Tlevels  {	
			    if "`asmptn'" == "wc" {
					qui scalar LB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymin'*(1-${pr_`t'})
					qui scalar UB_`t'_wc = ${m_y`t'}*${pr_`t'} + `ymax'*(1-${pr_`t'})
				}
				else if "`asmptn'" == "mts" {
					qui scalar LB_`t'_mts = `ymin'*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
					qui scalar UB_`t'_mts = `ymax'*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
				}
				else if "`asmptn'" == "mtr" {
					qui scalar LB_`t'_mtr = ${m_yLeq`t'}*${prTleq`t'} + `ymin'*${prTgt`t'}
					qui scalar UB_`t'_mtr = ${m_yGeq`t'}*${prTgeq`t'} + `ymax'*${prTlt`t'}
				}
				else if "`asmptn'" == "mtsr" {
					qui scalar LB_`t'_mtsr = ${m_yL`t'}*${prTlt`t'} + ${m_y`t'}*${prTgeq`t'}
					qui scalar UB_`t'_mtsr = ${m_yG`t'}*${prTgt`t'} + ${m_y`t'}*${prTleq`t'}
				}
			}
		}
		
		/*BOUNDS ON THE MEAN POTENTIAL OUTCOMES, SEPARATELY FROM ATE, NO MIV*/
			/*Non-positive MTS + non-negative MTR*/

			if ("`nmts'" != "" & "`nmtr'" == "") {
				foreach t of local Tlevels  {
					if "`asmptn'" == "wc" {
						qui matrix LB_`t'_wc = LB_`t'_wc
						qui matrix UB_`t'_wc = UB_`t'_wc
					}
					else if "`asmptn'" == "mts" {
						qui matrix LB_`t'_mts = LB_`t'_mts
						qui matrix UB_`t'_mts = UB_`t'_mts
					}
					else if "`asmptn'" == "mtr" {
						qui matrix LB_`t'_mtr = LB_`t'_mtr
						qui matrix UB_`t'_mtr = UB_`t'_mtr
					}
					else if "`asmptn'" == "mtsr" {
						qui matrix LB_`t'_mtsr = LB_`t'_mtsr
						qui matrix UB_`t'_mtsr = UB_`t'_mtsr
					}
				}
			}

			/*Non-negative MTS + non-positive MTR*/

			else if ("`nmts'" == "" & "`nmtr'" != "") {
				foreach t of local Tlevels  {
					if "`asmptn'" == "wc" {
						qui matrix LB_`t'_wc = LB_`t'_wc
						qui matrix UB_`t'_wc = UB_`t'_wc
					}
					else if "`asmptn'" == "mts" {
						qui matrix LB_`t'_mts = LB_`t'_mts
						qui matrix UB_`t'_mts = UB_`t'_mts
					}
					else if "`asmptn'" == "mtr" {
						qui matrix LB_`t'_mtr = LB_`t'_mtr
						qui matrix UB_`t'_mtr = UB_`t'_mtr
					}
					else if "`asmptn'" == "mtsr" {
						qui matrix LB_`t'_mtsr = LB_`t'_mtsr
						qui matrix UB_`t'_mtsr = UB_`t'_mtsr
					}
				}
			}

			/*Non-positive MTS + non-positive MTR*/

			else if ("`nmts'" != "" & "`nmtr'" != "") {
				foreach t of local Tlevels  {		
					if "`asmptn'" == "wc" {
						qui matrix LB_`t'_wc = LB_`t'_wc
						qui matrix UB_`t'_wc = UB_`t'_wc
					}
					else if "`asmptn'" == "mts" {
						qui matrix LB_`t'_mts = LB_`t'_mts
						qui matrix UB_`t'_mts = UB_`t'_mts
					}
					else if "`asmptn'" == "mtr" {
						qui matrix LB_`t'_mtr = LB_`t'_mtr
						qui matrix UB_`t'_mtr = UB_`t'_mtr
					}
					else if "`asmptn'" == "mtsr" {
						qui matrix LB_`t'_mtsr = LB_`t'_mtsr
						qui matrix UB_`t'_mtsr = UB_`t'_mtsr
					}
				}
			}

			/*Non-negative MTS + non-negative MTR*/

			else {
				foreach t of local Tlevels  {	
					if "`asmptn'" == "wc" {
						qui matrix LB_`t'_wc = LB_`t'_wc
						qui matrix UB_`t'_wc = UB_`t'_wc
					}
					else if "`asmptn'" == "mts" {
						qui matrix LB_`t'_mts = LB_`t'_mts
						qui matrix UB_`t'_mts = UB_`t'_mts
					}
					else if "`asmptn'" == "mtr" {
						qui matrix LB_`t'_mtr = LB_`t'_mtr
						qui matrix UB_`t'_mtr = UB_`t'_mtr
					}
					else if "`asmptn'" == "mtsr" {
						qui matrix LB_`t'_mtsr = LB_`t'_mtsr
						qui matrix UB_`t'_mtsr = UB_`t'_mtsr
					}
				}
			}

	restore 
	
	}
	
	/*Second, the bounds which do involve the MIV*/
	
	if "`miv'" != "" {
		/*Again preserve and set the estimation sample*/
		
	    preserve
		qui keep if `touse'
		qui keep if inrange(`y',(`ymin'-0.0000001),(`ymax'+0.0000001))
		
		/*Create the bins for the MIV and estimate bin probabilities*/
		if "`discretemiv'" != "" {
			tempvar mivbin
			gen `mivbin' = `miv'
		}
		else {
			tempvar mivbin
			xtile `mivbin' = `miv', nq(`bins')
		}
		
		tempvar binj
		gen `binj' = 0
		forval j=1/`bins' {
			qui replace `binj' = (`mivbin'==`j')
			qui count if `binj'==1
			global nmiv`j' = r(N)
			qui reg `binj' [`weight' `exp']
			qui mat `eb'=e(b)
			global p_`j' = `eb'[1,1]
		}
		
		restore
		
		/*Stack bin probabilities into a vector for later*/
		
		matrix Pj=[${p_1}]
		forval j=2/`bins' {
			matrix Pj = [Pj \ ${p_`j'}]
		}
		
		/*Estimate the bin-specific treatment probabilities and conditional means
		from above*/
		
		forval j=1/`bins' {
			preserve
			qui keep if `touse'
			qui keep if inrange(`y',(`ymin'-0.0000001),(`ymax'+0.0000001))
			
			tempvar i_treat itreat
			qui ta `T', gen(`i_treat')
	
			local i=1
			foreach t of local Tlevels {
				gen `itreat'`t' = `i_treat'`i'
				local ++i
			}
			
			if "`discretemiv'" != "" {
				tempvar mivbin
				gen `mivbin' = `miv'
			}
			else {
				tempvar mivbin
				xtile `mivbin' = `miv', nq(`bins')
			}

			qui keep if `mivbin'==`j'
			
			foreach t of local Tlevels {
				qui reg `itreat'`t' [`weight' `exp']
				qui mat `eb'=e(b)
				global pr_`t'_`j' = `eb'[1,1]
			}
			
			foreach t of local Tlevels {
				tempvar Tgt`t'_`j' Tlt`t'_`j' Tgeq`t'_`j' Tleq`t'_`j'
				qui gen `Tgt`t'_`j'' = (`T'>`t')
				qui gen `Tlt`t'_`j'' = (`T'<`t')
				qui gen `Tgeq`t'_`j'' = (`T'>=`t')
				qui gen `Tleq`t'_`j'' = (`T'<=`t')
			}
	
			foreach t of local Tlevels  {
				qui reg `Tgt`t'_`j'' [`weight' `exp']
				qui mat `eb'=e(b)
				global prTgt`t'_`j' = `eb'[1,1]
		
				qui reg `Tlt`t'_`j'' [`weight' `exp']
				qui mat `eb'=e(b)
				global prTlt`t'_`j' = `eb'[1,1]
		
				qui reg `Tgeq`t'_`j'' [`weight' `exp']
				qui mat `eb'=e(b)
				global prTgeq`t'_`j' = `eb'[1,1]
		
				qui reg `Tleq`t'_`j'' [`weight' `exp']
				qui mat `eb'=e(b)
				global prTleq`t'_`j' = `eb'[1,1]
			}
			

				foreach t of local Tlevels  {
					cap qui reg `y' if `T'==`t' [`weight' `exp']
						if _rc==2000 {
							global m_y`t'_`j' = 0
						}
						else {
							qui mat `eb'=e(b)
							global m_y`t'_`j' = `eb'[1,1]
						}

					cap qui reg `y' if `T' < `t' [`weight' `exp']
						if _rc==2000 {
							global m_yL`t'_`j' = 0
						}
						else {
							qui mat `eb'=e(b)
							global m_yL`t'_`j' = `eb'[1,1]
						}

					cap qui reg `y' if `T' > `t' [`weight' `exp']
						if _rc==2000 {
							global m_yG`t'_`j' = 0
						}
						else {
							qui mat `eb'=e(b)
							global m_yG`t'_`j' = `eb'[1,1]
						}

					cap qui reg `y' if `T' <= `t' [`weight' `exp']
						if _rc==2000 {
							global m_yLeq`t'_`j' = 0
						}
						else {
							qui mat `eb'=e(b)
							global m_yLeq`t'_`j' = `eb'[1,1]
						}

					cap qui reg `y' if `T' >= `t' [`weight' `exp']
						if _rc==2000 {
							global m_yGeq`t'_`j' = 0
						}
						else {
							qui mat `eb'=e(b)
							global m_yGeq`t'_`j' = `eb'[1,1]
						}
				}

				global m_yL`minT'_`j' = 0
				global m_yG`maxT'_`j' = 0
				
				foreach t of local Tlevels {
					if ${pr_`t'_`j'} == 0 {
						global m_y`t'_`j' = ${m_y`t'}
						di in red "Warning: MIV bin `j' for `T'==`t' is empty. The bin `j' mean has been replaced by the overall sample mean for `T'==`t'."
					}
				}
		
		/*Estimate the WC, MTS, MTR, and MTS+MTR bounds within each bin*/
				
		local subminT = `minT' + 1
		local submaxT = `maxT' - 1
		
			/*Non-positive MTS + non-negative MTR*/
		if "`att'" != "" {
				local subminT = `minT' + 1
				local submaxT = `maxT' - 1
				/*Bounds on the ATT, no MIV*/

				/*Non-positive MTS + non-negative MTR*/
				if ("`nmts'" != "" & "`nmtr'" == "") {
					forval r=`subminT'/`maxT' {
						local w=`r'-1
						forval s=`minT'/`w' {
						if "`asmptn'" == "wcv" {
							qui scalar LB_`s'_`r'_wcv_`j' = `ymin'
							qui scalar UB_`s'_`r'_wcv_`j' = `ymax'
							qui scalar LB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsv" {
							qui scalar LB_`s'_`r'_mtsv_`j' = `ymin'
							qui scalar UB_`s'_`r'_mtsv_`j' = ${m_y`s'_`j'}
							qui scalar LB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtrv" {
							qui scalar LB_`s'_`r'_mtrv_`j' = `ymin'
							qui scalar UB_`s'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar LB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsrv" {
							qui scalar LB_`s'_`r'_mtsrv_`j' = `ymin'
							qui scalar UB_`s'_`r'_mtsrv_`j' = min(${m_y`r'_`j'},${m_y`s'_`j'})
							qui scalar LB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
						}
						}
					}
				}

				/*Non-negative MTS + non-positive MTR*/

				else if ("`nmts'" == "" & "`nmtr'" != "") {
					forval r=`subminT'/`maxT' {
						local w=`r'-1
						forval s=`minT'/`w' {
						if "`asmptn'" == "wcv" {
							qui scalar LB_`s'_`r'_wcv_`j' = `ymin'
							qui scalar UB_`s'_`r'_wcv_`j' = `ymax'
							qui scalar LB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsv" {
							qui scalar LB_`s'_`r'_mtsv_`j' = ${m_y`s'_`j'}
							qui scalar UB_`s'_`r'_mtsv_`j' = `ymax'
							qui scalar LB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtrv" {
							qui scalar LB_`s'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`s'_`r'_mtrv_`j' = `ymax'
							qui scalar LB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsrv" {
							qui scalar LB_`s'_`r'_mtsrv_`j' = max(${m_y`r'_`j'},${m_y`s'_`j'})
							qui scalar UB_`s'_`r'_mtsrv_`j' = `ymax'
							qui scalar LB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
						}
						}
					}
				}

				/*Non-positive MTS + non-positive MTR*/

				else if ("`nmts'" != "" & "`nmtr'" != "") {
					forval r=`subminT'/`maxT' {
						local w=`r'-1	
						forval s=`minT'/`w' {
						if "`asmptn'" == "wcv" {
							qui scalar LB_`s'_`r'_wcv_`j' = `ymin'
							qui scalar UB_`s'_`r'_wcv_`j' = `ymax'
							qui scalar LB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsv" {
							qui scalar LB_`s'_`r'_mtsv_`j' = `ymin'
							qui scalar UB_`s'_`r'_mtsv_`j' = ${m_y`s'_`j'}
							qui scalar LB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtrv" {
							qui scalar LB_`s'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`s'_`r'_mtrv_`j' = `ymax'
							qui scalar LB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsrv" {
							qui scalar LB_`s'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`s'_`r'_mtsrv_`j' = ${m_y`s'_`j'}
							qui scalar LB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
						}
						}
					}
				}

				/*Non-negative MTS + non-negative MTR*/

				else {
					forval r=`subminT'/`maxT' {
						local w=`r'-1
						forval s=`minT'/`w' {
						if "`asmptn'" == "wcv" {
							qui scalar LB_`s'_`r'_wcv_`j' = `ymin'
							qui scalar UB_`s'_`r'_wcv_`j' = `ymax'
							qui scalar LB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_wcv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsv" {
							qui scalar LB_`s'_`r'_mtsv_`j' = ${m_y`s'_`j'}
							qui scalar UB_`s'_`r'_mtsv_`j' = `ymax'
							qui scalar LB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtrv" {
							qui scalar LB_`s'_`r'_mtrv_`j' = `ymin'
							qui scalar UB_`s'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar LB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtrv_`j' = ${m_y`r'_`j'}
						}
						else if "`asmptn'" == "mtsrv" {
							qui scalar LB_`s'_`r'_mtsrv_`j' = ${m_y`s'_`j'}
							qui scalar UB_`s'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar LB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
							qui scalar UB_`r'_`r'_mtsrv_`j' = ${m_y`r'_`j'}
						}
						}
					}
				}			
			
		}	
		
		else {
			if ("`nmts'" != "" & "`nmtr'" == "") {
				foreach t of local Tlevels  {	
				    if "`asmptn'" == "wcv" {
						qui scalar LB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymin'*(1-${pr_`t'_`j'})
						qui scalar UB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymax'*(1-${pr_`t'_`j'})
					}
					else if "`asmptn'" == "mtsv" {
						qui scalar LB_`t'_mtsv_`j' = `ymin'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
						qui scalar UB_`t'_mtsv_`j' = `ymax'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
					}
					else if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = ${m_yLeq`t'_`j'}*${prTleq`t'_`j'} + `ymin'*${prTgt`t'_`j'}
						qui scalar UB_`t'_mtrv_`j' = ${m_yGeq`t'_`j'}*${prTgeq`t'_`j'} + `ymax'*${prTlt`t'_`j'}
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = `ymin'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
						qui scalar UB_`t'_mtsrv_`j' = `ymax'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
					}
				}

				forval t=`subminT'/`maxT' {
					local k=`t'-1
					if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = max(LB_`t'_mtrv_`j',UB_`k'_mtrv_`j')
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = max(LB_`t'_mtsrv_`j',UB_`k'_mtsrv_`j')
					}
				}
			}

			/*Non-negative MTS + non-negative MTR*/

			else if ("`nmts'" == "" & "`nmtr'" != "") {
				foreach t of local Tlevels  {
				    if "`asmptn'" == "wcv" {
						qui scalar LB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymin'*(1-${pr_`t'_`j'})
						qui scalar UB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymax'*(1-${pr_`t'_`j'})
					}
					else if "`asmptn'" == "mtsv" {
						qui scalar LB_`t'_mtsv_`j' = `ymin'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
						qui scalar UB_`t'_mtsv_`j' = `ymax'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
					}
					else if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = ${m_yGeq`t'_`j'}*${prTgeq`t'_`j'} + `ymin'*${prTlt`t'_`j'}
						qui scalar UB_`t'_mtrv_`j' = ${m_yLeq`t'_`j'}*${prTleq`t'_`j'} + `ymax'*${prTgt`t'_`j'}
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = `ymin'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
						qui scalar UB_`t'_mtsrv_`j' = `ymax'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
					}
				}

				forval t=`subminT'/`maxT' {
					local k=`t'-1
					if "`asmptn'" == "mtrv" {
						qui scalar UB_`t'_mtrv_`j' = min(UB_`t'_mtrv_`j',LB_`k'_mtrv_`j')
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar UB_`t'_mtsrv_`j' = min(UB_`t'_mtsrv_`j',LB_`k'_mtsrv_`j')
					}
				}
			}

			/*Non-positive MTS + non-negative MTR*/

			else if ("`nmts'" != "" & "`nmtr'" != "") {
				foreach t of local Tlevels  {	
				    if "`asmptn'" == "wcv" {
						qui scalar LB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymin'*(1-${pr_`t'_`j'})
						qui scalar UB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymax'*(1-${pr_`t'_`j'})
					}
					else if "`asmptn'" == "mtsv" {
						qui scalar LB_`t'_mtsv_`j' = `ymin'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
						qui scalar UB_`t'_mtsv_`j' = `ymax'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
					}
					else if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = ${m_yGeq`t'_`j'}*${prTgeq`t'_`j'} + `ymin'*${prTlt`t'_`j'}
						qui scalar UB_`t'_mtrv_`j' = ${m_yLeq`t'_`j'}*${prTleq`t'_`j'} + `ymax'*${prTgt`t'_`j'}
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = ${m_yG`t'_`j'}*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
						qui scalar UB_`t'_mtsrv_`j' = ${m_yL`t'_`j'}*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
					}
				}

				forval t=`subminT'/`maxT' {
					local k=`t'-1
					if "`asmptn'" == "mtrv" {
						qui scalar UB_`t'_mtrv_`j' = min(UB_`t'_mtrv_`j',LB_`k'_mtrv_`j')
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar UB_`t'_mtsrv_`j' = min(UB_`t'_mtsrv_`j',LB_`k'_mtsrv_`j')
					}
				}
			}

			/*Non-negative MTS + non-negative MTR*/

			else {
				foreach t of local Tlevels  {	
				    if "`asmptn'" == "wcv" {
						qui scalar LB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymin'*(1-${pr_`t'_`j'})
						qui scalar UB_`t'_wcv_`j' = ${m_y`t'_`j'}*${pr_`t'_`j'} + `ymax'*(1-${pr_`t'_`j'})
					}
					else if "`asmptn'" == "mtsv" {
						qui scalar LB_`t'_mtsv_`j' = `ymin'*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
						qui scalar UB_`t'_mtsv_`j' = `ymax'*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
					}
					else if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = ${m_yLeq`t'_`j'}*${prTleq`t'_`j'} + `ymin'*${prTgt`t'_`j'}
						qui scalar UB_`t'_mtrv_`j' = ${m_yGeq`t'_`j'}*${prTgeq`t'_`j'} + `ymax'*${prTlt`t'_`j'}
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = ${m_yL`t'_`j'}*${prTlt`t'_`j'} + ${m_y`t'_`j'}*${prTgeq`t'_`j'}
						qui scalar UB_`t'_mtsrv_`j' = ${m_yG`t'_`j'}*${prTgt`t'_`j'} + ${m_y`t'_`j'}*${prTleq`t'_`j'}
					}
				}

				forval t=`subminT'/`maxT' {
					local k=`t'-1
					if "`asmptn'" == "mtrv" {
						qui scalar LB_`t'_mtrv_`j' = max(LB_`t'_mtrv_`j',UB_`k'_mtrv_`j')
					}
					else if "`asmptn'" == "mtsrv" {
						qui scalar LB_`t'_mtsrv_`j' = max(LB_`t'_mtsrv_`j',UB_`k'_mtsrv_`j')
					}
				}
			}
		}
	
		restore
		}
	}
	
	/*This section produces the CLR bounding functions by treating each bound
	as its own 'dataset', doing some tricks with sorting to produce the string
	of possible evaluations of the MIV estimator max/min operators and then
	store into a matrix. More information on what is being done is available
	from Giuseppe*/
	
	if "`miv'" != "" {
	    /*non-positive MIV*/
	if "`nmiv'" != "" {
		if "`att'" != "" {
		if "`ul'" == "U" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
			
				qui gen id=_n
				qui gen pos0 = UB_`contrs'_`contrt'_`asmptn'_1
				global sorterU pos0
						
				forval j=1/`bins' {
					local i=`j'-1
					qui bysort ${sorterU}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = UB_`contrs'_`contrt'_`asmptn'_`j' if half`j'==2
					
					qui drop half`j'
					global sorterU ${sorterU} pos`j'
				}
						
				mkmat pos1-pos`bins', mat(U_`contrs'_`contrt'_`asmptn')
			restore
		}
		
		/*Prepare lower bound bounding functions*/
		
		else if "`ul'" == "L" {			
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
						
				qui gen id=_n
				local init = `bins' + 1
				qui gen pos`init' = -LB_`contrs'_`contrt'_`asmptn'_`bins'
				global sorterL pos`init'
						
				forval j=`bins' (-1) 1 {
					local i=`j'+1
					qui bysort ${sorterL}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = -LB_`contrs'_`contrt'_`asmptn'_`j' if half`j'==2
							
					qui drop half`j'
					global sorterL ${sorterL} pos`j'
				}
						
				mkmat pos`bins'-pos1, mat(L_`contrs'_`contrt'_`asmptn')
			restore
				
			matselrc L_`contrs'_`contrt'_`asmptn' L_`contrs'_`contrt'_`asmptn', c(`bins'/1)
			matrix L_`contrs'_`contrt'_`asmptn' = -L_`contrs'_`contrt'_`asmptn'
		}	
		}
		
		else {
	    /*Prepare upper bound bounding functions*/
		if "`ul'" == "U" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
			
				qui gen id=_n
				qui gen pos0 = UB_`contrt'_`asmptn'_1
				global sorterU pos0
						
				forval j=1/`bins' {
					local i=`j'-1
					qui bysort ${sorterU}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = UB_`contrt'_`asmptn'_`j' if half`j'==2
					
					qui drop half`j'
					global sorterU ${sorterU} pos`j'
				}
						
				mkmat pos1-pos`bins', mat(U_`contrt'_`asmptn')
			restore
		}
		
		/*Prepare lower bound bounding functions*/
		
		else if "`ul'" == "L" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
						
				qui gen id=_n
				local init = `bins' + 1
				qui gen pos`init' = -LB_`contrt'_`asmptn'_`bins'
				global sorterL pos`init'
						
				forval j=`bins' (-1) 1 {
					local i=`j'+1
					qui bysort ${sorterL}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = -LB_`contrt'_`asmptn'_`j' if half`j'==2
							
					qui drop half`j'
					global sorterL ${sorterL} pos`j'
				}
						
				mkmat pos`bins'-pos1, mat(L_`contrt'_`asmptn')
			restore
				
			matselrc L_`contrt'_`asmptn' L_`contrt'_`asmptn', c(`bins'/1)
			matrix L_`contrt'_`asmptn' = -L_`contrt'_`asmptn'
		}
		
	}
	}
	
	/*Non-negative MIV*/
	
	else {

		/*Prepare lower bound bounding functions*/
		
		if "`att'" != "" {
		if "`ul'" == "L" {
preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
			
				qui gen id=_n
				qui gen pos0 = LB_`contrs'_`contrt'_`asmptn'_1
				global sorterL pos0
						
				forval j=1/`bins' {
					local i=`j'-1
					qui bysort ${sorterL}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = LB_`contrs'_`contrt'_`asmptn'_`j' if half`j'==2
					
					qui drop half`j'
					global sorterL ${sorterL} pos`j'
				}
						
				mkmat pos1-pos`bins', mat(L_`contrs'_`contrt'_`asmptn')
			restore
		}
		
		/*Prepare upper bound bounding functions*/
		
		else if "`ul'" == "U" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
						
				qui gen id=_n
				local init = `bins' + 1
				qui gen pos`init' = -UB_`contrs'_`contrt'_`asmptn'_`bins'
				global sorterU pos`init'
						
				forval j=`bins' (-1) 1 {
					local i=`j'+1
					qui bysort ${sorterU}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = -UB_`contrs'_`contrt'_`asmptn'_`j' if half`j'==2
							
					qui drop half`j'
					global sorterU ${sorterU} pos`j'
				}
						
				mkmat pos`bins'-pos1, mat(U_`contrs'_`contrt'_`asmptn')
			restore
				
			matselrc U_`contrs'_`contrt'_`asmptn' U_`contrs'_`contrt'_`asmptn', c(`bins'/1)
			matrix U_`contrs'_`contrt'_`asmptn' = -U_`contrs'_`contrt'_`asmptn'
		}			
		}
	
	
		else { /*ATE form*/
		if "`ul'" == "L" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
			
				qui gen id=_n
				qui gen pos0 = LB_`contrt'_`asmptn'_1
				global sorterL pos0
						
				forval j=1/`bins' {
					local i=`j'-1
					qui bysort ${sorterL}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = LB_`contrt'_`asmptn'_`j' if half`j'==2
					
					qui drop half`j'
					global sorterL ${sorterL} pos`j'
				}
						
				mkmat pos1-pos`bins', mat(L_`contrt'_`asmptn')
			restore
		}
		
		/*Prepare upper bound bounding functions*/
		
		else if "`ul'" == "U" {
			preserve
				clear
				local tterms = 2^(`bins'-1)
				qui set obs `tterms'
						
				qui gen id=_n
				local init = `bins' + 1
				qui gen pos`init' = -UB_`contrt'_`asmptn'_`bins'
				global sorterU pos`init'
						
				forval j=`bins' (-1) 1 {
					local i=`j'+1
					qui bysort ${sorterU}: replace id=_n
					qui egen half`j' = xtile(id), by(pos`i') nq(2)
					qui gen pos`j' = pos`i' if half`j'==1
					qui replace pos`j' = -UB_`contrt'_`asmptn'_`j' if half`j'==2
							
					qui drop half`j'
					global sorterU ${sorterU} pos`j'
				}
						
				mkmat pos`bins'-pos1, mat(U_`contrt'_`asmptn')
			restore
				
			matselrc U_`contrt'_`asmptn' U_`contrt'_`asmptn', c(`bins'/1)
			matrix U_`contrt'_`asmptn' = -U_`contrt'_`asmptn'
		}
		}
	}
	
	/*Multiply by bin-probabilities vector in next step of building bounding functions*/
	
	if "`att'" != "" {
		if "`ul'" == "L" {
			matrix LL_`contrs'_`contrt'_`asmptn' = L_`contrs'_`contrt'_`asmptn'*Pj
		}
		if "`ul'" == "U" {
			matrix UU_`contrs'_`contrt'_`asmptn' = U_`contrs'_`contrt'_`asmptn'*Pj
		}		
	}
	
	else {
		if "`ul'" == "L" {
			matrix LL_`contrt'_`asmptn' = L_`contrt'_`asmptn'*Pj
		}
		if "`ul'" == "U" {
			matrix UU_`contrt'_`asmptn' = U_`contrt'_`asmptn'*Pj
		}
	}

	/*Need to use Mata for the larger matrix size when we do the subtraction step
	making these into ATE bounding functions*/
			
	scalar tterms = 2^(`bins'-1)
	mata: tterms = st_numscalar("tterms")
	mata: sqterms = tterms^2
		
	if "`att'" != "" {
		if "`ul'" == "L" {
			mata: LL_`contrs'_`contrt'_`asmptn' = st_matrix("LL_`contrs'_`contrt'_`asmptn'")
			
			mata: LB_`contrs'_`contrt'_`asmptn' = (LL_`contrs'_`contrt'_`asmptn')'
		
			mata: LB_`contrs'_`contrt'_`asmptn' = (LB_`contrs'_`contrt'_`asmptn', max(LB_`contrs'_`contrt'_`asmptn'))
		
			mata: st_matrix("LB_`contrs'_`contrt'_`asmptn'",LB_`contrs'_`contrt'_`asmptn')
		}
	
		if "`ul'" == "U" {
			mata: UU_`contrs'_`contrt'_`asmptn' = st_matrix("UU_`contrs'_`contrt'_`asmptn'")
			
			mata: UB_`contrs'_`contrt'_`asmptn' = (UU_`contrs'_`contrt'_`asmptn')'
		
			mata: UB_`contrs'_`contrt'_`asmptn' = (UB_`contrs'_`contrt'_`asmptn', min(UB_`contrs'_`contrt'_`asmptn'))
			
			mata: st_matrix("UB_`contrs'_`contrt'_`asmptn'",UB_`contrs'_`contrt'_`asmptn')
		}
	}

	else {
		if "`ul'" == "L" {
			mata: LL_`contrt'_`asmptn' = st_matrix("LL_`contrt'_`asmptn'")
			
			mata: LB_`contrt'_`asmptn' = (LL_`contrt'_`asmptn')'
		
			mata: LB_`contrt'_`asmptn' = (LB_`contrt'_`asmptn', max(LB_`contrt'_`asmptn'))
		
			mata: st_matrix("LB_`contrt'_`asmptn'",LB_`contrt'_`asmptn')
		}
	
		if "`ul'" == "U" {
			mata: UU_`contrt'_`asmptn' = st_matrix("UU_`contrt'_`asmptn'")
			
			mata: UB_`contrt'_`asmptn' = (UU_`contrt'_`asmptn')'
		
			mata: UB_`contrt'_`asmptn' = (UB_`contrt'_`asmptn', min(UB_`contrt'_`asmptn'))
		
			mata: st_matrix("UB_`contrt'_`asmptn'",UB_`contrt'_`asmptn')
		}
	}
	}

	/*Return the vectors of bounding functions for the main 'mpclr' command to
	pass into the CLR program*/
	
	if "`att'" != "" {
		if "`ul'" == "L" {	
			ereturn post LB_`contrs'_`contrt'_`asmptn', esample(`touse') properties("b")
		}
	
		if "`ul'" == "U" {
			ereturn post UB_`contrs'_`contrt'_`asmptn', esample(`touse') properties("b")
		}
	}
	
	else {
		if "`ul'" == "L" {	
			ereturn post LB_`contrt'_`asmptn', esample(`touse') properties("b")
		}

		if "`ul'" == "U" {
			ereturn post UB_`contrt'_`asmptn', esample(`touse') properties("b")
		}
	}



	
end
	
	
