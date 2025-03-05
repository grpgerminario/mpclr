program pompclr, eclass properties(svyb)
	version 10.0
	#delimit ;
	syntax varlist (min=2 max=2 fv) [if] [in] [pw iw],  YMIN(string) YMAX(string) [ATT BOUNDS(string) NMTS NMTR MIV(varname) MIVBOUNDS(string) bins(integer 5) DISCRETEMIV NMIV reps(integer 25) seed(integer 22176) level(integer 95) NOISYCLR LATEX(string) NPSU(int -1) SURVEY BRR UNCORRECTED];
	marksample touse;
	#delimit cr
	
	/*Gather variables/options for use in the program*/

	gettoken y T : varlist
	
	tokenize "`ymin'"
	tokenize "`ymax'"
	tokenize "`bounds'"
	tokenize "`mivbounds'"
	tokenize "`bins'"
	tokenize "`reps'"
	tokenize "`seed'"
	tokenize "`level'"
	tokenize "`moret'"
	tokenize "`latex'"
	tokenize "`npsu'"
	
	/*Preliminary checks*/
	cap _svy_newrule
	if "`survey'" != "" & _rc==119 {
	    di in red "Please svyset the data to use the survey options."
		exit 119
	}
	
	cap _svy_newrule
	if "`brr'" != "" & _rc==119 {
		di in red "Please svyset the data with BRR weights to use the brr option."
	}
	
	/*Gather survey info*/
	
	if "`survey'" != "" {
	    qui svyset
		global svysettings = r(settings)
		
		qui svydes
		local nstrata = r(N_strata)
		if `nstrata' > 1 {
		    qui bsweights bsw, reps(`reps') n(`npsu') seed(`seed') replace
		}
		else {
		    qui bsweights bsw, reps(`reps') n(`npsu') seed(`seed') replace nosvy
		}
		
		qui svyset ${svysettings} bsrweight(bsw*)
	}
/*Gather the min/max treatment levels; subminT is one above min and helpful in loops*/
	
	qui sum `T' if `touse'
	local minT = r(min)
	local maxT = r(max)
	local subminT = `minT' + 1
	local submaxT = `maxT' - 1
	
	/*Store all levels of the treatment*/
	
	qui levelsof `T' if `touse', local(Tlevels)
		
	/*Begin to compute the non-MIV bounds, separate by lower/upper and by assumption used;
	(mtsr = mts+mtr)*/
	
	
	if "`bounds'" != "" {
		foreach m of local bounds {
		    if "`m'" == "wc" {
				local bset1 wc
			}
			else if "`m'" == "mts" {
				local bset2 mts
			}
			else if "`m'" == "mtr" {
				local bset3 mtr
			}
			else if "`m'" == "mtsr" {
				local bset4 mtsr
			}
			else if "`m'" == "none" {
			    local bset5 none
			}
		}
		
		global bset "`bset1' `bset2' `bset3' `bset4' `bset5'"
	}
	else {
		global bset "wc mts mtr mtsr"
	}
	
	if "${bset}" == "    none" {
		di ""
		di ""
	    di "No non-MIV bounds requested"
	}
	else {
	
		local asmpt ${bset}
		local lowup "L U"
		
		if "`att'" != "" {

		if "`survey'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						svy bootstrap _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')
							
						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)

						matrix se`x'_`k'_`t'_`a' = e(se)
							
					}
				}	
			}
		}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						svy brr _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						
						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)

					}
				}		
			}
		}
		}


		else {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)

						matrix se`x'_`k'_`t'_`a' = e(se)
					}
				}		
			}
		}
	}
	} /*end of ATT if*/
	
	else {
		if "`survey'" != "" {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						svy bootstrap _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')
							
						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)

						matrix se`x'_`t'_`a' = e(se)
							
					}
				}
			}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						svy brr _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

							
						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)
						
					}	
				}
			}
		}


		else {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'

						/*This bootstraps the 'pompbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)

						matrix se`x'_`t'_`a' = e(se)
					}	
				}
			}
		}
	}
	}
	
	
	/*Next are the bounds which use the MIV, again by lower/upper and by assumption
	(the 'v' at the end of each denotes 'with MIV')*/
	
	if "`mivbounds'" != "" {
		foreach m of local mivbounds {
		    if "`m'" == "wcv" {
				local mbset1 wcv
			}
			else if "`m'" == "mtsv" {
				local mbset2 mtsv
			}
			else if "`m'" == "mtrv" {
				local mbset3 mtrv
			}
			else if "`m'" == "mtsrv" {
				local mbset4 mtsrv
			}
		}
		
		global mbset "`mbset1' `mbset2' `mbset3' `mbset4'"
	}
	else {
		global mbset "wcv mtsv mtrv mtsrv"
	}
	
	if "`miv'" == "" {
		di ""
		di ""
		di "No MIV bounds requested"
	}
	
	else {
	
		local asmpt ${mbset}
		local lowup "L U"

		if "`att'" != "" {
		if "`survey'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'

						di ""
						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						svy bootstrap _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)

						matrix se`x'_`k'_`t'_`a' = e(se)
					}
				}		
			}
		}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'

						di ""
						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						svy brr _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)
					}
				}		
			}
		}
		}

		else {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'-1

						di ""
						di ""
						di "Begin bootstrap for `x'B_`k'_`t'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`k'_`t'_`a' = e(b)'

						matrix Vb`x'_`k'_`t'_`a' = e(V)

						matrix se`x'_`k'_`t'_`a' = e(se)
					}
				}		
			}
		}
	}
	} /*end of ATT if*/
	
	else {
		if "`survey'" != "" {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'

						di ""
						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						svy bootstrap _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)

						matrix se`x'_`t'_`a' = e(se)
					}
				}		
			}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'

						di ""
						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						svy brr _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)
					}
				}		
			}
		}

		else {
			foreach a of local asmpt {
				foreach t of local Tlevels {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'

						di ""
						di ""
						di "Begin bootstrap for `x'B_`t'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  pompbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `discretemiv' `nmiv' ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'_`t'_`a' = e(b)'

						matrix Vb`x'_`t'_`a' = e(V)

						matrix se`x'_`t'_`a' = e(se)
					}
				}		
			}
		}
	}
	}
	
	/*Get estimation sample size for use in tables and in CLR*/
	
	qui count if `touse' & inrange(`y',(`ymin'-0.0000001),(`ymax'+0.0000001))
	qui scalar Nob = r(N)
	qui ereturn scalar Nob = Nob
	
	if "`survey'" != "" {
		qui svyset ${svysettings}
		drop bsw*
	}

	
	
********************************************************************************
********************************************************************************
/*Now move on to running the CLR procedure for bias correction and the CIs*/

	/*Indicate to user that CLR has begun to run*/
	
	di ""
	di "Beginning CLR Procedure..."
	di "{hline 30}"
	di "{hline 30}"
	
	/*set (1-\alpha)% CIs from user input; default being level=95*/
	
	qui scalar alpha=1-0.`level'
	
	if "${bset}" == "    none" & "`miv'" != "" {
		local asmpt "${mbset}"
	}
	else if "${bset}" != "    none" & "`miv'" == "" {
	    local asmpt "${bset}"
	}
	else {
		local asmpt "${bset} ${mbset}"
	}
	
	/*Display some heading for which bounds are being computed when 'noisyclr' is used*/
	
	if "${bset}" == "    none" & "`miv'" == "" {
		di ""
	}
	else {
		if "`noisyclr'" != "" {
			
			if "`att'" != "" {
				scalar tterms = 2^(`bins'-1)
				scalar sqterms = tterms^2
				scalar uncterms = sqterms+1
				local sqterms = sqterms
				local uncterms = uncterms
			
				foreach m of global mbset {
					forval t=`subminT'/`maxT' { /*note: for ATT forms, t is the condition*/
						forval k=`minT'/`t' {
						matselrc bL_`k'_`t'_`m' unc_bL_`k'_`t'_`m', r(`uncterms')
						matselrc bL_`k'_`t'_`m' bL_`k'_`t'_`m', r(1/`sqterms')
					
						matselrc VbL_`k'_`t'_`m' unc_VbL_`k'_`t'_`m', r(`uncterms') c(`uncterms')
						matselrc VbL_`k'_`t'_`m' VbL_`k'_`t'_`m', r(1/`sqterms') c(1/`sqterms')
					
						matselrc bU_`k'_`t'_`m' unc_bU_`k'_`t'_`m', r(`uncterms')
						matselrc bU_`k'_`t'_`m' bU_`k'_`t'_`m', r(1/`sqterms')
					
						matselrc VbU_`k'_`t'_`m' unc_VbU_`k'_`t'_`m', r(`uncterms') c(`uncterms')
						matselrc VbU_`k'_`t'_`m' VbU_`k'_`t'_`m', r(1/`sqterms') c(1/`sqterms')
						}
					}
				}

				foreach a of local asmpt {
					forval t=`subminT'/`maxT' {
						forval k=`minT'/`t' {
						matrix bL = bL_`k'_`t'_`a' 
						matrix VbL = VbL_`k'_`t'_`a'
						matrix bU = bU_`k'_`t'_`a' 
						matrix VbU = VbU_`k'_`t'_`a' 
						

					/*define the 'current assumption' being used*/

						qui scalar crntasmpt = "`a'"

						if crntasmpt=="wc" {
							di ""
							di "Begin CLR for Worst-Case for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "Begin CLR for MTS for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "Begin CLR for MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "Begin CLR for MTS+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "Begin CLR for MIV-Only for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "Begin CLR for MIV+MTS for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "Begin CLR for MIV+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "Begin CLR for MIV+MTS+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}


					/*Run the CLR command from Flores & Wang;
					   100000 is the number of 'simulations' R in the CLR algorithm*/


						CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar LB_`k'_`t'_`a' = eLstar[1,1]	
						qui scalar UB_`k'_`t'_`a' = eUstar[1,1]
						qui scalar ciLB_`k'_`t'_`a' = eCIL[1,1]
						qui scalar ciUB_`k'_`t'_`a' = eCIU[1,1]

						if crntasmpt=="wc" {
							di ""
							di "End CLR for Worst-Case for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "End CLR for MTS for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "End CLR for MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "End CLR for MTS+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "End CLR for MIV-Only for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "End CLR for MIV+MTS for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "End CLR for MIV+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "End CLR for MIV+MTS+MTR for E[Y(`k') | T=`t']"
							di "{hline 30}"
						}
						}
					}
				}
			} /*end of att if*/
			
			else {
			
				scalar tterms = 2^(`bins'-1)
				scalar sqterms = tterms^2
				scalar uncterms = sqterms+1
				local sqterms = sqterms
				local uncterms = uncterms
			
				foreach m of global mbset {
					foreach t of local Tlevels {
						local k=`t'-1
						matselrc bL_`t'_`m' unc_bL_`t'_`m', r(`uncterms')
						matselrc bL_`t'_`m' bL_`t'_`m', r(1/`sqterms')
					
						matselrc VbL_`t'_`m' unc_VbL_`t'_`m', r(`uncterms') c(`uncterms')
						matselrc VbL_`t'_`m' VbL_`t'_`m', r(1/`sqterms') c(1/`sqterms')
					
						matselrc bU_`t'_`m' unc_bU_`t'_`m', r(`uncterms')
						matselrc bU_`t'_`m' bU_`t'_`m', r(1/`sqterms')
					
						matselrc VbU_`t'_`m' unc_VbU_`t'_`m', r(`uncterms') c(`uncterms')
						matselrc VbU_`t'_`m' VbU_`t'_`m', r(1/`sqterms') c(1/`sqterms')
					}
				}

				foreach a of local asmpt {
					foreach t of local Tlevels {
						local k=`t'-1
						matrix bL = bL_`t'_`a' 
						matrix VbL = VbL_`t'_`a'
						matrix bU = bU_`t'_`a' 
						matrix VbU = VbU_`t'_`a' 

					/*define the 'current assumption' being used*/

						qui scalar crntasmpt = "`a'"

						if crntasmpt=="wc" {
							di ""
							di "Begin CLR for Worst-Case for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "Begin CLR for MTS for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "Begin CLR for MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "Begin CLR for MTS+MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "Begin CLR for MIV-Only for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "Begin CLR for MIV+MTS for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "Begin CLR for MIV+MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "Begin CLR for MIV+MTS+MTR for E[Y(`t')]"
							di "{hline 30}"
						}


					/*Run the CLR command from Flores & Wang;
					   100000 is the number of 'simulations' R in the CLR algorithm*/


						CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar LB_`t'_`a' = eLstar[1,1]	
						qui scalar UB_`t'_`a' = eUstar[1,1]
						qui scalar ciLB_`t'_`a' = eCIL[1,1]
						qui scalar ciUB_`t'_`a' = eCIU[1,1]

						if crntasmpt=="wc" {
							di ""
							di "End CLR for Worst-Case for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "End CLR for MTS for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "End CLR for MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "End CLR for MTS+MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "End CLR for MIV-Only for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "End CLR for MIV+MTS for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "End CLR for MIV+MTR for E[Y(`t')]"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "End CLR for MIV+MTS+MTR for E[Y(`t')]"
							di "{hline 30}"
						}

					}
				}

			} /*end ate if*/
		} /*end noisyclr if*/
	

		/*If 'noisyclr' is NOT requested, then everything is done quietly...
			   this makes the Viewer a lot more readable; if noisy is requested, we
		   suggest starting a log file before running the 'mpclr' command*/

		else {	
			if "`att'" != "" {
			scalar tterms = 2^(`bins'-1)
			scalar sqterms = tterms^2
			*scalar uncterms = sqterms+1
			*local sqterms = sqterms
			*local uncterms = uncterms
			
			foreach m of global mbset {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					
					scalar unctermsL = rowsof(bL_`k'_`t'_`m')
					scalar ftermsL = unctermsL-1
					local unctermsL = unctermsL
					local ftermsL = ftermsL
					
					matselrc bL_`k'_`t'_`m' unc_bL_`k'_`t'_`m', r(`unctermsL')
					matselrc bL_`k'_`t'_`m' bL_`k'_`t'_`m', r(1/`ftermsL')
					
					matselrc VbL_`k'_`t'_`m' unc_VbL_`k'_`t'_`m', r(`unctermsL') c(`unctermsL')
					matselrc VbL_`k'_`t'_`m' VbL_`k'_`t'_`m', r(1/`ftermsL') c(1/`ftermsL')
					
					scalar unctermsU = rowsof(bU_`k'_`t'_`m')
					scalar ftermsU = unctermsU-1
					local unctermsU = unctermsU
					local ftermsU = ftermsU
					
					matselrc bU_`k'_`t'_`m' unc_bU_`k'_`t'_`m', r(`unctermsU')
					matselrc bU_`k'_`t'_`m' bU_`k'_`t'_`m', r(1/`ftermsU')
					
					matselrc VbU_`k'_`t'_`m' unc_VbU_`k'_`t'_`m', r(`unctermsU') c(`unctermsU')
					matselrc VbU_`k'_`t'_`m' VbU_`k'_`t'_`m', r(1/`ftermsU') c(1/`ftermsU')
				}
				}
			}

			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
						
					matrix bL = bL_`k'_`t'_`a' 
					matrix VbL = VbL_`k'_`t'_`a'
					matrix bU = bU_`k'_`t'_`a' 
					matrix VbU = VbU_`k'_`t'_`a' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar LB_`k'_`t'_`a' = eLstar[1,1]	
					qui scalar UB_`k'_`t'_`a' = eUstar[1,1]
					qui scalar ciLB_`k'_`t'_`a' = eCIL[1,1]
					qui scalar ciUB_`k'_`t'_`a' = eCIU[1,1]
				}
				}
			}
			
			} /*end att if*/
			
			else {
			scalar tterms = 2^(`bins'-1)
			scalar sqterms = tterms^2
			*scalar uncterms = sqterms+1
			*local sqterms = sqterms
			*local uncterms = uncterms
			
			foreach m of global mbset {
				foreach t of local Tlevels {
					
					scalar unctermsL = rowsof(bL_`t'_`m')
					scalar ftermsL = unctermsL-1
					local unctermsL = unctermsL
					local ftermsL = ftermsL
					
					matselrc bL_`t'_`m' unc_bL_`t'_`m', r(`unctermsL')
					matselrc bL_`t'_`m' bL_`t'_`m', r(1/`ftermsL')
					
					matselrc VbL_`t'_`m' unc_VbL_`t'_`m', r(`unctermsL') c(`unctermsL')
					matselrc VbL_`t'_`m' VbL_`t'_`m', r(1/`ftermsL') c(1/`ftermsL')
					
					scalar unctermsU = rowsof(bU_`t'_`m')
					scalar ftermsU = unctermsU-1
					local unctermsU = unctermsU
					local ftermsU = ftermsU
					
					matselrc bU_`t'_`m' unc_bU_`t'_`m', r(`unctermsU')
					matselrc bU_`t'_`m' bU_`t'_`m', r(1/`ftermsU')
					
					matselrc VbU_`t'_`m' unc_VbU_`t'_`m', r(`unctermsU') c(`unctermsU')
					matselrc VbU_`t'_`m' VbU_`t'_`m', r(1/`ftermsU') c(1/`ftermsU')
				}
			}

			foreach a of local asmpt {
				foreach t of local Tlevels {
					
					matrix bL = bL_`t'_`a' 
					matrix VbL = VbL_`t'_`a'
					matrix bU = bU_`t'_`a' 
					matrix VbU = VbU_`t'_`a' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar LB_`t'_`a' = eLstar[1,1]	
					qui scalar UB_`t'_`a' = eUstar[1,1]
					qui scalar ciLB_`t'_`a' = eCIL[1,1]
					qui scalar ciUB_`t'_`a' = eCIU[1,1]
				}
			}

		} /*end att's else (aka end ate)*/
		
		if "`uncorrected'" != "" {
			if "`att'" != "" {
				foreach m of global mbset {
				forval t=`subminT'/`maxT' {
					forval k=`minT'/`t' {
					matrix bL = unc_bL_`k'_`t'_`m' 
					matrix VbL = unc_VbL_`k'_`t'_`m'
					matrix bU = unc_bU_`k'_`t'_`m' 
					matrix VbU = unc_VbU_`k'_`t'_`m' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar unc_LB_`k'_`t'_`m' = eLstar[1,1]	
					qui scalar unc_UB_`k'_`t'_`m' = eUstar[1,1]
					qui scalar unc_ciLB_`k'_`t'_`m' = eCIL[1,1]
					qui scalar unc_ciUB_`k'_`t'_`m' = eCIU[1,1]
				}
				}
			}
			}
			else {
			foreach m of global mbset {
				foreach t of local Tlevels {
					local k=`t'-1
					matrix bL = unc_bL_`t'_`m' 
					matrix VbL = unc_VbL_`t'_`m'
					matrix bU = unc_bU_`t'_`m' 
					matrix VbU = unc_VbU_`t'_`m' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar unc_LB_`t'_`m' = eLstar[1,1]	
					qui scalar unc_UB_`t'_`m' = eUstar[1,1]
					qui scalar unc_ciLB_`t'_`m' = eCIL[1,1]
					qui scalar unc_ciUB_`t'_`m' = eCIU[1,1]
				}
			}
			}
		}

	} /*close else for noisyclr*/
	} /*close else for bound sets*/
	

	/*re-store the bounds as estimates*/
	
	if "`att'" != "" {
		forval t=`subminT'/`maxT' {
		forval k=`minT'/`t' {
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar LB_`k'_`t'_wc = LB_`k'_`t'_wc
			qui ereturn scalar UB_`k'_`t'_wc = UB_`k'_`t'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar LB_`k'_`t'_mts = LB_`k'_`t'_mts
			qui ereturn scalar UB_`k'_`t'_mts = UB_`k'_`t'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar LB_`k'_`t'_mts = LB_`k'_`t'_mts
			qui ereturn scalar UB_`k'_`t'_mts = UB_`k'_`t'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar LB_`k'_`t'_mtr = LB_`k'_`t'_mtr
			qui ereturn scalar UB_`k'_`t'_mtr = UB_`k'_`t'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar LB_`k'_`t'_mtsr = LB_`k'_`t'_mtsr
			qui ereturn scalar UB_`k'_`t'_mtsr = UB_`k'_`t'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
			di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {
				qui ereturn scalar LB_`k'_`t'_mivo = LB_`k'_`t'_wcv
				qui ereturn scalar UB_`k'_`t'_mivo = UB_`k'_`t'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`k'_`t'_mivo = unc_LB_`k'_`t'_wcv
					qui ereturn scalar unc_UB_`k'_`t'_mivo = unc_UB_`k'_`t'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar LB_`k'_`t'_mivs = LB_`k'_`t'_mtsv
				qui ereturn scalar UB_`k'_`t'_mivs = UB_`k'_`t'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`k'_`t'_mivs = unc_LB_`k'_`t'_mtsv
					qui ereturn scalar unc_UB_`k'_`t'_mivs = unc_UB_`k'_`t'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar LB_`k'_`t'_mivr = LB_`k'_`t'_mtrv
				qui ereturn scalar UB_`k'_`t'_mivr = UB_`k'_`t'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`k'_`t'_mivr = unc_LB_`k'_`t'_mtrv
					qui ereturn scalar unc_UB_`k'_`t'_mivr = unc_UB_`k'_`t'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar LB_`k'_`t'_mivsr = LB_`k'_`t'_mtsrv
				qui ereturn scalar UB_`k'_`t'_mivsr = UB_`k'_`t'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`k'_`t'_mivsr = unc_LB_`k'_`t'_mtsrv
					qui ereturn scalar unc_UB_`k'_`t'_mivsr = unc_UB_`k'_`t'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}
	}
	}
	
	else {
		
	foreach t of local Tlevels {
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar LB_`t'_wc = LB_`t'_wc
			qui ereturn scalar UB_`t'_wc = UB_`t'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar LB_`t'_mts = LB_`t'_mts
			qui ereturn scalar UB_`t'_mts = UB_`t'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar LB_`t'_mts = LB_`t'_mts
			qui ereturn scalar UB_`t'_mts = UB_`t'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar LB_`t'_mtr = LB_`t'_mtr
			qui ereturn scalar UB_`t'_mtr = UB_`t'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar LB_`t'_mtsr = LB_`t'_mtsr
			qui ereturn scalar UB_`t'_mtsr = UB_`t'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
			di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {
				qui ereturn scalar LB_`t'_mivo = LB_`t'_wcv
				qui ereturn scalar UB_`t'_mivo = UB_`t'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`t'_mivo = unc_LB_`t'_wcv
					qui ereturn scalar unc_UB_`t'_mivo = unc_UB_`t'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar LB_`t'_mivs = LB_`t'_mtsv
				qui ereturn scalar UB_`t'_mivs = UB_`t'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`t'_mivs = unc_LB_`t'_mtsv
					qui ereturn scalar unc_UB_`t'_mivs = unc_UB_`t'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar LB_`t'_mivr = LB_`t'_mtrv
				qui ereturn scalar UB_`t'_mivr = UB_`t'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`t'_mivr = unc_LB_`t'_mtrv
					qui ereturn scalar unc_UB_`t'_mivr = unc_UB_`t'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar LB_`t'_mivsr = LB_`t'_mtsrv
				qui ereturn scalar UB_`t'_mivsr = UB_`t'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB_`t'_mivsr = unc_LB_`t'_mtsrv
					qui ereturn scalar unc_UB_`t'_mivsr = unc_UB_`t'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}
	}

		
	/*And the same for the CI bounds*/
	
	if "`att'" != "" {
		forval t=`subminT'/`maxT' {
		forval k=`minT'/`t' {
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar ciLB_`k'_`t'_wc = ciLB_`k'_`t'_wc
			qui ereturn scalar ciUB_`k'_`t'_wc = ciUB_`k'_`t'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar ciLB_`k'_`t'_mts = ciLB_`k'_`t'_mts
			qui ereturn scalar ciUB_`k'_`t'_mts = ciUB_`k'_`t'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar ciLB_`k'_`t'_mts = ciLB_`k'_`t'_mts
			qui ereturn scalar ciUB_`k'_`t'_mts = ciUB_`k'_`t'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar ciLB_`k'_`t'_mtr = ciLB_`k'_`t'_mtr
			qui ereturn scalar ciUB_`k'_`t'_mtr = ciUB_`k'_`t'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar ciLB_`k'_`t'_mtsr = ciLB_`k'_`t'_mtsr
			qui ereturn scalar ciUB_`k'_`t'_mtsr = ciUB_`k'_`t'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
		    qui di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {	
				qui ereturn scalar ciLB_`k'_`t'_mivo = ciLB_`k'_`t'_wcv
				qui ereturn scalar ciUB_`k'_`t'_mivo = ciUB_`k'_`t'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`k'_`t'_mivo = unc_ciLB_`k'_`t'_wcv
					qui ereturn scalar unc_ciUB_`k'_`t'_mivo = unc_ciUB_`k'_`t'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar ciLB_`k'_`t'_mivs = ciLB_`k'_`t'_mtsv
				qui ereturn scalar ciUB_`k'_`t'_mivs = ciUB_`k'_`t'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`k'_`t'_mivs = unc_ciLB_`k'_`t'_mtsv
					qui ereturn scalar unc_ciUB_`k'_`t'_mivs = unc_ciUB_`k'_`t'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar ciLB_`k'_`t'_mivr = ciLB_`k'_`t'_mtrv
				qui ereturn scalar ciUB_`k'_`t'_mivr = ciUB_`k'_`t'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`k'_`t'_mivr = unc_ciLB_`k'_`t'_mtrv
					qui ereturn scalar unc_ciUB_`k'_`t'_mivr = unc_ciUB_`k'_`t'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar ciLB_`k'_`t'_mivsr = ciLB_`k'_`t'_mtsrv
				qui ereturn scalar ciUB_`k'_`t'_mivsr = ciUB_`k'_`t'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`k'_`t'_mivsr = unc_ciLB_`k'_`t'_mtsrv
					qui ereturn scalar unc_ciUB_`k'_`t'_mivsr = unc_ciUB_`k'_`t'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
		}
		}
	}
	else {
		
	foreach t of local Tlevels {
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar ciLB_`t'_wc = ciLB_`t'_wc
			qui ereturn scalar ciUB_`t'_wc = ciUB_`t'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar ciLB_`t'_mts = ciLB_`t'_mts
			qui ereturn scalar ciUB_`t'_mts = ciUB_`t'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar ciLB_`t'_mts = ciLB_`t'_mts
			qui ereturn scalar ciUB_`t'_mts = ciUB_`t'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar ciLB_`t'_mtr = ciLB_`t'_mtr
			qui ereturn scalar ciUB_`t'_mtr = ciUB_`t'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar ciLB_`t'_mtsr = ciLB_`t'_mtsr
			qui ereturn scalar ciUB_`t'_mtsr = ciUB_`t'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
		    qui di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {	
				qui ereturn scalar ciLB_`t'_mivo = ciLB_`t'_wcv
				qui ereturn scalar ciUB_`t'_mivo = ciUB_`t'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`t'_mivo = unc_ciLB_`t'_wcv
					qui ereturn scalar unc_ciUB_`t'_mivo = unc_ciUB_`t'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar ciLB_`t'_mivs = ciLB_`t'_mtsv
				qui ereturn scalar ciUB_`t'_mivs = ciUB_`t'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`t'_mivs = unc_ciLB_`t'_mtsv
					qui ereturn scalar unc_ciUB_`t'_mivs = unc_ciUB_`t'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar ciLB_`t'_mivr = ciLB_`t'_mtrv
				qui ereturn scalar ciUB_`t'_mivr = ciUB_`t'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`t'_mivr = unc_ciLB_`t'_mtrv
					qui ereturn scalar unc_ciUB_`t'_mivr = unc_ciUB_`t'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar ciLB_`t'_mivsr = ciLB_`t'_mtsrv
				qui ereturn scalar ciUB_`t'_mivsr = ciUB_`t'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB_`t'_mivsr = unc_ciLB_`t'_mtsrv
					qui ereturn scalar unc_ciUB_`t'_mivsr = unc_ciUB_`t'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}
	}

	
/*Display results in the Stata window*/

if "`att'" != "" {
		if strpos("${bset}","wc") != 0 {	
		di ""
		di "{col 15}Worst-Case Selection"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
			di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_wc)' "," %12.3f `e(UB_`k'_`t'_wc)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_wc)' "," %12.3f `e(ciUB_`k'_`t'_wc)' ")"
		di ""
		}
		}
	}
	else if strpos("${bset}","wc") == 0 {
		qui di ""
	}
	*************************
	if "`miv'" == "" {
		qui di ""
	}
	else {	
		if strpos("${mbset}","wcv") != 0 {
			di ""
			di "{col 15}MIV-Only"
			di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivo)' "," %12.3f `e(UB_`k'_`t'_mivo)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivo)' "," %12.3f `e(ciUB_`k'_`t'_mivo)' ")"
				di ""
			}
			}
		}
		else if strpos("${mbset}","wcv") == 0 {
			qui di ""
		}
	}
		
	************************
	if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
			di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mts)' "," %12.3f `e(UB_`k'_`t'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mts)' "," %12.3f `e(ciUB_`k'_`t'_mts)' ")"
			di ""
		}
		}
	}
	else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
			di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mts)' "," %12.3f `e(UB_`k'_`t'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mts)' "," %12.3f `e(ciUB_`k'_`t'_mts)' ")"
			di ""
		}
		}
	}
	else if strpos("${bset}","mts") == 0 {
		qui di ""
	}
		
	***********************
	/*The display of bounds involving MTR will depend on which direction is used;
	i.e., which bound is equal to zero*/

	if ("`nmts'" != "" & "`nmtr'" == "") {
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtr)' "," %12.3f `e(UB_`k'_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtr)' "," %12.3f `e(ciUB_`k'_`t'_mtr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtsr)' "," %12.3f `e(UB_`k'_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtsr)' "," %12.3f `e(ciUB_`k'_`t'_mtsr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
			if strpos("${mbset}","mtsv") != 0 {
				di ""
				di "{col 15}MIV+MTS"				
				di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

					di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivs)' "," %12.3f `e(UB_`k'_`t'_mivs)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivs)' "," %12.3f `e(ciUB_`k'_`t'_mivs)' ")"
					di ""
				}
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				di ""
				di "{col 15}MIV+MTR"
				di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

					di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivr)' "," %12.3f `e(UB_`k'_`t'_mivr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivr)' "," %12.3f `e(ciUB_`k'_`t'_mivr)' ")"
					di ""
				}
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				di ""
				di "{col 15}MIV+MTS+MTR"
				di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

					di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivsr)' "," %12.3f `e(UB_`k'_`t'_mivsr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivsr)' "," %12.3f `e(ciUB_`k'_`t'_mivsr)' ")"
					di ""
				}
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}


	******************
	else if ("`nmts'" == "" & "`nmtr'" != "") {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtr)' "," %12.3f min(`e(UB_`k'_`t'_mtr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtr)' "," %12.3f min(`e(ciUB_`k'_`t'_mtr)',0) ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtsr)' "," %12.3f min(`e(UB_`k'_`t'_mtsr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtsr)' "," %12.3f min(`e(ciUB_`k'_`t'_mtsr)',0) ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivs)' "," %12.3f `e(UB_`k'_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivs)' "," %12.3f `e(ciUB_`k'_`t'_mivs)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivr)' "," %12.3f `e(UB_`k'_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivr)' "," %12.3f `e(ciUB_`k'_`t'_mivr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivsr)' "," %12.3f `e(UB_`k'_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivsr)' "," %12.3f `e(ciUB_`k'_`t'_mivsr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

	**********************
		
	else if ("`nmts'" != "" & "`nmtr'" != "") {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtr)' "," %12.3f `e(UB_`k'_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtr)' "," %12.3f `e(ciUB_`k'_`t'_mtr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtsr)' "," %12.3f `e(UB_`k'_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtsr)' "," %12.3f `e(ciUB_`k'_`t'_mtsr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivs)' "," %12.3f `e(UB_`k'_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivs)' "," %12.3f `e(ciUB_`k'_`t'_mivs)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivr)' "," %12.3f `e(UB_`k'_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivr)' "," %12.3f `e(ciUB_`k'_`t'_mivr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivsr)' "," %12.3f `e(UB_`k'_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivsr)' "," %12.3f `e(ciUB_`k'_`t'_mivsr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

	
	else {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtr)' "," %12.3f `e(UB_`k'_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtr)' "," %12.3f `e(ciUB_`k'_`t'_mtr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}	
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {
		
				di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mtsr)' "," %12.3f `e(UB_`k'_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mtsr)' "," %12.3f `e(ciUB_`k'_`t'_mtsr)' ")"
				di ""
			}
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivs)' "," %12.3f `e(UB_`k'_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivs)' "," %12.3f `e(ciUB_`k'_`t'_mivs)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivr)' "," %12.3f `e(UB_`k'_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivr)' "," %12.3f `e(ciUB_`k'_`t'_mivr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

		forval t=`subminT'/`maxT' {
			forval k=`minT'/`t' {

						di "E[Y(`k') | T=`t']" "{col 15}[" %-12.3f `e(LB_`k'_`t'_mivsr)' "," %12.3f `e(UB_`k'_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`k'_`t'_mivsr)' "," %12.3f `e(ciUB_`k'_`t'_mivsr)' ")"
						di ""
					}
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

}

	*********************
	***ATE**************

else {

	if strpos("${bset}","wc") != 0 {	
		di ""
		di "{col 15}Worst-Case Selection"
		di "{hline 30}"

		foreach t of local Tlevels {
			local k=`t'-1
		
			di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_wc)' "," %12.3f `e(UB_`t'_wc)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_wc)' "," %12.3f `e(ciUB_`t'_wc)' ")"
		di ""
		}
	}
	else if strpos("${bset}","wc") == 0 {
		qui di ""
	}
	*************************
	if "`miv'" == "" {
		qui di ""
	}
	else {	
		if strpos("${mbset}","wcv") != 0 {
			di ""
			di "{col 15}MIV-Only"
			di "{hline 30}"

			foreach t of local Tlevels {
				local k=`t'-1

				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivo)' "," %12.3f `e(UB_`t'_mivo)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivo)' "," %12.3f `e(ciUB_`t'_mivo)' ")"
				di ""
			}
		}
		else if strpos("${mbset}","wcv") == 0 {
			qui di ""
		}
	}
		
	************************
	if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		foreach t of local Tlevels {
			local k=`t'-1
		
			di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mts)' "," %12.3f `e(UB_`t'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mts)' "," %12.3f `e(ciUB_`t'_mts)' ")"
			di ""
		}
	}
	else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		foreach t of local Tlevels {
			local k=`t'-1
		
			di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mts)' "," %12.3f `e(UB_`t'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mts)' "," %12.3f `e(ciUB_`t'_mts)' ")"
			di ""
		}
	}
	else if strpos("${bset}","mts") == 0 {
		qui di ""
	}
		
	***********************
	/*The display of bounds involving MTR will depend on which direction is used;
	i.e., which bound is equal to zero*/

	if ("`nmts'" != "" & "`nmtr'" == "") {
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtr)' "," %12.3f `e(UB_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtr)' "," %12.3f `e(ciUB_`t'_mtr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtsr)' "," %12.3f `e(UB_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtsr)' "," %12.3f `e(ciUB_`t'_mtsr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
			if strpos("${mbset}","mtsv") != 0 {
				di ""
				di "{col 15}MIV+MTS"				
				di "{hline 30}"

				foreach t of local Tlevels {
					local k=`t'-1

					di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivs)' "," %12.3f `e(UB_`t'_mivs)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivs)' "," %12.3f `e(ciUB_`t'_mivs)' ")"
					di ""
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				di ""
				di "{col 15}MIV+MTR"
				di "{hline 30}"

				foreach t of local Tlevels {
					local k=`t'-1

					di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivr)' "," %12.3f `e(UB_`t'_mivr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivr)' "," %12.3f `e(ciUB_`t'_mivr)' ")"
					di ""
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				di ""
				di "{col 15}MIV+MTS+MTR"
				di "{hline 30}"

				foreach t of local Tlevels {
					local k=`t'-1

					di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivsr)' "," %12.3f `e(UB_`t'_mivsr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivsr)' "," %12.3f `e(ciUB_`t'_mivsr)' ")"
					di ""
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}


	******************
	else if ("`nmts'" == "" & "`nmtr'" != "") {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtr)' "," %12.3f `e(UB_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtr)' "," %12.3f `e(ciUB_`t'_mtr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtsr)' "," %12.3f `e(UB_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtsr)' "," %12.3f `e(ciUB_`t'_mtsr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivs)' "," %12.3f `e(UB_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivs)' "," %12.3f `e(ciUB_`t'_mivs)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivr)' "," %12.3f `e(UB_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivr)' "," %12.3f `e(ciUB_`t'_mivr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivsr)' "," %12.3f `e(UB_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivsr)' "," %12.3f `e(ciUB_`t'_mivsr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

	**********************
		
	else if ("`nmts'" != "" & "`nmtr'" != "") {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtr)' "," %12.3f `e(UB_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtr)' "," %12.3f `e(ciUB_`t'_mtr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtsr)' "," %12.3f `e(UB_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtsr)' "," %12.3f `e(ciUB_`t'_mtsr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivs)' "," %12.3f `e(UB_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivs)' "," %12.3f `e(ciUB_`t'_mivs)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivr)' "," %12.3f `e(UB_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivr)' "," %12.3f `e(ciUB_`t'_mivr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivsr)' "," %12.3f `e(UB_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivsr)' "," %12.3f `e(ciUB_`t'_mivsr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

	*********************

	else {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtr)' "," %12.3f `e(UB_`t'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtr)' "," %12.3f `e(ciUB_`t'_mtr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtr") == 0 {
			qui di ""
		}	
		
		if strpos("${bset}","mtsr") != 0 {
			di ""
			di "{col 15}MTS+MTR"
			di "{hline 30}"
		
			foreach t of local Tlevels {
				local k=`t'-1
		
				di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mtsr)' "," %12.3f `e(UB_`t'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mtsr)' "," %12.3f `e(ciUB_`t'_mtsr)' ")"
				di ""
			}
		}
		else if strpos("${bset}","mtsr") == 0 {
			qui di ""
		}
		
		if "`miv'" == "" {
			qui di ""
		}
		else {
				if strpos("${mbset}","mtsv") != 0 {
					di ""
					di "{col 15}MIV+MTS"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivs)' "," %12.3f `e(UB_`t'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivs)' "," %12.3f `e(ciUB_`t'_mivs)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					di ""
					di "{col 15}MIV+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivr)' "," %12.3f `e(UB_`t'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivr)' "," %12.3f `e(ciUB_`t'_mivr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					di ""
					di "{col 15}MIV+MTS+MTR"
					di "{hline 30}"

					foreach t of local Tlevels {
						local k=`t'-1

						di "E[Y(`t')]" "{col 15}[" %-12.3f `e(LB_`t'_mivsr)' "," %12.3f `e(UB_`t'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB_`t'_mivsr)' "," %12.3f `e(ciUB_`t'_mivsr)' ")"
						di ""
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}
}

/*Create the code for a LaTeX table of results when the option is used*/

if "`latex'" != "" {
	tempname textable
	file open `textable' using "`latex'.tex", write replace
	file write `textable' "\begin{table}[h!]" _n ///
	"\centering" _n ///
	"\scalebox{0.7125}{" _n ///
	"\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" _n ///
	"\begin{tabular}{l*{8}{c}}" _n ///
	"\toprule" _n ///
					"&\multicolumn{1}{c}{(1)}&\multicolumn{1}{c}{(2)}&\multicolumn{1}{c}{(3)}&\multicolumn{1}{c}{(4)}&\multicolumn{1}{c}{(5)}&\multicolumn{1}{c}{(6)}&\multicolumn{1}{c}{(7)}&\multicolumn{1}{c}{(8)}\\" _n ///
									"&\multicolumn{1}{c}{Worst-Case}&\multicolumn{1}{c}{MIV}&\multicolumn{1}{c}{MTS}&\multicolumn{1}{c}{MTR}&\multicolumn{1}{c}{MTS+MTR}&\multicolumn{1}{c}{MIV+MTS}&\multicolumn{1}{c}{MIV+MTR}&\multicolumn{1}{c}{MIV+MTS+MTR}\\" _n ///
	"\midrule" _n
	
	
	if ("`nmts'" != "" & "`nmtr'" == "") {
	forval t=`subminT'/`maxT' {
		local k=`t'-1
		file write `textable' "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtr)',0)) "," %12.3f (`e(UB`t'`k'_mtr)') "] & [" %-12.3f (max((`e(LB`t'`k'_mtsr)',0)) "," %12.3f (`e(UB`t'`k'_mtsr)') "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (max(`e(LB`t'`k'_mivr)',0)) "," %12.3f (`e(UB`t'`k'_mivr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mivsr)',0)) "," %12.3f (`e(UB`t'`k'_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtsr)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (max(`e(ciLB`t'`k'_mivr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mivsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write `textable' "\midrule" _n ///
			"\addlinespace" _n ///
			"$\mathit{ATE(`t`i'',`s`i'')}$ & [" %-12.3f (`e(LB`t`i''`s`i''_wc)') "," %12.3f (`e(UB`t`i''`s`i''_wc)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivo)') "," %12.3f (`e(UB`t`i''`s`i''_mivo)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mts)') "," %12.3f (`e(UB`t`i''`s`i''_mts)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtsr)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivs)') "," %12.3f (`e(UB`t`i''`s`i''_mivs)') "] & [" %12.3f (max(`e(LB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t`i''`s`i''_wc)') "," %12.3f (`e(ciUB`t`i''`s`i''_wc)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivo)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivo)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mts)') "," %12.3f (`e(ciUB`t`i''`s`i''_mts)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtsr)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivs)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivs)') ") & (" %12.3f (max(`e(ciLB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	}
	}
	
	else if ("`nmts'" == "" & "`nmtr'" != "") {
		forval t=`subminT'/`maxT' {
		local k=`t'-1
		file write `textable' "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (`e(LB`t'`k'_mtr)') "," %12.3f (min(`e(UB`t'`k'_mtr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mtsr)') "," %12.3f (min(`e(UB`t'`k'_mtsr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (`e(LB`t'`k'_mivr)') "," %12.3f (min(`e(UB`t'`k'_mivr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivsr)') "," %12.3f (min(`e(UB`t'`k'_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (`e(ciLB`t'`k'_mtr)') "," %12.3f (min(`e(ciUB`t'`k'_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mtsr)') "," %12.3f (min(`e(ciUB`t'`k'_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (`e(ciLB`t'`k'_mivr)') "," %12.3f (min(`e(ciUB`t'`k'_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivsr)') "," %12.3f (min(`e(ciUB`t'`k'_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write `textable' "\midrule" _n ///
			"\addlinespace" _n ///
			"$\mathit{ATE(`t`i'',`s`i'')}$ & [" %-12.3f (`e(LB`t`i''`s`i''_wc)') "," %12.3f (`e(UB`t`i''`s`i''_wc)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivo)') "," %12.3f (`e(UB`t`i''`s`i''_mivo)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mts)') "," %12.3f (`e(UB`t`i''`s`i''_mts)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mtr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mtr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mtsr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mtsr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivs)') "," %12.3f (`e(UB`t`i''`s`i''_mivs)') "] & [" %12.3f (`e(LB`t`i''`s`i''_mivr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mivr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivsr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t`i''`s`i''_wc)') "," %12.3f (`e(ciUB`t`i''`s`i''_wc)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivo)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivo)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mts)') "," %12.3f (`e(ciUB`t`i''`s`i''_mts)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mtr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mtsr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivs)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivs)') ") & (" %12.3f (`e(ciLB`t`i''`s`i''_mivr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivsr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	}
	}
	
	
	else if ("`nmts'" != "" & "`nmtr'" != "") {
	
	forval t=`subminT'/`maxT' {
		local k=`t'-1
		file write `textable' "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (`e(LB`t'`k'_mtr)') "," %12.3f (min(`e(UB`t'`k'_mtr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mtsr)') "," %12.3f (min(`e(UB`t'`k'_mtsr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (`e(LB`t'`k'_mivr)') "," %12.3f (min(`e(UB`t'`k'_mivr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivsr)') "," %12.3f (min(`e(UB`t'`k'_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (`e(ciLB`t'`k'_mtr)') "," %12.3f (min(`e(ciUB`t'`k'_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mtsr)') "," %12.3f (min(`e(ciUB`t'`k'_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (`e(ciLB`t'`k'_mivr)') "," %12.3f (min(`e(ciUB`t'`k'_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivsr)') "," %12.3f (min(`e(ciUB`t'`k'_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write `textable' "\midrule" _n ///
			"\addlinespace" _n ///
			"$\mathit{ATE(`t`i'',`s`i'')}$ & [" %-12.3f (`e(LB`t`i''`s`i''_wc)') "," %12.3f (`e(UB`t`i''`s`i''_wc)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivo)') "," %12.3f (`e(UB`t`i''`s`i''_mivo)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mts)') "," %12.3f (`e(UB`t`i''`s`i''_mts)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mtr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mtr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mtsr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mtsr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivs)') "," %12.3f (`e(UB`t`i''`s`i''_mivs)') "] & [" %12.3f (`e(LB`t`i''`s`i''_mivr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mivr)',0)) "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivsr)') "," %12.3f (min(`e(UB`t`i''`s`i''_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t`i''`s`i''_wc)') "," %12.3f (`e(ciUB`t`i''`s`i''_wc)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivo)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivo)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mts)') "," %12.3f (`e(ciUB`t`i''`s`i''_mts)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mtr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mtsr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivs)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivs)') ") & (" %12.3f (`e(ciLB`t`i''`s`i''_mivr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivsr)') "," %12.3f (min(`e(ciUB`t`i''`s`i''_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	}
	}
	
	else {
	
	forval t=`subminT'/`maxT' {
		local k=`t'-1
		file write `textable' "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtr)',0)) "," %12.3f (`e(UB`t'`k'_mtr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtsr)',0)) "," %12.3f (`e(UB`t'`k'_mtsr)') "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (max(`e(LB`t'`k'_mivr)',0)) "," %12.3f (`e(UB`t'`k'_mivr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mivsr)',0)) "," %12.3f (`e(UB`t'`k'_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtsr)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (max(`e(ciLB`t'`k'_mivr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mivsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write `textable' "\midrule" _n ///
			"\addlinespace" _n ///
			"$\mathit{ATE(`t`i'',`s`i'')}$ & [" %-12.3f (`e(LB`t`i''`s`i''_wc)') "," %12.3f (`e(UB`t`i''`s`i''_wc)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivo)') "," %12.3f (`e(UB`t`i''`s`i''_mivo)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mts)') "," %12.3f (`e(UB`t`i''`s`i''_mts)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtsr)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivs)') "," %12.3f (`e(UB`t`i''`s`i''_mivs)') "] & [" %12.3f (max(`e(LB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t`i''`s`i''_wc)') "," %12.3f (`e(ciUB`t`i''`s`i''_wc)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivo)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivo)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mts)') "," %12.3f (`e(ciUB`t`i''`s`i''_mts)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtsr)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivs)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivs)') ") & (" %12.3f (max(`e(ciLB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	}
	}
	
	
	file write `textable' "\midrule" _n ///
	"\addlinespace" _n ///
	"Observations & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' \\" _n ///
	"\addlinespace" _n ///
	"\bottomrule" _n ///
	"\multicolumn{8}{l}{\footnotesize  [$\cdot$]: half-median unbiased MIV estimates \ ; \ ($\cdot$): `level'\% CLR confidence intervals}" _n ///
	"\end{tabular}" _n ///
	"}" _n ///
	"\caption{\small Bounds on the ATE: `miv' as MIV, `bins' Bins, with \$K_0=`ymin'$}" _n ///
	"\end{table}"
	

file close `textable'

}

end
