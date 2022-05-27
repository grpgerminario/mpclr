program mpclr, eclass properties(svyb)
	version 10.0
	#delimit ;
	syntax varlist (min=2 max=2 fv) [if] [in] [pw iw],  YMIN(string) YMAX(string) [ATT BOUNDS(string) NMTS NMTR MIV(varname) MIVBOUNDS(string) bins(integer 5) NMIV moret(numlist) mores(numlist) reps(integer 25) seed(integer 22176) level(integer 95) NOISYCLR LATEX(string) NPSU(int -1) SURVEY BRR UNCORRECTED];
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
	
	qui sum `T'
	local minT = r(min)
	local maxT = r(max)
	local subminT = `minT' + 1
	local submaxT = `maxT' - 1
	
	/*Store all levels of the treatment*/
	
	qui levelsof `T', local(Tlevels)
	
	/*If options for "more" ATEs is used, this gathers the info*/

	if "`moret'" != "" {
		local conta: word count `moret'
		local i=1
		foreach x of numlist `moret' {
			local t`i'=`x'
			local ++i
		}
			
		local i=1
		foreach x of numlist `mores' {
			local s`i'=`x'
			local ++i
		}
	}
	
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

		if "`survey'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'-1

						/*This bootstraps the 'mpbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						svy bootstrap _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)

						matrix se`x'`t'`k'_`a' = e(se)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosam_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							svy bootstrap _b, saving(`cosam_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a') 

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
							matrix se`x'`t`i''`s`i''_`a' = e(se)
						}
					}
				}		
			}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'-1

						/*This bootstraps the 'mpbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						svy brr _b, saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosam_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							svy brr _b, saving(`cosam_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a') 

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
						}
					}
				}		
			}
		}


		else {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosa_`x'_`t'_`a'
						local k=`t'-1

						/*This bootstraps the 'mpbounds' subcommand to compute the variance-covariance
						   matrix for the CLR bias correction and inference*/

						di ""
						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosa_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)

						matrix se`x'`t'`k'_`a' = e(se)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosam_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							bootstrap _b, reps(`reps') seed(`seed') saving(`cosam_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a')

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
							matrix se`x'`t`i''`s`i''_`a' = e(se)
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

		if "`survey'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'-1

						di ""
						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						svy bootstrap _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a') 

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)

						matrix se`x'`t'`k'_`a' = e(se)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosavm_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							svy bootstrap _b, saving(`cosavm_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a') 

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
							matrix se`x'`t`i''`s`i''_`a' = e(se)
						}
					}
				}		
			}
		}
		
		else if "`brr'" != "" {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'-1

						di ""
						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						svy brr _b, saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a') 

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosavm_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							svy brr _b, saving(`cosavm_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a') 

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
						}
					}
				}		
			}
		}

		else {
			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					foreach x of local lowup {
						tempfile cosav_`x'_`t'_`a'
						local k=`t'-1

						di ""
						di ""
						di "Begin bootstrap for `x'B`t'`k'_`a' ... "

						bootstrap _b, reps(`reps') seed(`seed') saving(`cosav_`x'_`t'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t') contrs(`k') asmptn(`a')

						matrix b`x'`t'`k'_`a' = e(b)'

						matrix Vb`x'`t'`k'_`a' = e(V)

						matrix se`x'`t'`k'_`a' = e(se)
					}
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						foreach x of local lowup {
							tempfile cosavm_`x'_`i'_`a'

							di ""
							di ""
							di "Begin bootstrap for `x'B`t`i''`s`i''_`a' ... "

							bootstrap _b, reps(`reps') seed(`seed') saving(`cosavm_`x'_`i'_`a''.dta, double replace) noheader notable:  mpbounds `y' `T' if `touse', ymin(`ymin') ymax(`ymax') `att' `nmts' `nmtr' miv(`miv') bins(`bins') `nmiv' moret(`moret') mores(`mores') ul(`x') contrt(`t`i'') contrs(`s`i'') asmptn(`a')

							matrix b`x'`t`i''`s`i''_`a' = e(b)'
							matrix Vb`x'`t`i''`s`i''_`a' = e(V)
							matrix se`x'`t`i''`s`i''_`a' = e(se)
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
					forval t=`subminT'/`maxT' {
						local k=`t'-1
						matselrc bL`t'`k'_`m' unc_bL`t'`k'_`m', r(`uncterms')
						matselrc bL`t'`k'_`m' bL`t'`k'_`m', r(1/`sqterms')
					
						matselrc VbL`t'`k'_`m' unc_VbL`t'`k'_`m', r(`uncterms') c(`uncterms')
						matselrc VbL`t'`k'_`m' VbL`t'`k'_`m', r(1/`sqterms') c(1/`sqterms')
					
						matselrc bU`t'`k'_`m' unc_bU`t'`k'_`m', r(`uncterms')
						matselrc bU`t'`k'_`m' bU`t'`k'_`m', r(1/`sqterms')
					
						matselrc VbU`t'`k'_`m' unc_VbU`t'`k'_`m', r(`uncterms') c(`uncterms')
						matselrc VbU`t'`k'_`m' VbU`t'`k'_`m', r(1/`sqterms') c(1/`sqterms')
					}
					if "`moret'" != "" {		
						forval i=1/`conta' {
							matselrc bL`t`i''`s`i''_`m' unc_bL`t`i''`s`i''_`m', r(`uncterms')
							matselrc bL`t`i''`s`i''_`m' bL`t`i''`s`i''_`m', r(1/`sqterms')
					
							matselrc VbL`t`i''`s`i''_`m' unc_VbL`t`i''`s`i''_`m', r(`uncterms') c(`uncterms')
							matselrc VbL`t`i''`s`i''_`m' VbL`t`i''`s`i''_`m', r(1/`sqterms') c(1/`sqterms')
					
							matselrc bU`t`i''`s`i''_`m' unc_bU`t`i''`s`i''_`m', r(`uncterms')
							matselrc bU`t`i''`s`i''_`m' bU`t`i''`s`i''_`m', r(1/`sqterms')
					
							matselrc VbU`t`i''`s`i''_`m' unc_VbU`t`i''`s`i''_`m', r(`uncterms') c(`uncterms')
							matselrc VbU`t`i''`s`i''_`m' VbU`t`i''`s`i''_`m', r(1/`sqterms') c(1/`sqterms')
						}
					}
				}

				foreach a of local asmpt {
					forval t=`subminT'/`maxT' {
						local k=`t'-1
						matrix bL = bL`t'`k'_`a' 
						matrix VbL = VbL`t'`k'_`a'
						matrix bU = bU`t'`k'_`a' 
						matrix VbU = VbU`t'`k'_`a' 

					/*define the 'current assumption' being used*/

						qui scalar crntasmpt = "`a'"

						if crntasmpt=="wc" {
							di ""
							di "Begin CLR for Worst-Case for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "Begin CLR for MTS for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "Begin CLR for MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "Begin CLR for MTS+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "Begin CLR for MIV-Only for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "Begin CLR for MIV+MTS for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "Begin CLR for MIV+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "Begin CLR for MIV+MTS+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}


					/*Run the CLR command from Flores & Wang;
					   100000 is the number of 'simulations' R in the CLR algorithm*/


						CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar LB`t'`k'_`a' = eLstar[1,1]	
						qui scalar UB`t'`k'_`a' = eUstar[1,1]
						qui scalar ciLB`t'`k'_`a' = eCIL[1,1]
						qui scalar ciUB`t'`k'_`a' = eCIU[1,1]

						if crntasmpt=="wc" {
							di ""
							di "End CLR for Worst-Case for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "End CLR for MTS for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "End CLR for MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "End CLR for MTS+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "End CLR for MIV-Only for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "End CLR for MIV+MTS for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "End CLR for MIV+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "End CLR for MIV+MTS+MTR for ATT(`t',`k')"
							di "{hline 30}"
						}

					}

				/*Repeat the process for the 'more' ATTs if requested*/

					if "`moret'" != "" {
					
					
						forval i=1/`conta' {
							matrix bL = bL`t`i''`s`i''_`a'
							matrix VbL = VbL`t`i''`s`i''_`a'
							matrix bU = bU`t`i''`s`i''_`a'
							matrix VbU = VbU`t`i''`s`i''_`a' 

							if crntasmpt=="wc" {
								di ""
								di "Begin CLR for Worst-Case for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mts" {
								di ""
								di "Begin CLR for MTS for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtr" {
								di ""
								di "Begin CLR for MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsr" {
								di ""
								di "Begin CLR for MTS+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="wcv" {
								di ""
								di "Begin CLR for MIV-Only for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsv" {
								di ""
								di "Begin CLR for MIV+MTS for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtrv" {
								di ""
								di "Begin CLR for MIV+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsrv" {
								di ""
								di "Begin CLR for MIV+MTS+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							CLR bL VbL bU VbU Nob 100000 alpha
							qui scalar LB`t`i''`s`i''_`a' = eLstar[1,1]
							qui scalar UB`t`i''`s`i''_`a' = eUstar[1,1]
							qui scalar ciLB`t`i''`s`i''_`a' = eCIL[1,1]
							qui scalar ciUB`t`i''`s`i''_`a' = eCIU[1,1]

							if crntasmpt=="wc" {
								di ""
								di "End CLR for Worst-Case for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mts" {
								di ""
								di "End CLR for MTS for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtr" {
								di ""
								di "End CLR for MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsr" {
								di ""
								di "End CLR for MTS+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="wcv" {
								di ""
								di "End CLR for MIV-Only for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsv" {
								di ""
								di "End CLR for MIV+MTS for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtrv" {
								di ""
								di "End CLR for MIV+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsrv" {
								di ""
								di "End CLR for MIV+MTS+MTR for ATT(`t`i'',`s`i'')"
								di "{hline 30}"
							}
						}
					}
				}
			}
			
			else {
			
				scalar tterms = 2^(`bins'-1)
				scalar sqterms = tterms^2
				scalar uncterms = sqterms+1
				local sqterms = sqterms
				local uncterms = uncterms
			
				foreach m of global mbset {
					forval t=`subminT'/`maxT' {
						local k=`t'-1
						matselrc bL`t'`k'_`m' unc_bL`t'`k'_`m', r(`uncterms')
						matselrc bL`t'`k'_`m' bL`t'`k'_`m', r(1/`sqterms')
					
						matselrc VbL`t'`k'_`m' unc_VbL`t'`k'_`m', r(`uncterms') c(`uncterms')
						matselrc VbL`t'`k'_`m' VbL`t'`k'_`m', r(1/`sqterms') c(1/`sqterms')
					
						matselrc bU`t'`k'_`m' unc_bU`t'`k'_`m', r(`uncterms')
						matselrc bU`t'`k'_`m' bU`t'`k'_`m', r(1/`sqterms')
					
						matselrc VbU`t'`k'_`m' unc_VbU`t'`k'_`m', r(`uncterms') c(`uncterms')
						matselrc VbU`t'`k'_`m' VbU`t'`k'_`m', r(1/`sqterms') c(1/`sqterms')
					}
					if "`moret'" != "" {		
						forval i=1/`conta' {
							matselrc bL`t`i''`s`i''_`m' unc_bL`t`i''`s`i''_`m', r(`uncterms')
							matselrc bL`t`i''`s`i''_`m' bL`t`i''`s`i''_`m', r(1/`sqterms')
					
							matselrc VbL`t`i''`s`i''_`m' unc_VbL`t`i''`s`i''_`m', r(`uncterms') c(`uncterms')
							matselrc VbL`t`i''`s`i''_`m' VbL`t`i''`s`i''_`m', r(1/`sqterms') c(1/`sqterms')
					
							matselrc bU`t`i''`s`i''_`m' unc_bU`t`i''`s`i''_`m', r(`uncterms')
							matselrc bU`t`i''`s`i''_`m' bU`t`i''`s`i''_`m', r(1/`sqterms')
					
							matselrc VbU`t`i''`s`i''_`m' unc_VbU`t`i''`s`i''_`m', r(`uncterms') c(`uncterms')
							matselrc VbU`t`i''`s`i''_`m' VbU`t`i''`s`i''_`m', r(1/`sqterms') c(1/`sqterms')
						}
					}
				}

				foreach a of local asmpt {
					forval t=`subminT'/`maxT' {
						local k=`t'-1
						matrix bL = bL`t'`k'_`a' 
						matrix VbL = VbL`t'`k'_`a'
						matrix bU = bU`t'`k'_`a' 
						matrix VbU = VbU`t'`k'_`a' 

					/*define the 'current assumption' being used*/

						qui scalar crntasmpt = "`a'"

						if crntasmpt=="wc" {
							di ""
							di "Begin CLR for Worst-Case for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "Begin CLR for MTS for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "Begin CLR for MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "Begin CLR for MTS+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "Begin CLR for MIV-Only for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "Begin CLR for MIV+MTS for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "Begin CLR for MIV+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "Begin CLR for MIV+MTS+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}


					/*Run the CLR command from Flores & Wang;
					   100000 is the number of 'simulations' R in the CLR algorithm*/


						CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar LB`t'`k'_`a' = eLstar[1,1]	
						qui scalar UB`t'`k'_`a' = eUstar[1,1]
						qui scalar ciLB`t'`k'_`a' = eCIL[1,1]
						qui scalar ciUB`t'`k'_`a' = eCIU[1,1]

						if crntasmpt=="wc" {
							di ""
							di "End CLR for Worst-Case for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mts" {
							di ""
							di "End CLR for MTS for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtr" {
							di ""
							di "End CLR for MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsr" {
							di ""
							di "End CLR for MTS+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="wcv" {
							di ""
							di "End CLR for MIV-Only for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsv" {
							di ""
							di "End CLR for MIV+MTS for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtrv" {
							di ""
							di "End CLR for MIV+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

						else if crntasmpt=="mtsrv" {
							di ""
							di "End CLR for MIV+MTS+MTR for ATE(`t',`k')"
							di "{hline 30}"
						}

					}

				/*Repeat the process for the 'more' ATEs if requested*/

					if "`moret'" != "" {
					
					
						forval i=1/`conta' {
							matrix bL = bL`t`i''`s`i''_`a'
							matrix VbL = VbL`t`i''`s`i''_`a'
							matrix bU = bU`t`i''`s`i''_`a'
							matrix VbU = VbU`t`i''`s`i''_`a' 

							if crntasmpt=="wc" {
								di ""
								di "Begin CLR for Worst-Case for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mts" {
								di ""
								di "Begin CLR for MTS for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtr" {
								di ""
								di "Begin CLR for MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsr" {
								di ""
								di "Begin CLR for MTS+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="wcv" {
								di ""
								di "Begin CLR for MIV-Only for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsv" {
								di ""
								di "Begin CLR for MIV+MTS for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtrv" {
								di ""
								di "Begin CLR for MIV+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsrv" {
								di ""
								di "Begin CLR for MIV+MTS+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							CLR bL VbL bU VbU Nob 100000 alpha
							qui scalar LB`t`i''`s`i''_`a' = eLstar[1,1]
							qui scalar UB`t`i''`s`i''_`a' = eUstar[1,1]
							qui scalar ciLB`t`i''`s`i''_`a' = eCIL[1,1]
							qui scalar ciUB`t`i''`s`i''_`a' = eCIU[1,1]

							if crntasmpt=="wc" {
								di ""
								di "End CLR for Worst-Case for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mts" {
								di ""
								di "End CLR for MTS for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtr" {
								di ""
								di "End CLR for MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsr" {
								di ""
								di "End CLR for MTS+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="wcv" {
								di ""
								di "End CLR for MIV-Only for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsv" {
								di ""
								di "End CLR for MIV+MTS for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtrv" {
								di ""
								di "End CLR for MIV+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}

							else if crntasmpt=="mtsrv" {
								di ""
								di "End CLR for MIV+MTS+MTR for ATE(`t`i'',`s`i'')"
								di "{hline 30}"
							}
						}
					}
				}

			}
		}

		/*If 'noisyclr' is NOT requested, then everything is done quietly...
			   this makes the Viewer a lot more readable; if noisy is requested, we
		   suggest starting a log file before running the 'mpclr' command*/

		else {	
			
			scalar tterms = 2^(`bins'-1)
			scalar sqterms = tterms^2
			*scalar uncterms = sqterms+1
			*local sqterms = sqterms
			*local uncterms = uncterms
			
			foreach m of global mbset {
				forval t=`subminT'/`maxT' {
					local k=`t'-1
					
					scalar unctermsL = rowsof(bL`t'`k'_`m')
					scalar ftermsL = unctermsL-1
					local unctermsL = unctermsL
					local ftermsL = ftermsL
					
					matselrc bL`t'`k'_`m' unc_bL`t'`k'_`m', r(`unctermsL')
					matselrc bL`t'`k'_`m' bL`t'`k'_`m', r(1/`ftermsL')
					
					matselrc VbL`t'`k'_`m' unc_VbL`t'`k'_`m', r(`unctermsL') c(`unctermsL')
					matselrc VbL`t'`k'_`m' VbL`t'`k'_`m', r(1/`ftermsL') c(1/`ftermsL')
					
					scalar unctermsU = rowsof(bU`t'`k'_`m')
					scalar ftermsU = unctermsU-1
					local unctermsU = unctermsU
					local ftermsU = ftermsU
					
					matselrc bU`t'`k'_`m' unc_bU`t'`k'_`m', r(`unctermsU')
					matselrc bU`t'`k'_`m' bU`t'`k'_`m', r(1/`ftermsU')
					
					matselrc VbU`t'`k'_`m' unc_VbU`t'`k'_`m', r(`unctermsU') c(`unctermsU')
					matselrc VbU`t'`k'_`m' VbU`t'`k'_`m', r(1/`ftermsU') c(1/`ftermsU')
				}
				if "`moret'" != "" {		
					forval i=1/`conta' {
						
						scalar unctermsL = rowsof(bL`t`i''`s`i''_`m')
						scalar ftermsL = unctermsL-1
						local unctermsL = unctermsL
						local ftermsL = ftermsL
						
						matselrc bL`t`i''`s`i''_`m' unc_bL`t`i''`s`i''_`m', r(`unctermsL')
						matselrc bL`t`i''`s`i''_`m' bL`t`i''`s`i''_`m', r(1/`ftermsL')
					
						matselrc VbL`t`i''`s`i''_`m' unc_VbL`t`i''`s`i''_`m', r(`unctermsL') c(`unctermsL')
						matselrc VbL`t`i''`s`i''_`m' VbL`t`i''`s`i''_`m', r(1/`ftermsL') c(1/`ftermsL')
						
						scalar unctermsU = rowsof(bU`t`i''`s`i''_`m')
						scalar ftermsU = unctermsU-1
						local unctermsU = unctermsU
						local ftermsU = ftermsU
					
						matselrc bU`t`i''`s`i''_`m' unc_bU`t`i''`s`i''_`m', r(`unctermsU')
						matselrc bU`t`i''`s`i''_`m' bU`t`i''`s`i''_`m', r(1/`ftermsU')
					
						matselrc VbU`t`i''`s`i''_`m' unc_VbU`t`i''`s`i''_`m', r(`unctermsU') c(`unctermsU')
						matselrc VbU`t`i''`s`i''_`m' VbU`t`i''`s`i''_`m', r(1/`ftermsU') c(1/`ftermsU')
					}
				}
			}

			foreach a of local asmpt {
				forval t=`subminT'/`maxT' {
					local k=`t'-1
					matrix bL = bL`t'`k'_`a' 
					matrix VbL = VbL`t'`k'_`a'
					matrix bU = bU`t'`k'_`a' 
					matrix VbU = VbU`t'`k'_`a' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar LB`t'`k'_`a' = eLstar[1,1]	
					qui scalar UB`t'`k'_`a' = eUstar[1,1]
					qui scalar ciLB`t'`k'_`a' = eCIL[1,1]
					qui scalar ciUB`t'`k'_`a' = eCIU[1,1]
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						matrix bL = bL`t`i''`s`i''_`a'
						matrix VbL = VbL`t`i''`s`i''_`a'
						matrix bU = bU`t`i''`s`i''_`a'
						matrix VbU = VbU`t`i''`s`i''_`a' 

						qui CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar LB`t`i''`s`i''_`a' = eLstar[1,1]
						qui scalar UB`t`i''`s`i''_`a' = eUstar[1,1]
						qui scalar ciLB`t`i''`s`i''_`a' = eCIL[1,1]
						qui scalar ciUB`t`i''`s`i''_`a' = eCIU[1,1]
					}
				}
			}

		} /*close else for noisyclr*/
		
		if "`uncorrected'" != "" {
			foreach m of global mbset {
				forval t=`subminT'/`maxT' {
					local k=`t'-1
					matrix bL = unc_bL`t'`k'_`m' 
					matrix VbL = unc_VbL`t'`k'_`m'
					matrix bU = unc_bU`t'`k'_`m' 
					matrix VbU = unc_VbU`t'`k'_`m' 

					qui CLR bL VbL bU VbU Nob 100000 alpha
					qui scalar unc_LB`t'`k'_`m' = eLstar[1,1]	
					qui scalar unc_UB`t'`k'_`m' = eUstar[1,1]
					qui scalar unc_ciLB`t'`k'_`m' = eCIL[1,1]
					qui scalar unc_ciUB`t'`k'_`m' = eCIU[1,1]
				}
				if "`moret'" != "" {
					forval i=1/`conta' {
						matrix bL = unc_bL`t`i''`s`i''_`m'
						matrix VbL = unc_VbL`t`i''`s`i''_`m'
						matrix bU = unc_bU`t`i''`s`i''_`m'
						matrix VbU = unc_VbU`t`i''`s`i''_`m' 

						qui CLR bL VbL bU VbU Nob 100000 alpha
						qui scalar unc_LB`t`i''`s`i''_`m' = eLstar[1,1]
						qui scalar unc_UB`t`i''`s`i''_`m' = eUstar[1,1]
						qui scalar unc_ciLB`t`i''`s`i''_`m' = eCIL[1,1]
						qui scalar unc_ciUB`t`i''`s`i''_`m' = eCIU[1,1]
					}
				}
			}
		}

	}
	

	/*re-store the bounds as estimates*/
	
	forval t=`subminT'/`maxT' {
		local k=`t'-1
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar LB`t'`k'_wc = LB`t'`k'_wc
			qui ereturn scalar UB`t'`k'_wc = UB`t'`k'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar LB`t'`k'_mts = LB`t'`k'_mts
			qui ereturn scalar UB`t'`k'_mts = UB`t'`k'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar LB`t'`k'_mts = LB`t'`k'_mts
			qui ereturn scalar UB`t'`k'_mts = UB`t'`k'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar LB`t'`k'_mtr = LB`t'`k'_mtr
			qui ereturn scalar UB`t'`k'_mtr = UB`t'`k'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar LB`t'`k'_mtsr = LB`t'`k'_mtsr
			qui ereturn scalar UB`t'`k'_mtsr = UB`t'`k'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
			di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {
				qui ereturn scalar LB`t'`k'_mivo = LB`t'`k'_wcv
				qui ereturn scalar UB`t'`k'_mivo = UB`t'`k'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB`t'`k'_mivo = unc_LB`t'`k'_wcv
					qui ereturn scalar unc_UB`t'`k'_mivo = unc_UB`t'`k'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar LB`t'`k'_mivs = LB`t'`k'_mtsv
				qui ereturn scalar UB`t'`k'_mivs = UB`t'`k'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB`t'`k'_mivs = unc_LB`t'`k'_mtsv
					qui ereturn scalar unc_UB`t'`k'_mivs = unc_UB`t'`k'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar LB`t'`k'_mivr = LB`t'`k'_mtrv
				qui ereturn scalar UB`t'`k'_mivr = UB`t'`k'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB`t'`k'_mivr = unc_LB`t'`k'_mtrv
					qui ereturn scalar unc_UB`t'`k'_mivr = unc_UB`t'`k'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar LB`t'`k'_mivsr = LB`t'`k'_mtsrv
				qui ereturn scalar UB`t'`k'_mivsr = UB`t'`k'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_LB`t'`k'_mivsr = unc_LB`t'`k'_mtsrv
					qui ereturn scalar unc_UB`t'`k'_mivsr = unc_UB`t'`k'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}
	
	if "`moret'" != "" {
		forval i=1/`conta' {
			
			if strpos("${bset}","wc") != 0 {
				qui ereturn scalar LB`t`i''`s`i''_wc = LB`t`i''`s`i''_wc
				qui ereturn scalar UB`t`i''`s`i''_wc = UB`t`i''`s`i''_wc
			}
			else if strpos("${bset}","wc") == 0 {
				qui di ""
			}

			if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
				qui ereturn scalar LB`t`i''`s`i''_mts = LB`t`i''`s`i''_mts
				qui ereturn scalar UB`t`i''`s`i''_mts = UB`t`i''`s`i''_mts
			}
			else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
				qui ereturn scalar LB`t`i''`s`i''_mts = LB`t`i''`s`i''_mts
				qui ereturn scalar UB`t`i''`s`i''_mts = UB`t`i''`s`i''_mts
			}
			else if strpos("${bset}","mts") == 0 {
				qui di ""
			}

			if strpos("${bset}","mtr") != 0 {
				qui ereturn scalar LB`t`i''`s`i''_mtr = LB`t`i''`s`i''_mtr
				qui ereturn scalar UB`t`i''`s`i''_mtr = UB`t`i''`s`i''_mtr
			}
			else if strpos("${bset}","mtr") == 0 {
				qui di ""
			}

			if strpos("${bset}","mtsr") != 0 {
				qui ereturn scalar LB`t`i''`s`i''_mtsr = LB`t`i''`s`i''_mtsr
				qui ereturn scalar UB`t`i''`s`i''_mtsr = UB`t`i''`s`i''_mtsr
			}
			else if strpos("${bset}","mtsr") == 0 {
				qui di ""
			}



			if "`miv'" == "" {
				qui di ""
			}
			else {
				if strpos("${mbset}","wcv") != 0 {
					qui ereturn scalar LB`t`i''`s`i''_mivo = LB`t`i''`s`i''_wcv
					qui ereturn scalar UB`t`i''`s`i''_mivo = UB`t`i''`s`i''_wcv
				
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_LB`t`i''`s`i''_mivo = unc_LB`t`i''`s`i''_wcv
						qui ereturn scalar unc_UB`t`i''`s`i''_mivo = unc_UB`t`i''`s`i''_wcv
					}
				}
				else if strpos("${mbset}","wcv") == 0 {
					qui di ""
				}
				
				if strpos("${mbset}","mtsv") != 0 {
					qui ereturn scalar LB`t`i''`s`i''_mivs = LB`t`i''`s`i''_mtsv
					qui ereturn scalar UB`t`i''`s`i''_mivs = UB`t`i''`s`i''_mtsv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_LB`t`i''`s`i''_mivs = unc_LB`t`i''`s`i''_mtsv
						qui ereturn scalar unc_UB`t`i''`s`i''_mivs = unc_UB`t`i''`s`i''_mtsv
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					qui ereturn scalar LB`t`i''`s`i''_mivr = LB`t`i''`s`i''_mtrv
					qui ereturn scalar UB`t`i''`s`i''_mivr = UB`t`i''`s`i''_mtrv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_LB`t`i''`s`i''_mivr = unc_LB`t`i''`s`i''_mtrv
						qui ereturn scalar unc_UB`t`i''`s`i''_mivr = unc_UB`t`i''`s`i''_mtrv
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					qui ereturn scalar LB`t`i''`s`i''_mivsr = LB`t`i''`s`i''_mtsrv
					qui ereturn scalar UB`t`i''`s`i''_mivsr = UB`t`i''`s`i''_mtsrv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_LB`t`i''`s`i''_mivsr = unc_LB`t`i''`s`i''_mtsrv
						qui ereturn scalar unc_UB`t`i''`s`i''_mivsr = unc_UB`t`i''`s`i''_mtsrv
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
			}
		}
	}
		
	/*And the same for the CI bounds*/
	
	forval t=`subminT'/`maxT' {
		local k=`t'-1
		
		if strpos("${bset}","wc") != 0 {
			qui ereturn scalar ciLB`t'`k'_wc = ciLB`t'`k'_wc
			qui ereturn scalar ciUB`t'`k'_wc = ciUB`t'`k'_wc
		}
		else if strpos("${bset}","wc") == 0 {
		    qui di ""
		}
		
		if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
			qui ereturn scalar ciLB`t'`k'_mts = ciLB`t'`k'_mts
			qui ereturn scalar ciUB`t'`k'_mts = ciUB`t'`k'_mts
		}
		else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
			qui ereturn scalar ciLB`t'`k'_mts = ciLB`t'`k'_mts
			qui ereturn scalar ciUB`t'`k'_mts = ciUB`t'`k'_mts
		}
		else if strpos("${bset}","mts") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtr") != 0 {
			qui ereturn scalar ciLB`t'`k'_mtr = ciLB`t'`k'_mtr
			qui ereturn scalar ciUB`t'`k'_mtr = ciUB`t'`k'_mtr
		}
		else if strpos("${bset}","mtr") == 0 {
		    qui di ""
		}
		
		if strpos("${bset}","mtsr") != 0 {
			qui ereturn scalar ciLB`t'`k'_mtsr = ciLB`t'`k'_mtsr
			qui ereturn scalar ciUB`t'`k'_mtsr = ciUB`t'`k'_mtsr
		}
		else if strpos("${bset}","mtsr") == 0 {
		    qui di ""
		}
		
		if "`miv'" == "" {
		    qui di ""
		}
		else {
			if strpos("${mbset}","wcv") != 0 {	
				qui ereturn scalar ciLB`t'`k'_mivo = ciLB`t'`k'_wcv
				qui ereturn scalar ciUB`t'`k'_mivo = ciUB`t'`k'_wcv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB`t'`k'_mivo = unc_ciLB`t'`k'_wcv
					qui ereturn scalar unc_ciUB`t'`k'_mivo = unc_ciUB`t'`k'_wcv
				}
			}
			else if strpos("${mbset}","wcv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsv") != 0 {
				qui ereturn scalar ciLB`t'`k'_mivs = ciLB`t'`k'_mtsv
				qui ereturn scalar ciUB`t'`k'_mivs = ciUB`t'`k'_mtsv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB`t'`k'_mivs = unc_ciLB`t'`k'_mtsv
					qui ereturn scalar unc_ciUB`t'`k'_mivs = unc_ciUB`t'`k'_mtsv
				}
			}
			else if strpos("${mbset}","mtsv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtrv") != 0 {
				qui ereturn scalar ciLB`t'`k'_mivr = ciLB`t'`k'_mtrv
				qui ereturn scalar ciUB`t'`k'_mivr = ciUB`t'`k'_mtrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB`t'`k'_mivr = unc_ciLB`t'`k'_mtrv
					qui ereturn scalar unc_ciUB`t'`k'_mivr = unc_ciUB`t'`k'_mtrv
				}
			}
			else if strpos("${mbset}","mtrv") == 0 {
				qui di ""
			}

			if strpos("${mbset}","mtsrv") != 0 {
				qui ereturn scalar ciLB`t'`k'_mivsr = ciLB`t'`k'_mtsrv
				qui ereturn scalar ciUB`t'`k'_mivsr = ciUB`t'`k'_mtsrv
				
				if "`uncorrected'" != "" {
					qui ereturn scalar unc_ciLB`t'`k'_mivsr = unc_ciLB`t'`k'_mtsrv
					qui ereturn scalar unc_ciUB`t'`k'_mivsr = unc_ciUB`t'`k'_mtsrv
				}
			}
			else if strpos("${mbset}","mtsrv") == 0 {
				qui di ""
			}
		}
	}
	
	if "`moret'" != "" {
		forval i=1/`conta' {
			
			if strpos("${bset}","wc") != 0 {
				qui ereturn scalar ciLB`t`i''`s`i''_wc = ciLB`t`i''`s`i''_wc
				qui ereturn scalar ciUB`t`i''`s`i''_wc = ciUB`t`i''`s`i''_wc
			}
			else if strpos("${bset}","wc") == 0 {
				qui di ""
			}

			if (strpos("${bset}","mts") != 0 & strpos("${bset}","mts") < strpos("${bset}","mtsr")) {
				qui ereturn scalar ciLB`t`i''`s`i''_mts = ciLB`t`i''`s`i''_mts
				qui ereturn scalar ciUB`t`i''`s`i''_mts = ciUB`t`i''`s`i''_mts
			}
			else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {
				qui ereturn scalar ciLB`t`i''`s`i''_mts = ciLB`t`i''`s`i''_mts
				qui ereturn scalar ciUB`t`i''`s`i''_mts = ciUB`t`i''`s`i''_mts
			}
			else if strpos("${bset}","mts") == 0 {
				qui di ""
			}

			if strpos("${bset}","mtr") != 0 {
				qui ereturn scalar ciLB`t`i''`s`i''_mtr = ciLB`t`i''`s`i''_mtr
				qui ereturn scalar ciUB`t`i''`s`i''_mtr = ciUB`t`i''`s`i''_mtr
			}
			else if strpos("${bset}","mtr") == 0 {
				qui di ""
			}

			if strpos("${bset}","mtsr") != 0 {
				qui ereturn scalar ciLB`t`i''`s`i''_mtsr = ciLB`t`i''`s`i''_mtsr
				qui ereturn scalar ciUB`t`i''`s`i''_mtsr = ciUB`t`i''`s`i''_mtsr
			}
			else if strpos("${bset}","mtsr") == 0 {
				qui di ""
			}

			if "`miv'" == "" {
			    qui di ""
			}
			else {
				if strpos("${mbset}","wcv") != 0 {
					qui ereturn scalar ciLB`t`i''`s`i''_mivo = ciLB`t`i''`s`i''_wcv
					qui ereturn scalar ciUB`t`i''`s`i''_mivo = ciUB`t`i''`s`i''_wcv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_ciLB`t`i''`s`i''_mivo = unc_ciLB`t`i''`s`i''_wcv
						qui ereturn scalar unc_ciUB`t`i''`s`i''_mivo = unc_ciUB`t`i''`s`i''_wcv
					}
				}
				else if strpos("${mbset}","wcv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsv") != 0 {
					qui ereturn scalar ciLB`t`i''`s`i''_mivs = ciLB`t`i''`s`i''_mtsv
					qui ereturn scalar ciUB`t`i''`s`i''_mivs = ciUB`t`i''`s`i''_mtsv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_ciLB`t`i''`s`i''_mivs = unc_ciLB`t`i''`s`i''_mtsv
						qui ereturn scalar unc_ciUB`t`i''`s`i''_mivs = unc_ciUB`t`i''`s`i''_mtsv
					}
				}
				else if strpos("${mbset}","mtsv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtrv") != 0 {
					qui ereturn scalar ciLB`t`i''`s`i''_mivr = ciLB`t`i''`s`i''_mtrv
					qui ereturn scalar ciUB`t`i''`s`i''_mivr = ciUB`t`i''`s`i''_mtrv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_ciLB`t`i''`s`i''_mivr = unc_ciLB`t`i''`s`i''_mtrv
						qui ereturn scalar unc_ciUB`t`i''`s`i''_mivr = unc_ciUB`t`i''`s`i''_mtrv
					}
				}
				else if strpos("${mbset}","mtrv") == 0 {
					qui di ""
				}

				if strpos("${mbset}","mtsrv") != 0 {
					qui ereturn scalar ciLB`t`i''`s`i''_mivsr = ciLB`t`i''`s`i''_mtsrv
					qui ereturn scalar ciUB`t`i''`s`i''_mivsr = ciUB`t`i''`s`i''_mtsrv
					
					if "`uncorrected'" != "" {
						qui ereturn scalar unc_ciLB`t`i''`s`i''_mivsr = unc_ciLB`t`i''`s`i''_mtsrv
						qui ereturn scalar unc_ciUB`t`i''`s`i''_mivsr = unc_ciUB`t`i''`s`i''_mtsrv
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
			local k=`t'-1
		
			di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_wc)' "," %12.3f `e(UB`t'`k'_wc)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_wc)' "," %12.3f `e(ciUB`t'`k'_wc)' ")"
		di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_wc)' "," %12.3f `e(UB`t`i''`s`i''_wc)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_wc)' "," %12.3f `e(ciUB`t`i''`s`i''_wc)' ")"
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
				local k=`t'-1

				di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivo)' "," %12.3f `e(UB`t'`k'_mivo)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivo)' "," %12.3f `e(ciUB`t'`k'_mivo)' ")"
				di ""
			}

			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivo)' "," %12.3f `e(UB`t`i''`s`i''_mivo)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivo)' "," %12.3f `e(ciUB`t`i''`s`i''_mivo)' ")"
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
			local k=`t'-1
		
			di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mts)' "," %12.3f `e(UB`t'`k'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mts)' "," %12.3f `e(ciUB`t'`k'_mts)' ")"
			di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mts)' "," %12.3f `e(UB`t`i''`s`i''_mts)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mts)' "," %12.3f `e(ciUB`t`i''`s`i''_mts)' ")"
				di ""
			}
		}
	}
	else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			local k=`t'-1
		
			di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mts)' "," %12.3f `e(UB`t'`k'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mts)' "," %12.3f `e(ciUB`t'`k'_mts)' ")"
			di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mts)' "," %12.3f `e(UB`t`i''`s`i''_mts)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mts)' "," %12.3f `e(ciUB`t`i''`s`i''_mts)' ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtr)',0) "," %12.3f `e(UB`t'`k'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtr)',0) "," %12.3f `e(ciUB`t'`k'_mtr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtr)' ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtsr)',0) "," %12.3f `e(UB`t'`k'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtsr)',0) "," %12.3f `e(ciUB`t'`k'_mtsr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtsr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtsr)' ")"
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
					local k=`t'-1

					di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
					local k=`t'-1

					di "ATT(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivr)',0) "," %12.3f `e(UB`t'`k'_mivr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivr)',0) "," %12.3f `e(ciUB`t'`k'_mivr)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivr)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivr)' ")"
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
					local k=`t'-1

					di "ATT(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivsr)',0) "," %12.3f `e(UB`t'`k'_mivsr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivsr)',0) "," %12.3f `e(ciUB`t'`k'_mivsr)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivsr)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivsr)' ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtr)' "," %12.3f min(`e(UB`t'`k'_mtr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtr)' "," %12.3f min(`e(ciUB`t'`k'_mtr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtr)',0) ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtsr)' "," %12.3f min(`e(UB`t'`k'_mtsr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtsr)' "," %12.3f min(`e(ciUB`t'`k'_mtsr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtsr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtsr)',0) ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivr)' "," %12.3f min(`e(UB`t'`k'_mivr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivr)' "," %12.3f min(`e(ciUB`t'`k'_mivr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivr)',0) ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivsr)' "," %12.3f min(`e(UB`t'`k'_mivsr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivsr)' "," %12.3f min(`e(ciUB`t'`k'_mivsr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivsr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivsr)',0) ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtr)' "," %12.3f min(`e(UB`t'`k'_mtr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtr)' "," %12.3f min(`e(ciUB`t'`k'_mtr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtr)',0) ")"
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
				local k=`t'-1
		
				di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtsr)' "," %12.3f min(`e(UB`t'`k'_mtsr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtsr)' "," %12.3f min(`e(ciUB`t'`k'_mtsr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtsr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtsr)',0) ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivr)' "," %12.3f min(`e(UB`t'`k'_mivr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivr)' "," %12.3f min(`e(ciUB`t'`k'_mivr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivr)',0) ")"
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
						local k=`t'-1

						di "ATT(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivsr)' "," %12.3f min(`e(UB`t'`k'_mivsr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivsr)' "," %12.3f min(`e(ciUB`t'`k'_mivsr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATT(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivsr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivsr)',0) ")"
							di ""
						}
					}
				}
				else if strpos("${mbset}","mtsrv") == 0 {
					qui di ""
				}
		}
	}

	*********************
	***ATE**************
	
	else {
		
		if strpos("${bset}","mtr") != 0 {
			di ""
			di "{col 15}MTR"
			di "{hline 30}"
		
			forval t=`subminT'/`maxT' {
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtr)',0) "," %12.3f `e(UB`t'`k'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtr)',0) "," %12.3f `e(ciUB`t'`k'_mtr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtr)' ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtsr)',0) "," %12.3f `e(UB`t'`k'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtsr)',0) "," %12.3f `e(ciUB`t'`k'_mtsr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtsr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtsr)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivr)',0) "," %12.3f `e(UB`t'`k'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivr)',0) "," %12.3f `e(ciUB`t'`k'_mivr)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivr)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivr)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivsr)',0) "," %12.3f `e(UB`t'`k'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivsr)',0) "," %12.3f `e(ciUB`t'`k'_mivsr)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivsr)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivsr)' ")"
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

else {

	if strpos("${bset}","wc") != 0 {	
		di ""
		di "{col 15}Worst-Case Selection"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			local k=`t'-1
		
			di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_wc)' "," %12.3f `e(UB`t'`k'_wc)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_wc)' "," %12.3f `e(ciUB`t'`k'_wc)' ")"
		di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_wc)' "," %12.3f `e(UB`t`i''`s`i''_wc)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_wc)' "," %12.3f `e(ciUB`t`i''`s`i''_wc)' ")"
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
				local k=`t'-1

				di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivo)' "," %12.3f `e(UB`t'`k'_mivo)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivo)' "," %12.3f `e(ciUB`t'`k'_mivo)' ")"
				di ""
			}

			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivo)' "," %12.3f `e(UB`t`i''`s`i''_mivo)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivo)' "," %12.3f `e(ciUB`t`i''`s`i''_mivo)' ")"
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
			local k=`t'-1
		
			di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mts)' "," %12.3f `e(UB`t'`k'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mts)' "," %12.3f `e(ciUB`t'`k'_mts)' ")"
			di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mts)' "," %12.3f `e(UB`t`i''`s`i''_mts)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mts)' "," %12.3f `e(ciUB`t`i''`s`i''_mts)' ")"
				di ""
			}
		}
	}
	else if (strpos("${bset}","mts") != 0 & strpos("${bset}","mtsr") == 0) {	
		di ""
		di "{col 15}MTS"
		di "{hline 30}"

		forval t=`subminT'/`maxT' {
			local k=`t'-1
		
			di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mts)' "," %12.3f `e(UB`t'`k'_mts)' "]" 
			di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mts)' "," %12.3f `e(ciUB`t'`k'_mts)' ")"
			di ""
		}

		if "`moret'" != "" {
			forval i=1/`conta' {
				di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mts)' "," %12.3f `e(UB`t`i''`s`i''_mts)' "]"
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mts)' "," %12.3f `e(ciUB`t`i''`s`i''_mts)' ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtr)',0) "," %12.3f `e(UB`t'`k'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtr)',0) "," %12.3f `e(ciUB`t'`k'_mtr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtr)' ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtsr)',0) "," %12.3f `e(UB`t'`k'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtsr)',0) "," %12.3f `e(ciUB`t'`k'_mtsr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtsr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtsr)' ")"
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
					local k=`t'-1

					di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
					local k=`t'-1

					di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivr)',0) "," %12.3f `e(UB`t'`k'_mivr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivr)',0) "," %12.3f `e(ciUB`t'`k'_mivr)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivr)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivr)' ")"
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
					local k=`t'-1

					di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivsr)',0) "," %12.3f `e(UB`t'`k'_mivsr)' "]" 
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivsr)',0) "," %12.3f `e(ciUB`t'`k'_mivsr)' ")"
					di ""
				}

				if "`moret'" != "" {
					forval i=1/`conta' {
						di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivsr)' "]"
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivsr)' ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtr)' "," %12.3f min(`e(UB`t'`k'_mtr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtr)' "," %12.3f min(`e(ciUB`t'`k'_mtr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtr)',0) ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtsr)' "," %12.3f min(`e(UB`t'`k'_mtsr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtsr)' "," %12.3f min(`e(ciUB`t'`k'_mtsr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtsr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtsr)',0) ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivr)' "," %12.3f min(`e(UB`t'`k'_mivr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivr)' "," %12.3f min(`e(ciUB`t'`k'_mivr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivr)',0) ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivsr)' "," %12.3f min(`e(UB`t'`k'_mivsr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivsr)' "," %12.3f min(`e(ciUB`t'`k'_mivsr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivsr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivsr)',0) ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtr)' "," %12.3f min(`e(UB`t'`k'_mtr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtr)' "," %12.3f min(`e(ciUB`t'`k'_mtr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtr)',0) ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mtsr)' "," %12.3f min(`e(UB`t'`k'_mtsr)',0) "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mtsr)' "," %12.3f min(`e(ciUB`t'`k'_mtsr)',0) ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mtsr)',0) "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mtsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mtsr)',0) ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivr)' "," %12.3f min(`e(UB`t'`k'_mivr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivr)' "," %12.3f min(`e(ciUB`t'`k'_mivr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivr)',0) ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivsr)' "," %12.3f min(`e(UB`t'`k'_mivsr)',0) "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivsr)' "," %12.3f min(`e(ciUB`t'`k'_mivsr)',0) ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(UB`t`i''`s`i''_mivsr)',0) "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivsr)' "," %12.3f min(`e(ciUB`t`i''`s`i''_mivsr)',0) ")"
							di ""
						}
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
		
			forval t=`subminT'/`maxT' {
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtr)',0) "," %12.3f `e(UB`t'`k'_mtr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtr)',0) "," %12.3f `e(ciUB`t'`k'_mtr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtr)' ")"
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
				local k=`t'-1
		
				di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mtsr)',0) "," %12.3f `e(UB`t'`k'_mtsr)' "]" 
				di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mtsr)',0) "," %12.3f `e(ciUB`t'`k'_mtsr)' ")"
				di ""
			}
		
			if "`moret'" != "" {
				forval i=1/`conta' {
					di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mtsr)' "]"
					di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mtsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mtsr)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f `e(LB`t'`k'_mivs)' "," %12.3f `e(UB`t'`k'_mivs)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t'`k'_mivs)' "," %12.3f `e(ciUB`t'`k'_mivs)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f `e(LB`t`i''`s`i''_mivs)' "," %12.3f `e(UB`t`i''`s`i''_mivs)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f `e(ciLB`t`i''`s`i''_mivs)' "," %12.3f `e(ciUB`t`i''`s`i''_mivs)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivr)',0) "," %12.3f `e(UB`t'`k'_mivr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivr)',0) "," %12.3f `e(ciUB`t'`k'_mivr)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivr)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivr)' ")"
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
						local k=`t'-1

						di "ATE(`t',`k')" "{col 15}[" %-12.3f max(`e(LB`t'`k'_mivsr)',0) "," %12.3f `e(UB`t'`k'_mivsr)' "]" 
						di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t'`k'_mivsr)',0) "," %12.3f `e(ciUB`t'`k'_mivsr)' ")"
						di ""
					}

					if "`moret'" != "" {
						forval i=1/`conta' {
							di "ATE(`t`i'',`s`i'')" "{col 15}[" %-12.3f max(`e(LB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(UB`t`i''`s`i''_mivsr)' "]"
							di "`level'% CLR CI" "{col 15}(" %-12.3f max(`e(ciLB`t`i''`s`i''_mivsr)',0) "," %12.3f `e(ciUB`t`i''`s`i''_mivsr)' ")"
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

/*Create the code for a LaTeX table of results when the option is used*/

if "`latex'" != "" {
	file open textable using "`latex'.tex", write replace
	file write textable "\begin{table}[h!]" _n ///
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
		file write textable "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtr)',0)) "," %12.3f (`e(UB`t'`k'_mtr)') "] & [" %-12.3f (max((`e(LB`t'`k'_mtsr)',0)) "," %12.3f (`e(UB`t'`k'_mtsr)') "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (max(`e(LB`t'`k'_mivr)',0)) "," %12.3f (`e(UB`t'`k'_mivr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mivsr)',0)) "," %12.3f (`e(UB`t'`k'_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtsr)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (max(`e(ciLB`t'`k'_mivr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mivsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write textable "\midrule" _n ///
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
		file write textable "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (`e(LB`t'`k'_mtr)') "," %12.3f (min(`e(UB`t'`k'_mtr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mtsr)') "," %12.3f (min(`e(UB`t'`k'_mtsr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (`e(LB`t'`k'_mivr)') "," %12.3f (min(`e(UB`t'`k'_mivr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivsr)') "," %12.3f (min(`e(UB`t'`k'_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (`e(ciLB`t'`k'_mtr)') "," %12.3f (min(`e(ciUB`t'`k'_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mtsr)') "," %12.3f (min(`e(ciUB`t'`k'_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (`e(ciLB`t'`k'_mivr)') "," %12.3f (min(`e(ciUB`t'`k'_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivsr)') "," %12.3f (min(`e(ciUB`t'`k'_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write textable "\midrule" _n ///
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
		file write textable "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (`e(LB`t'`k'_mtr)') "," %12.3f (min(`e(UB`t'`k'_mtr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mtsr)') "," %12.3f (min(`e(UB`t'`k'_mtsr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (`e(LB`t'`k'_mivr)') "," %12.3f (min(`e(UB`t'`k'_mivr)',0)) "] & [" %-12.3f (`e(LB`t'`k'_mivsr)') "," %12.3f (min(`e(UB`t'`k'_mivsr)',0)) "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (`e(ciLB`t'`k'_mtr)') "," %12.3f (min(`e(ciUB`t'`k'_mtr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mtsr)') "," %12.3f (min(`e(ciUB`t'`k'_mtsr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (`e(ciLB`t'`k'_mivr)') "," %12.3f (min(`e(ciUB`t'`k'_mivr)',0)) ") & (" %-12.3f (`e(ciLB`t'`k'_mivsr)') "," %12.3f (min(`e(ciUB`t'`k'_mivsr)',0)) ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write textable "\midrule" _n ///
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
		file write textable "$\mathit{ATE(`t',`k')}$ & [" %-12.3f (`e(LB`t'`k'_wc)') "," %12.3f (`e(UB`t'`k'_wc)') "] & [" %-12.3f (`e(LB`t'`k'_mivo)') "," %12.3f (`e(UB`t'`k'_mivo)') "] & [" %-12.3f (`e(LB`t'`k'_mts)') "," %12.3f (`e(UB`t'`k'_mts)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtr)',0)) "," %12.3f (`e(UB`t'`k'_mtr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mtsr)',0)) "," %12.3f (`e(UB`t'`k'_mtsr)') "] & [" %-12.3f (`e(LB`t'`k'_mivs)') "," %12.3f (`e(UB`t'`k'_mivs)') "] & [" %12.3f (max(`e(LB`t'`k'_mivr)',0)) "," %12.3f (`e(UB`t'`k'_mivr)') "] & [" %-12.3f (max(`e(LB`t'`k'_mivsr)',0)) "," %12.3f (`e(UB`t'`k'_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t'`k'_wc)') "," %12.3f (`e(ciUB`t'`k'_wc)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivo)') "," %12.3f (`e(ciUB`t'`k'_mivo)') ") & (" %-12.3f (`e(ciLB`t'`k'_mts)') "," %12.3f (`e(ciUB`t'`k'_mts)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mtsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mtsr)') ") & (" %-12.3f (`e(ciLB`t'`k'_mivs)') "," %12.3f (`e(ciUB`t'`k'_mivs)') ") & (" %12.3f (max(`e(ciLB`t'`k'_mivr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivr)') ") & (" %-12.3f (max(`e(ciLB`t'`k'_mivsr)',0)) "," %12.3f (`e(ciUB`t'`k'_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	if "`moret'" != "" {
		forval i=1/`conta' {
			file write textable "\midrule" _n ///
			"\addlinespace" _n ///
			"$\mathit{ATE(`t`i'',`s`i'')}$ & [" %-12.3f (`e(LB`t`i''`s`i''_wc)') "," %12.3f (`e(UB`t`i''`s`i''_wc)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivo)') "," %12.3f (`e(UB`t`i''`s`i''_mivo)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mts)') "," %12.3f (`e(UB`t`i''`s`i''_mts)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mtsr)') "] & [" %-12.3f (`e(LB`t`i''`s`i''_mivs)') "," %12.3f (`e(UB`t`i''`s`i''_mivs)') "] & [" %12.3f (max(`e(LB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivr)') "] & [" %-12.3f (max(`e(LB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(UB`t`i''`s`i''_mivsr)') "]\\" _n ///
		"\addlinespace" _n ///
		" & (" %-12.3f (`e(ciLB`t`i''`s`i''_wc)') "," %12.3f (`e(ciUB`t`i''`s`i''_wc)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivo)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivo)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mts)') "," %12.3f (`e(ciUB`t`i''`s`i''_mts)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mtsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mtsr)') ") & (" %-12.3f (`e(ciLB`t`i''`s`i''_mivs)') "," %12.3f (`e(ciUB`t`i''`s`i''_mivs)') ") & (" %12.3f (max(`e(ciLB`t`i''`s`i''_mivr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivr)') ") & (" %-12.3f (max(`e(ciLB`t`i''`s`i''_mivsr)',0)) "," %12.3f (`e(ciUB`t`i''`s`i''_mivsr)') ")\\" _n ///
		"\addlinespace" _n
	}
	}
	}
	
	
	file write textable "\midrule" _n ///
	"\addlinespace" _n ///
	"Observations & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' & `e(Nob)' \\" _n ///
	"\addlinespace" _n ///
	"\bottomrule" _n ///
	"\multicolumn{8}{l}{\footnotesize  [$\cdot$]: half-median unbiased MIV estimates \ ; \ ($\cdot$): `level'\% CLR confidence intervals}" _n ///
	"\end{tabular}" _n ///
	"}" _n ///
	"\caption{\small Bounds on the ATE: `miv' as MIV, `bins' Bins, with \$K_0=`ymin'$}" _n ///
	"\end{table}"
	

file close textable

}
	
end
