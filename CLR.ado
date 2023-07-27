* This Stata function is a tentative translation of Prof. Carlos Flores's Conf_Int_CLR, 12/26/11
* to Stata. 

* Xintong 1/26/2015

* All "my notes" indicates the notes "Implementing_CLR_2011_122111" by Prof. Carlos Flores

* Conf_Int_CLR. Based on the paper by Chernozhukov, Lee and Rosen (2013), 
* this function provides: (1) half-median unbiased estimators of lower and 
* upper bounds that may involve "min/max" operators; and, (2) 1-alfa
* percent confidence intervals for the parameter of interest (like those 
* in Imbens and Manski (2004), but allowing for min/max operators).

* SETUP. See notes "Implementing_CLR_2011_122111" pages 7 and 8. 
* Let theta_star be the parameter of interest, which lies between 
* (is bounded by) [Lstar, Ustar], where Lstar=max(L1,L2,...,LkL) 
* and Ustar=min(U1,U2,...,UkU). As mentioned above, we may have kL 
* and/or kU equal to one.

* INPUTS: (bL, VbL, bU, VbU, n, R, alfa). bL is the kL by 1 vector of 
* estimated LOWER bounds, and VbL is its kL by kL variance-covariance 
* matrix. bU is the kU by 1 vector of estimated UPPER bounds, and VbU is 
* its kU by kU var-cov matrix. n is the sample size used to generate those 
* estimates. R is the number of repetitions to be used to compute 
* half-unbiased estimators and confidence intervals (Recommended: at least 
* 100,000. See "Implementing_CLR_2011_122111" page 15). 
* alfa is the significance level, so that 
* 1-alfa is the confidence level, e.g., for a 95 percent CI use alfa=0.05.

* OUTPUT: [Lstar, Ustar, CIL, CIU]. Lstar and Ustar are the half-median 
* unbiased estimators of the lower and upper bounds, respectively. CIL and 
* CIU are the lower and upper ends of the (1-alfa) percent confidence 
* interval, respectively.

qui program define CLR, eclass
 version 10.0
 marksample touse
 set more off

 matrix bL=`1'
 matrix VbL=`2'
 matrix bU=`3'
 matrix VbU=`4'
 sca n=`5'
 sca R=`6'
 sca alfa=`7' 

 mata: clr("'bL'","'VbL'","'bU'","'VbU'","'n'","'R'","'alfa'")
 mat Lstar=e(Lstar)
 mat Ustar=e(Ustar)
 mat CIL=e(CIL)
 mat CIU=e(CIU)
 
 ereturn mat Lstar Lstar
 ereturn mat Ustar Ustar
 ereturn mat CIL CIL
 ereturn mat CIU CIU
 
end 

quietly mata:

function clr(matrix bL, matrix VbL, matrix bU, matrix VbU, scalar n, scalar R, scalar alfa)
 {
 bL=st_matrix("bL")
 VbL=st_matrix("VbL")
 bU=st_matrix("bU")
 VbU=st_matrix("VbU")
 n=st_numscalar("n")
 R=st_numscalar("R")
 alfa=st_numscalar("alfa")
 
 //****************** ***** I. GET HALF-MEDIAN UNBIASED ESTIMATE OF LOWER BOUND, Lstar *****
 ""
 "**************************************************************"
 ""
 "----- Consider LOWER BOUND: Lstar=max(L1,L2,...,LkL) -----"
 "--- I. Get half-median unbiased estimator of Lstar ---"
 ""
 
 //Note: Steps 1 to 4(a) are the SAME for lower and upper bound.
 
 //Step 1: Get random draws from normals and values of constants. 
 // Set value of constant, gamma tilda in CLR (13)
 C=1-(0.1/log(n)) 
 // Set the number of possible LOWER bounds, #rows in bL.
 kL=nonmissing(bL) 
 rseed(1989)
 //It generates R draws from a kL-variate normal 
 ZL=rnormal(kL, R, 0, 1) 
 //distribution mean zero and kL by kL identity matrix. KEY: For us, each
 //COLUMN is a draw.
 

 //STEP 2: Get a consistent estimate of the asymptotic distribution of 
 //sqrt(n)[bL-bLo].
 OmegaL=VbL*n
 
 
 //STEP 3: Here we compute some objects to be used below. 
 //Note: Here is a little different from the Matlab version. As Mata does not have 
 //the sqrtm function in Matlab that will calculates matrix A's sqaure root X, such 
 //that X*X=A, so we need to calculated it "manually" in Stata. 
  //The unique square root of OmegaL
 if (nonmissing(bL)==1) SqrtOmegaL=(OmegaL)^(1/2)
 else {
 eigensystem(OmegaL,X,L)
 sqrtLvec=sqrt(L)
 sqrtL=diag(sqrtLvec)
 SqrtOmegaL=Re(X*sqrtL*luinv(X))
 } 

 
 //sL is a kUx1 vector containing all the s(v) for v=1,...,kL
 sL=J(kL, 1, 0)
 for (i=1; i<=kL; i++){
     sL[i]=sqrt(SqrtOmegaL[i,.]*SqrtOmegaL[i,.]')/sqrt(n)
 }
 "Value of the vector of lower bound estimates, bL=[L1 L2 ... LkL]"
 bL
 "Value of sL (std. error of each lower bound estimate)"
 sL
 
 
 //STEP4: Get V hat. This set includes those estimators of the LOWER bound 
 //that can affect the asymptotic distribution of the estimator of the 
 //MAXIMUM, Lstar. (e.g., see "my notes" page 13).
 
 for (k=1; k<=kL; k++){
 if (sL[k,1]==0){
 sL[k,1]=epsilon(1)
 }
 }
 //Note: The above step is for the case when bL[i]=0, therefore sL[i]==0.
 //      While Matlab will assign infinite large to the result of ZstarL[v,.]/0,
 //      Stata will assign empty value "." in such a case, 
 //      which will lead to error in step 4(a)(ii).
 //      So I replace zero sL[i] with the value of approximately 2.22045e–16. This 
 //      is also applied to sU of the UPPER bound. 
 
 
 //4(a): Get Critical Value. 
 //(i): Get matrix ZstarL.  
 ZstarL=SqrtOmegaL*ZL
 for (v=1; v<=kL; v++){
 ZstarL[v,.]=ZstarL[v,.]/(sqrt(n)*sL[v])
 }
 //ZL can be a huge matrix and consume a lot of space.
 //In the original Matlab code, ZL is "cleared" from the memory.
 //There is no corresponding commend in mata, so I set ZL to be an empty matrix here,
 //which would serve the same purpose of saving space. 
 ZL=.

 //(ii): Get maxZstarL, which is a 1 by R row vector. 
 if (kL==1) maxZstarL=ZstarL
 else maxZstarL=colmax(ZstarL)
 
 ///(iii): Get the c-th quantile of maxZstarL. This is the critical value in 
 //Step 4a. 
 C=1-(0.1/log(n))
 cv1L=mm_quantile(maxZstarL',1,C,1)
 //Note: "moremata" package is needed to implement mm_quantile 
 "Critial value used to determine V hat"
 cv1L
 maxZstarL=.
 
 //4(b): ***THIS STEP IS LOWER-BOUND SPECIFIC***. Here we get the set V hat.
  
 //(i): Get MAXIMUM precision-adjusted estimator.
 maxv1=max(bL-cv1L*sL)
 //(ii): Create a kLx1 vector indicating whether the estimator belongs to 
 // Vhat (=1) or not (=0). 
 indVhatL=(bL:>=-2*cv1L*sL :+ maxv1)
 "Number of elements in Vhat (i.e., number of elements affecting the asymptotic distribution of Lstar):"
 sum(indVhatL) 
 "Elements belonging to Vhat (if =1)"
 indVhatL 
 //Vector can be long if kL is large (usually not our case).
 
 
 //STEP5: Get estimator of LOWER bound.
 
 //5(a): Get critical value.
 //(i):  We create the matrix Zstar_VhatL. It contains the rows from ZstarL
 //      corresponding to those estimators/bounds in Vhat.
 //Counter used to fill in the appropriate row of Zstar_VhatL.
 //ZstarL can be a huge matrix and consume a lot of space.

 Zstar_VhatL=J(sum(indVhatL), R, 0)
 i=1 
 for (v=1; v<=kL ; v++){
  if (indVhatL[v,1]==1){
     Zstar_VhatL[i,.]=ZstarL[v,.]
	 i=i+1
	 }
    }
 ZstarL=. 
 
 //(ii): Get a row vector containing the maximum element for each column of
 //Zstar_VhatL (i.e., the max value over "v IN Vhat" for each replication).
 //We use an "if" because if Zstar_VhatL is a vector, then the max function
 //returns a scalar. See notes "Implementing_CLR_2011_122111" page 11. 
 //maxZstar_VhatL is a 1 by R row vector. 
 //See also note (3) at the beginning of the m-file Preliminary_LB_122311.m.
 //KEY: This vector is also used below to get critical values for 
 //confidence intervals.
 if (sum(indVhatL)==1) maxZstar_VhatL=Zstar_VhatL
 else maxZstar_VhatL=colmax(Zstar_VhatL)
 Zstar_VhatL=. 
 //(iii): Get the p-th quantile of maxZstar_VhatL. This is the critical value 
 //in step 5a.
 //Use p=0.5 to get half-median unbiased 
 cv2L=mm_quantile(maxZstar_VhatL',1,0.5,1) 
 //estimator of lower bound.
 "Critical value used for half-median unbiased est. of lower bound"
 cv2L
 
 // 5(b): ***THIS STEP IS LOWER-BOUND SPECIFIC***. Get half-median unbiased 
 //estimator of the LOWER bound, Lstar. See notes "Implementing_CLR_2011_122111" page 13. 
 Lstar=max(bL-cv2L*sL)
 "***** Half-Median Unbiased Estimator of LOWER BOUND:"
 Lstar

 
 
 //***********II. GET HALF-MEDIAN UNBIASED ESTIMATE OF UPPER BOUND, Ustar***************
 
 "******************************************************************"
 ""
 "----- Consider UPPER BOUND: Ustar=min(U1,U2,...,UkU) -----" 
 "--- II. Get half-median unbiased estimator of Ustar ---" 
 ""
 //Note: Steps 1 to 4(a) are the SAME for lower and upper bound.
 
 //Step 1: Get random draws from normals and values of constants. 
 // Note that "c" is the same as for lower bound above.
 kU=nonmissing(bU)
 rseed(1989)
 ZU=rnormal(kU, R, 0, 1)
 
 //STEP 2: Get a consistent estimate of the asymptotic distribution of 
 //sqrt(n)[bU-bUo].
 OmegaU=VbU*n
 
 //STEP 3: Here we compute some objects to be used below. 
 //The unique square root of Omega
 if (nonmissing(bU)==1) SqrtOmegaU=(OmegaU)^(1/2)
 else {
 eigensystem(OmegaU,X,U)
 sqrtUvec=sqrt(U)
 sqrtU=diag(sqrtUvec)
 SqrtOmegaU=Re(X*sqrtU*luinv(X))
 }
 

 //sU is a kUx1 vector containing all the s(v) for v=1,...,kU
 sU=J(kU, 1, 0)
 for (i=1; i<=kU; i++){
 sU[i]=sqrt(SqrtOmegaU[i,.]*SqrtOmegaU[i,.]')/sqrt(n)
 }
 "Value of the vector of upper bound estimates, bU=[U1 U2 ... UkU]"
 bU
 "Value of sU (std. error of each upper bound estimate)"
 sU
 
 
 //STEP4: Get V hat. This set includes those estimators of the UPPER bound 
 //that can affect the asymptotic distribution of the estimator of the 
 //MAXIMUM, Ustar. 
 
 for (k=1; k<=kU; k++){
 if (sU[k,1]==0){
 sU[k,1]=epsilon(1)
 }
 }
 
 //4(a): Get Critical value. 
 //(i): Get matrix ZstarU. 
 ZstarU=SqrtOmegaU*ZU
 for (v=1; v<=kU; v++){
     ZstarU[v,.]=ZstarU[v,.]/(sqrt(n)*sU[v])
 }
 ZU=. 
 //ZU can be a huge matrix and consume a lot of space.
 
 //(ii): Get maxZstarU, which is a 1 by R row vector. 
 if (kU==1) maxZstarU=ZstarU
 else maxZstarU=colmax(ZstarU)
 
 ///(iii): Get the c-th quantile of maxZstarU. This is the critical value in 
 //Step 4a. 
 C=1-(0.1/log(n))
 cv1U=mm_quantile(maxZstarU',1,C,1)
 "Critial value used to determine V hat"
 cv1U 
 maxZstarU=.
 
 //4(b): ***THIS STEP IS UPPER-BOUND SPECIFIC***. Here we get the set V hat.
 //      See my notes page 11.
 //(i): Get MINIMUM precision-adjusted estimator.
 minv1=min(bU+cv1U*sU)
 //(ii): Create a kUx1 vector indicating whether the estimator belongs to 
 //      Vhat (=1) or not (=0). 
 indVhatU=(bU:<=2*cv1U*sU :+ minv1)
 "Number of elements in Vhat (i.e., number of elements affecting the asymptotic distribution of Ustar):"
 sum(indVhatU) 
 "Elements belonging to Vhat (if =1)"
 indVhatU //Vector can be long if kU is large (usually not our case).
 
 
 //STEP5: Get estimator of UPPER bound.
 
 //5(a): Get critical value.
 //(i): We create the matrix Zstar_VhatU. It contains the rows from ZstarU 
 //corresponding to those estimators/bounds in Vhat.
 
 Zstar_VhatU=J(sum(indVhatU), R, 0)
 //Counter used to fill in the appropriate row of Zstar_VhatU.
 i=1 
 for (v=1; v<=kU ; v++){
  if (indVhatU[v,1]==1){
     Zstar_VhatU[i,.]=ZstarU[v,.]
	 i=i+1
	 }
    }
 ZstarU=.
 
 //(ii): Get a row vector containing the maximum element for each column of
 //Zstar_VhatU (i.e., the max value over "v IN Vhat" for each replication).
 //We use an "if" because if Zstar_VhatU is a vector, then the max function
 //returns a scalar. See notes "Implementing_CLR_2011_122111" page 11. 
 //maxZstar_VhatU is a 1 by R row vector. 
 
 //See also note (4) at the beginning of the m-file 
 //Preliminary_122211.m.
 //KEY: This vector is also used below to get critical values for 
 //confidence intervals.
 if (sum(indVhatU)==1) maxZstar_VhatU=Zstar_VhatU
 else maxZstar_VhatU=colmax(Zstar_VhatU)
 Zstar_VhatU=. 
 
 //(iii): Get the p-th quantile of maxZstar_VhatU. This is the critical 
 //value in step 5a.
 cv2U=mm_quantile(maxZstar_VhatU',1,0.5,1)
"Critical value used for half-median unbiased est. of upper bound"
 cv2U 
 
 // 5(b): ***THIS STEP IS UPPER-BOUND SPECIFIC***. Get half-median unbiased 
 //estimator of the UPPER bound, Ustar. See notes "Implementing_CLR_2011_122111" page 4. 
 Ustar=min(bU+cv2U*sU)
 "***** Half-Median Unbiased Estimator of UPPER BOUND:"
 Ustar

 
 //***** III. GET (1-alfa) PERCENT CONFIDENCE INTERVAL FOR PARAMETER *****
 "**************************************************************"
 ""
 "--- III. Get (1-alfa) percent Confidence Interval for Parameter---"
 "--- Confidence Interval has the form: [CIL, CIU] ---"
 ""
 
 //Note: See my notes pages 17, 5 and CLR (09) pages 25 and 26.
 
 //STEP 1: Get the values of the variables used to determine the critical 
 //value to be used to construct the confidence interval.
 vector1= (0 \ Ustar-Lstar)
 Delta_hat_plus=max(vector1)
 //Get sigma, see "Implementing_CLR_2011_122111" pages 17 & 18.
 U75=min(bU+(mm_quantile(maxZstar_VhatU',1,0.75,1))*sU)
 U25=min(bU+(mm_quantile(maxZstar_VhatU',1,0.25,1))*sU)
 L75=max(bL-(mm_quantile(maxZstar_VhatL',1,0.75,1))*sL)
 L25=max(bL-(mm_quantile(maxZstar_VhatL',1,0.25,1))*sL)
 
 vector2=(U75-U25 \ L25-L75)
 sigma=max(vector2)
 tao=1/(sigma*log(n))
 phat=1-normal(tao*Delta_hat_plus)*alfa 
 //NOTE: This is the only thing 
 //that changes in the program when we want to change the confidence level.
 
 //For reference, we compute the phat we would use if instead we used 
 //tao=log(n), as discussed in CLR(09) p. 26, or as used in Nevo & Rosen 
 //(2010) page 23
 phat_alternative=1-normal(log(n)*Delta_hat_plus)*alfa
 
 "Quantities used to get tao & phat, see CRL(09) p. 26; notes Implementing_CLR_2011_122111 p.18"
 ""
 "Delta_hat_plus = max(0, Ustar-Lstar)"
 Delta_hat_plus
 "U75"
 U75
 "U25"
 U25
 "L75"
 L75
 "L25"
 L25
 "U75-U25"
 U75-U25
 "L25-L75"
 L25-L75
 "sigma"
 sigma
 "Value of tao when taking into account the variance of Delta hat."
 "This is the preferred choice of tao. See references above"
 tao
 "Specified Confidence Level (for reference, 1-significance level):"
 1-alfa
 "Value of p_hat. This is the one we use"
 phat
 
 "For reference, this is the phat we would use if instead we used"
 " tao=log(n). E.g., see CLR(09) p. 26, or Nevo & Rosen (2010) p. 23"
 phat_alternative
 
 
 //STEP 2: Get the (1-alfa) percent confidence interval.
 //Lower end of confidence interval.
 //Critical value usedxs
 cvCIL=mm_quantile(maxZstar_VhatL',1,phat,1) 
 ""
 ""
 "Critical value used for LOWER END of (1-alfa)% Confidence Interval"
 cvCIL
 "***** LOWER END of (1-alfa)% CONFIDENCE INTERVAL, CIL:"
 CIL=max(bL-cvCIL*sL)
 CIL
 
 //Upper end of confidence interval
 //Critical value used
 cvCIU=mm_quantile(maxZstar_VhatU',1, phat,1)
 "Critical value used for UPPER END of (1-alfa)% Confidence Interval"
 cvCIU
 "***** UPPER END of (1-alfa)% CONFIDENCE INTERVAL, CIU:"
 CIU=min(bU+cvCIU*sU)
 CIU
 
 //Finally, FOR REFERENCE, I report the confidence interval resulting from 
 //using tao=log(n). See comments above.
 ""
 ""
 "For reference, we provide the CI with tao=log(n). See notes above."
 CI_L_alternative=max(bL-(mm_quantile(maxZstar_VhatL',1,phat_alternative,1))*sL)
 "Lower end of (1-alfa)% CI with NON-PREFERRED tao=log(n)."
 CI_L_alternative

 CI_U_alternative=min(bU+(mm_quantile(maxZstar_VhatU',1,phat_alternative,1))*sU)
 "Upper end of (1-alfa)% CI with NON-PREFERRED tao=log(n)."
 CI_U_alternative

 displayas("txt")
 printf("finished running")
 
 st_matrix("eLstar", Lstar)
 st_matrix("eUstar", Ustar)
 st_matrix("eCIL", CIL)
 st_matrix("eCIU", CIU)
 }

end 

