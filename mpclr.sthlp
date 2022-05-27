{smcl}
{cmd:help mpclr}{right: ({placeholder, would be stata journal link})}
{hline}

{title:Title}

{p2colset 6 14 16 2}{...}
{p2col:{hi:mpclr} {c -}}Estimates for average treatment effects using Manski and Pepper (2000) bounds with bias correction{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:mpclr}
{depvar}
{it:treatvar}
{ifin}
{weight}{cmd:,} {opt ymin(#)} {opt ymax(#)}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opt ymin(#)}} assumed lower bound on the support of the outcome variable
{p_end}
{p2coldent:* {opt ymax(#)}} assumed upper bound on the support of the outcome variable
{p_end}
{synopt:{opt att}} requests estimates of bounds on the average treatment effect on the treated (ATT){p_end}
{synopt:{opth bounds(string)}} option to specify a subset of the available bounds (without an MIV) to estimate{p_end}
{synopt:{opt nmts}} specify non-positive monotone treatment selection (MTS) is to be used{p_end}
{synopt:{opt nmtr}} specify non-positive monotone treatment response (MTR) is to be used{p_end}
{synopt:{opt miv(varname)}} specify the monotone instrumental variable (MIV){p_end}
{synopt:{opth mivbounds(string)}} option to specify a subset of the available MIV bounds to estimate{p_end}
{synopt:{opt bins(#)}} specify the number of MIV bins to be used in the estimator{p_end}
{synopt:{opt nmiv}} specify that a non-positive MIV assumption is to be used{p_end}
{synopt:{opth moret(string)}} request additional treatment level comparisons; see below{p_end}
{synopt:{opth mores(string)}} request additional treatment level comparisons; see below{p_end}
{synopt:{opt level(#)}} set confidence level for the confidence intervals returned by the Chernozhukov, Lee, and Rosen (2013) procedure{p_end}
{synopt:{opt reps(#)}} specify number of bootstrap replications for the first step of the CLR procedure{p_end}
{synopt:{opt seed(#)}} set random number generator seed{p_end}
{synopt:{opt noisyclr}} display the intermediate CLR procedure output{p_end}
{synopt:{opth latex(string)}} specify {it:string}.tex to output a file containing the results in a LaTeX table{p_end}
{synopt:{opt survey}} declare the estimation sample as survey data{p_end}
{synopt:{opt npsu(#)}} set the number of PSUs per stratum for bootstrap samples when {opt survey} is used{p_end}
{synopt:{opt brr}} specify that the data includes balanced repeated replicate (BRR) weights for bootstrapping{p_end}
{synopt:{opt uncorrected}} request that the plug-in estimates of the bounds involving the MIV are stored in {cmd:e()} along with the unadjusted confidence intervals on the corresponding ATEs{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt ymin(#)} and {opt ymax(#)} are both required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mpclr} estimates bounds on the average treatment effect of different levels of {it:treatvar} on the outcome {depvar} under various combinations of the monotonicity assumptions of Manski and Pepper (2000).  

{pstd}
These are of the form {it:ATE(t,s)}, which is the effect of receiving treatment level {it:t} versus treatment level {it:s}; by default, {cmd:mpclr} produces all available estimates from setting {it:s = t-1}, but other comparisons may be added.

{pstd}
A total of eight sets of bounds may be estimated:

{p2colset 7 12 15 2}{...}
{p2col: (1)} no-assumption, worst-case selection bounds {p_end}
{p2col: (2)} bounds under a monotone instrumental variable assumption (MIV-only), if an MIV is specified{p_end}
{p2col: (3)} bounds under monotone treatment selection (MTS-only){p_end}
{p2col: (4)} bounds under monotone treatment response (MTR-only){p_end}
{p2col: (5)} bounds under monotone treatment selection and monotone treatment response (MTS+MTR){p_end}
{p2col: (6)} bounds under monotone treatment selection combined with a monotone instrumental variable (MIV+MTS), if an MIV is specified{p_end}
{p2col: (7)} bounds under monotone treatment response combined with a montone instrumental variable (MIV+MTR), if an MIV is specified{p_end}
{p2col: (8)} bounds under monotone treatment selection and monotone treatment response combined with a monotone instrumental variable (MIV+MTS+MTR), if an MIV is specified{p_end}

{pstd}
As noted in Manski and Pepper (2000), the MIV estimator suffers from bias in finite samples, such that the width of the estimated bounds is biased downwards.

{pstd}
To correct for this, {cmd:mpclr} implements the correction procedure of Chernozhukov, Lee, and Rosen (2013) to
 obtain half-median unbiased estimates of the upper and lower bounds under the MP-type assumptions. Thus,
 {cmd:mpclr} makes use of a secondary command {cmd:CLR} to perform this procedure.

{pstd}
The same procedure also yields confidence intervals for the true average treatment effect of interest at the
desired level.


{marker options}{...}
{title:Options}

{phang}
{opt ymin(#)} and {opt ymax(#)} set the respective lower and upper bounds for the outcome variable {depvar}, 
corresponding to K0 and K1 in Manski and Pepper (2000). {opt ymin()} and {opt ymax()} are required. Note that
observations of {it:depvar} outside of these values are excluded from the estimation. The user's dataset
is unchanged.

{phang}
{opt att} requests estimates of bounds on the average treatment effect on the treated (ATT). The default parameter is the population average treatment effect (ATE).

{phang}
{opt bounds(string)} allows the user to specify that only a subset of the non-MIV bounds are to be estimated;
supported inputs are {cmd:wc} (worst-case bounds), {cmd:mts} (MTS bounds), {cmd:mtr} (MTR bounds), and
{cmd:mtsr} (MTS+MTR bounds). The default is {cmd:bounds(wc mts mtr mtsr)}. If none of these are desired,
{cmd:bounds(none)} may be used.

{phang}
{opt nmts} requests that non-positive monotone treatment selection be assumed.  The default is the non-negative 
monotone treatment selection assumption, where those with weakly higher mean potential outcomes select into higher levels of treatment.

{phang}
{opt nmtr} requests that non-positive monotone treatment response be assumed.  The default is the non-negative 
monotone treatment response assumption, where potential outcomes are non-decreasing in the level of treatment.

{phang}
{opth miv(varname)} sets {varname} as the monotone instrumental variable.

{phang}
{opt mivbounds(string)} allows the user to specify that only a subset of the MIV bounds are to be estimated when
{cmd:miv()} is set. The supported inputs are {cmd:wcv} (MIV-only bounds), {cmd:mtsv} (MIV+MTS bounds), {cmd:mtrv} 
(MIV+MTR bounds), and {cmd:mtsrv} (MIV+MTS+MTR bounds). The default is {cmd:mivbounds(wcv mtsv mtrv mtsrv)}.

{phang}
{opt bins(#)} creates # bins of the MIV for the estimator, where these bins are defined by equally-spaced 
cutoffs at set quantiles of the empirical distribution of {varname}. The default is {cmd: bins(5)} and would 
demarcate bins using quantile intervals (0, 20], (20, 40], (40, 60], (60, 80], and (80, 100).
The {cmd:mpclr} command supports up to 8 bins using {cmd:bins(8)}.

{phang}
{opt nmiv} requests that {varname} be used as a non-positive monotone instrumental variable. The default is a 
non-negative MIV assumption, where mean potential outcomes are non-decreasing in the MIV.

{phang}
{opt moret(string)} and {opt mores(string)} requests the estimation of bounds on treatment effects other than {it:ATE(t, t-1)}.
For example, if there are 5 levels of treatment, then the default behavior of {cmd:mpclr} is to estimate bounds
on {it:ATE(2,1)}, {it:ATE(3,2)}, {it:ATE(4,3)}, and {it:ATE(5,4)}. But more general bounds on an {it:ATE(t,s)}
can be requested. In the example, if {it:ATE(4,2)} and {it:ATE(5,2)} are also of interest, then the user
should add the options {cmd:moret(4 5)} and {cmd:mores(2 2)}. Note that the order of input for additional
comparisons is essential, and that each option must have an equal number of elements ({it:i.e.}, the {cmd:2} in 
the example must be typed twice).

{phang}
{opt level(#)} requests that #% confidence intervals be returned by the Chernozhukov, Lee, and Rosen (2013) 
procedure. # should be an integer between 1 and 99. The default is {cmd:level(95)} for 95% confidence intervals 
around the average treatment effect of interest.

{phang}
{opt reps(#)} specifies that # bootstrap replications are to be performed in the first step of the CLR procedure. 
The default is set to {cmd:reps(25)} for quickness.

{phang}
{opt seed(#)} sets # as the random number generator seed for the bootstrap replications. The default is {cmd:seed(22176)}.

{phang}
{opt noisyclr} requests that the intermediate calculations from the CLR procedure are displayed. The default is 
that this portion is run quietly. A message is shown to note that the procedure has begun to run in either case.

{phang}
{opt latex(string)} saves a file {it:string}.tex in the active directory. This populates the estimation results
 into the code for a LaTeX table. Estimates of the lower and upper bounds are set to round to three decimal 
 places. The labels for treatment effects and the title given to the table are placeholders and may be modified 
 in the resulting .tex file according to the user's preferences. Note that there will be an error when trying 
 to compile a file with this table if the specified MIV has an underscore in the variable name. To remedy this, 
 remove the underscore in the table caption (or add the requisite backslash for the TeX syntax).
 
{phang}
{opt survey} declares that the data for estimation is survey data. Accordingly, the data must be {cmd:svyset} 
prior to using this option. When {cmd:survey} is invoked, the program executes the {helpb bsweights} command 
(Kolenikov 2010) to compute re-scaled weights for the number of bootstrap replications specified by {cmd:reps()}.

{phang}
{opt npsu(#)} can be used along with {cmd:survey} when specified. This tells the {cmd:bsweights} command how many 
PSUs to take per stratum for each bootstrap sample. The default is set to {cmd:npsu(-1)}, which says that, if 
there are n_h PSUs, the bootstrap samples will take m_h = n_h -1. This follows the discussion in Kolenikov (2010).

{phang}
{opt brr} is an option that should be used if the data includes balanced repeated replicate (BRR) weights for bootstrapping. When used, the data should be {cmd:svyset} according to the data's user guide using the {opth brrweight(varlist)} option.

{phang}
{opt uncorrected} requests that the plug-in estimates of the bounds (that is, bounds without the bias correction)
 involving the MIV are stored in {cmd:e()}. Note that the resulting estimates will therefore have the potential to be severely biased, in the sense that the bounds will be narrower than the true identified set. Unadjusted
 confidence intervals for the ATE will also be stored in {cmd:e()}.

{marker remarks}{...}
{title:Remarks}

{pstd}
For further details on the estimators used in {cmd:mpclr} and a discussion of the underlying monotonicity 
assumptions, see Manski and Pepper (2000) and Amin, Flores, Flores-Lagunes, and Germinario (2021). For details 
on the bias correction and inference methods employed here, see Chernozhukov, Lee, and Rosen (2013), Flores and
 Flores-Lagunes (2013), and Flores and Chen (2018).
 
{pstd}
We are grateful for the example of the {helpb tebounds} program code by McCarthy, Millimet, and Roy (2015), which
 was very helpful for implementing the {cmd:survey} option.
 
{pstd}
The {cmd:mpclr} command requires that {cmd:bs4rw}, {cmd:bsweights}, {cmd:egenmore}, {cmd:matselrc}, and
 {cmd:moremata} also be installed.


{marker references}{...}
{title:References}

{phang}
Chernozhukov, V., S. Lee, and A. Rosen. 2013. Intersection bounds: Estimation
and inference. {it:Econometrica} 81: 667-737.

{phang}
de Haan, M. 2011. The effect of parents' schooling on child's schooling: A nonparametric
bounds analysis. {it:Journal of Labor Economics} 29: 859-892.

{phang}
Flores, C., and A. Flores-Lagunes. 2013. Partial identification of local average
treatment effects with an invalid instrument.
{it: Journal of Business & Economic Statistics} 31: 534-545.

{phang}
Flores, C., and X. Chen. 2018. Average treatment effect bounds with an instrumental
variable: Theory and practice. Springer.

{phang}
Kolenikov, S. 2010. {browse "https://www.stata-journal.com/article.html?article=st0187":Resampling variance estimation for complex survey data.}
{it:Stata Journal} 10: 165-199.

{phang}
Manski, C.F. 1997. Monotone treatment response. {it:Econometrica} 65: 1311-1334.

{phang}
Manski, C. F., and J. V. Pepper.  2000.  Monotone instrumental
variables: With an application to the returns to schooling.
{it:Econometrica} 68: 997-1010.

{phang}
McCarthy, I., D. Millimet, and M. Roy. 2015.
{browse "http://www.stata-journal.com/article.html?article=st0386":Bounding treatment effects: A command for the partial identification of the average treatment effect with endogenous and misreported treatment assignment.}
{it:Stata Journal} 15: 411-436.


{marker authors}{...}
{title:Authors}

{pstd}Giuseppe Germinario{p_end}
{pstd}Syracuse University{p_end}
{pstd}Syracuse, NY{p_end}
{pstd}grgermin@syr.edu{p_end}

{pstd}Carlos Flores{p_end}
{pstd}California Polytechnic State University at San Luis Obispo{p_end}
{pstd}San Luis Obispo, CA{p_end}
{pstd}cflore32@calpoly.edu{p_end}

{pstd}Alfonso Flores-Lagunes{p_end}
{pstd}Syracuse University{p_end}
{pstd}Syracuse, NY{p_end}
{pstd}and IZA{p_end}
{pstd}Bonn, Germany{p_end}
{pstd}afloresl@maxwell.syr.edu{p_end}

{marker also_see}{...}
{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume xx, number x: {browse "http://www.stata-journal.com/":st0000}{p_end}

{p 7 14 2}Help: {helpb bsweights}, {helpb bs4rw}, {helpb egenmore},  {helpb matselrc}, {helpb moremata} (if installed){p_end}
