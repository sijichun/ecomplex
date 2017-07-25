*! version 14.0 Author : Jichun Si
* www.sijichun.pro
* The Economic Complexity index;
* Reference: The impact of services on Economic Complexity: Service sophistication as Route for Economic Growth
*    PLOS One, Stojkoski, Utkovski and Kocarev (2016).

cap program drop ecomplex_mr
program define ecomplex_mr
	version 14
	syntax varlist(max=1), country(varname) product(varname) time(varname) gen_country(name) gen_product(name)
	** generate second lagest eigen value for every year
	tempname time_list iter
	scalar `iter'=1
	quietly: levelsof `time', local(`time_list')
	preserve
	tempfile temp_c temp_p using_c using_p
	foreach t of local `time_list'{
		** firstly, country
		restore, preserve
		quietly: keep if `time'==`t'
		rename `varlist' __M__
		quietly: keep `country' `product' __M__ `time'
		quietly: reshape wide __M__, i(`country') j(`product')
		mata: ecomplex_eigen()
		quietly: svmat eigvec, names(`gen_country')
		keep `country' `time' `gen_country'*
		if `iter'==1{
			quietly: save `using_c', replace
		}
		else{
			quietly: save `temp_c', replace
			quietly: use `using_c', clear
			quietly: append using `temp_c'
			quietly: save`using_c', replace
		}
		** now, products
		restore, preserve
		quietly: keep if `time'==`t'
		rename `varlist' __M__
		quietly: keep `country' `product' __M__ `time'
		quietly: reshape wide __M__, i(`product') j(`country')
		mata: ecomplex_eigen()
		quietly: svmat eigvec, names(`gen_product')
		keep `product' `time' `gen_product'*
		if `iter'==1{
			quietly: save `using_p', replace
		}
		else{
			quietly: save `temp_p', replace
			quietly: use `using_p', clear
			quietly: append using `temp_p'
			quietly: save`using_p', replace
		}
		** iter+1
		scalar `iter'=`iter'+1
	}
	restore
	quietly: merge m:1 `country' `time' using `using_c'
	quietly: drop if _merge==2
	quietly: drop _merge
	quietly: merge m:1 `product' `time' using `using_p'
	quietly: drop if _merge==2
	quietly: drop _merge
end

version 14.0
cap mata : mata drop ecomplex_eigen()
mata:
mata set matastrict on
	void ecomplex_eigen(){
		real matrix M, fm, Mstar, MM
		complex matrix X, L
		real colvector Mrowsum
		real rowvector Mcolsum
		complex colvector v
    M=st_data(., "__M__*")
		Mrowsum=rowsum(M)
		Mcolsum=colsum(M)
		fm=Mrowsum*Mcolsum
		Mstar=M:/fm
    MM=M*Mstar'
		// eigenvector
		eigensystemselecti(MM, (1,2), X=., L=.)
		v=Re(X[.,2])
		st_matrix("eigvec",v)
	}
end

cap program drop ecomplex
program define ecomplex
	version 14
	syntax varlist(max=1), country(varname) product(varname) time(varname) [, mreig mr(integer 0) fcm mfcm gen(name) rca cv(real 1) ignore(integer 0) maxiter(integer 500) tol(real 1e-6) display]
	if `mr'<=0 & "`fcm'"=="" & "`mfcm'"=="" & "`mreig'"==""{
		di as err "Please give one or more methods, i.e. mr(d) or(and) mreig or(and) fcm or(and) mfcm."
		exit 198
	}
	if `mr'<=0{
		di as err "Positive integer 'd' should be given in the option mr(d)."
		exit 198
	}
	** if RCA index is given, rca should be set, then go on;
	** otherwise, the value of export should be given, and we will compute the rca first.
	tempname _rca
	if "`rca'"=="rca" {
		gen `_rca'=`varlist'
	}
	else{
		rcaindex `varlist', country(`country') product(`product') time(`time') gen(`_rca')
	}
	tempvar check_duplicates
	** check duplicates
	quietly: duplicates tag `country' `time' `product', gen(`check_duplicates')
	quietly: su `check_duplicates'
	if r(max)>0{
		di as err "Please check duplicates for country, time and products."
		exit 198
	}
	** balance the country*year*products
	tempfile temp_master temp_country temp_product temp_using
	** Cartesian product of country, product and time
	preserve
	quietly: keep `country'
	quietly: duplicates drop `country', force
	quietly: save `temp_country'
	restore, preserve
	quietly: keep `product'
	quietly: duplicates drop `product', force
	quietly: save `temp_product'
	restore, preserve
	quietly: keep `time'
	quietly: duplicates drop `time', force
	quietly: cross using `temp_country'
	quietly: cross using `temp_product'
	quietly: save `temp_using'
	restore, preserve
	quietly: merge 1:1 `country' `product' `time' using `temp_using'
	quietly: drop _merge
	quietly: replace `varlist'=0 if `varlist'==.
	** compute M
	tempvar M
	quietly: gen `M'=`_rca'>=`cv'
	** now, compute d
	tempvar d u
	quietly: egen `d'=total(`M'), by(`country' `time')
	** ignore countries with few comparative advantage goods
	quietly: drop if `d'<=`ignore'
	quietly: drop `d'
	** now, compute u
	quietly: egen `u'=total(`M'), by(`product' `time')
	** ignore products with no country
	quietly: drop if `u'==0
	quietly: drop `u'
	** recompute d and u since we droped some countries and products
	quietly: egen `d'=total(`M'), by(`country' `time')
	quietly: egen `u'=total(`M'), by(`product' `time')
	** num of country and num of products
	tempvar num_country num_product
	quietly: egen `num_country'=count(`time'), by(`time' `product')
	quietly: egen `num_product'=count(`time'), by(`time' `country')
	tempname keepvar_mreig keepvar_fcm keepvar_mfcm keepvar_mr
	********* Method of Reflections ********
	local `keepvar_mreig' ""
	tempvar mr_temp_country mr_temp_product
	if "`mreig'"=="mreig"{
		if "`display'"=="display"{
			disp "Compute Complexity Index with the Method of Reflections (MR), 2nd largest eigenvector is calculated."
		}
		quietly: ecomplex_mr `M', country(`country') product(`product') time(`time') gen_country(`mr_temp_country') gen_product(`mr_temp_product')
		tempvar mean_p std_p mean_c std_c
		quietly: egen `mean_c'=mean(`mr_temp_country'), by(`time' `product')
		quietly: egen `std_c' =  sd(`mr_temp_country'), by(`time' `product')
		quietly: egen `mean_p'=mean(`mr_temp_product'), by(`time' `country')
		quietly: egen `std_p' =  sd(`mr_temp_product'), by(`time' `country')
		quietly: gen `gen'_MREIG_country=(`mr_temp_country'-`mean_c')/`std_c'
		quietly: gen `gen'_MREIG_product=(`mr_temp_product'-`mean_p')/`std_p'
		local `keepvar_mreig' "`gen'_MREIG_country `gen'_MREIG_product"
	}
	*** init for iterations
	tempvar c_n_1 p_n_1 c_n p_n temp_sum_c temp_sum_p abs_error sum_c_n sum_p_n p_n_temp
	tempname distance iter
	quietly: gen `abs_error'=.
	quietly: gen `c_n_1'=1
	quietly: gen `p_n_1'=1
	********* Method of Reflections - one-step iteration ********
	local `keepvar_mr' ""
	if `mr'>0{
		if "`display'"=="display"{
			disp "Compute Complexity Index with the Method of Reflections (MR), with iteration times=" `mr'
		}
		** init **
		quietly: replace `c_n_1'=`d'
		quietly: replace `p_n_1'=`u'
		scalar `iter'=0
		** begin iteration
		while  `iter'<`mr'{
			quietly: gen `temp_sum_c'=`M'*`p_n_1'/`d'
			quietly: gen `temp_sum_p'=`M'*`c_n_1'/`u'
			quietly: egen `c_n'=total(`temp_sum_c'), by(`country' `time')
			quietly: egen `p_n'=total(`temp_sum_p'), by(`product' `time')
			scalar `iter'=`iter'+1
			quietly: drop `temp_sum_c' `temp_sum_p' `p_n_1' `c_n_1'
			rename `c_n' `c_n_1'
			rename `p_n' `p_n_1'
		}
		quietly: gen `gen'_MR`mr'_country=`c_n_1'
		quietly: gen `gen'_MR`mr'_product=`p_n_1'
		local `keepvar_mr' "`gen'_MR`mr'_country `gen'_MR`mr'_product"
	}
	********* Fitness-Complexity Method ********
	local `keepvar_fcm' ""
	if "`fcm'"=="fcm"{
		if "`display'"=="display"{
			disp "Compute Complexity Index with the Fitness-Complexity Method (FCM)"
		}
		** init **
		quietly: replace `abs_error'=.
		quietly: replace `c_n_1'=1
		quietly: replace `p_n_1'=1
		scalar `distance'=100000000.0
		scalar `iter'=0
		** begin iteration
		while `distance'>`tol' & `iter'<=`maxiter'{
			quietly: gen `temp_sum_c'=`M'*`p_n_1'
			quietly: gen `temp_sum_p'=`M'/`c_n_1'
			quietly: egen `c_n'     =total(`temp_sum_c'), by(`country' `time')
			quietly: egen `p_n_temp'=total(`temp_sum_p'), by(`product' `time')
			quietly: gen `p_n'=1/`p_n_temp'
			** normalize
			quietly: egen `sum_c_n'=total(`c_n'), by(`time')
			quietly: replace `c_n'=`c_n'/`sum_c_n'*`num_country'*`num_product'
			quietly: egen `sum_p_n'=total(`p_n'), by(`time')
			quietly: replace `p_n'=`p_n'/`sum_p_n'*`num_country'*`num_product'
			** distance between iterates
			quietly: replace `abs_error'=abs(`p_n'-`p_n_1')
			quietly: su `abs_error'
			scalar `distance'=r(max)
			if `iter'>0{
				if "`display'"=="display"{
					display `iter' "--" `distance'
				}
			}
			scalar `iter'=`iter'+1
			quietly: drop `temp_sum_c' `temp_sum_p' `p_n_1' `c_n_1' `sum_c_n' `sum_p_n' `p_n_temp'
			rename `c_n' `c_n_1'
			rename `p_n' `p_n_1'
		}
		quietly: gen `gen'_FCM_country=`c_n_1'
		quietly: gen `gen'_FCM_product=`p_n_1'
		local `keepvar_fcm' "`gen'_FCM_country `gen'_FCM_product"
	}
	********* Modified-Fitness-Complexity Method ********
	local `keepvar_mfcm' ""
	if "`mfcm'"=="mfcm"{
		if "`display'"=="display"{
			disp "Compute Complexity Index with the Modified-Fitness-Complexity Method (MFCM)"
		}
		** init **
		quietly: replace `abs_error'=.
		quietly: replace `c_n_1'=1
		quietly: replace `p_n_1'=1
		scalar `distance'=100000000.0
		scalar `iter'=0
		** begin iteration
		while `distance'>`tol' & `iter'<=`maxiter'{
			quietly: gen `temp_sum_c'=`M'*`p_n_1'
			quietly: gen `temp_sum_p'=`M'*(`num_product'-`c_n_1')
			quietly: egen `c_n'     =total(`temp_sum_c'), by(`country' `time')
			quietly: egen `p_n_temp'=total(`temp_sum_p'), by(`product' `time')
			quietly: gen `p_n'=1/`p_n_temp'
			** normalize
			quietly: egen `sum_c_n'=total(`c_n'), by(`time')
			quietly: replace `c_n'=`c_n'/`sum_c_n'*`num_country'*`num_product'
			quietly: egen `sum_p_n'=total(`p_n'), by(`time')
			quietly: replace `p_n'=`p_n'/`sum_p_n'*`num_country'*`num_product'
			** distance between iterates
			quietly: replace `abs_error'=abs(`p_n'-`p_n_1')
				quietly: su `abs_error'
			scalar `distance'=r(max)
			if `iter'>0{
				if "`display'"=="display"{
					display `iter' "--" `distance'
				}
			}
			scalar `iter'=`iter'+1
			quietly: drop `temp_sum_c' `temp_sum_p' `p_n_1' `c_n_1' `sum_c_n' `sum_p_n' `p_n_temp'
			rename `c_n' `c_n_1'
			rename `p_n' `p_n_1'
		}
		quietly: gen `gen'_MFCM_country=`c_n_1'
		quietly: gen `gen'_MFCM_product=`p_n_1'
		local `keepvar_mfcm' "`gen'_MFCM_country `gen'_MFCM_product"
	}
	// finally, merge file
	quietly: keep `country' `time' `product' ``keepvar_mreig'' ``keepvar_mr'' ``keepvar_fcm'' ``keepvar_mfcm''
	tempfile mergefile
	quietly: save `mergefile', replace
	restore
	quietly: merge 1:1 `country' `time' `product' using `mergefile'
	quietly: drop if _merge==2
	quietly: drop _merge
end
