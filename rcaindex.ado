*! version 14.0 Author : Jichun Si
* www.sijichun.pro
* The Revealed comparative advantage index

cap program drop rcaindex
program define rcaindex
	version 14
	syntax varlist(max=1), country(varname) product(varname) time(varname) [,gen(name)]
	tempvar _total_export _prod_export _country_export check_duplicates
	** check duplicates
	quietly: duplicates tag `country' `time' `product', gen(`check_duplicates')
	quietly: su `check_duplicates'
	if r(max)>0{
		di as err "Please check duplicates for country, time and products."
		exit 198
	}
	** compute
	egen `_total_export'=sum(`varlist'), by(`time')
	egen `_prod_export'=sum(`varlist'), by(`product' `time')
	egen `_country_export'=sum(`varlist'), by(`time' `country')
	tempvar _fm _fz
	gen `_fm'=`varlist'/`_country_export'
	gen `_fz'=`_prod_export'/`_total_export'
	if "`gen'"=="" {
		local gen "rcaindex"
	}
	gen `gen'=`_fm'/`_fz'
	label variable `gen' "revealed comparative advantage index"
end
