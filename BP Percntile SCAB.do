*Make sure the open STATA file to be analyzed has the following variables 
*      id, sex, age, height, sysbp, diask5.

 *      sex with value: 1 for male; 2 for female.
 *     age is in years.
 *     height is in CENTIMETERS.

 *    Note: The records may (a) list the above variables in any order;
 *                      and (b) contain missing data for 
 *                               either sysbp or diask5;
 *                      and (c) contain other variables.

 *           The percentiles for systolic and diastolic(?) are 
 *           provided as integer in the new variable names:
 *           syspct and diaspct.

 * This program calculates BP percentiles for children ages 0-17.99 years of age.


 
*ageflg: if age is an integer then                     */
*               set ageflg=1: +6.5 months when the age is     */
*               converted to months in the program;           */
*               if age is with decimal points then            */
*               set ageflg=2: +0.5 months when the age is     */
*               converted to months in the program  

use "C:\Users\ggorman\Documents\Research\BP Norms do File\sample_data.dta", clear

*make year a byte in case it's in integer form
recast byte age, force

gen agemon=round(age*12,1.0)+0.5

*merge main file with height file and calculate

*first merge with non infant height file and discard _merge==2
merge m:m agemon sex using "C:\Users\ggorman\Documents\Research\BP Norms do File\ht_zscore.dta"
drop if _merge==2
drop _merge
foreach var of varlist l m s {
	rename `var' `var'2
}

merge m:m agemon sex using "C:\Users\ggorman\Documents\Research\BP Norms do File\ht_zscore_inf.dta"
drop if _merge==2
drop _merge

*put all coefficients back into same variables
foreach var of varlist l m s {
	replace `var'=`var'2 if `var'==.
	drop `var'2
}

*greg, the exponents sign is prob wrong, substituted it fir **. check SAS code

gen zschtnw=(((height/m)^l)-1)/(l/s)
gen zHt=height/m
replace zHt=zHt^l
replace zHt=zHt-1
replace zHt=zHt/(l*s)
gen quant=normal(zHt)
gen quantile=(round(quant*100))/100
drop quant
gen ht1 =m*(-3.09*l*s+1)^(1/l)
gen ht2 =m*(3.09*l*s+1)^(1/l)

gen outlier=0
replace outlier=1 if height<ht1 | height>ht2

*set parameters of the 50th percentile SBP for males
*the nomeclature is t-then the knot # [5th, 27.5, 50th, 72.5th, 95th]

gen t1m=107.8 if sex==1
gen t2m=140.0 if sex==1
gen t3m=154.5 if sex==1
gen t4m=166.4 if sex==1
gen t5m=179.1 if sex==1
gen ta1m=5.06 if sex==1
gen ta2m=10.79 if sex==1
gen ta3m=13.22 if sex==1
gen ta4m=14.51 if sex==1
gen ta5m=17.30 if sex==1
gen tb1m=15.0 if sex==1
gen tb2m=8.9 if sex==1
gen tb3m=50.375 if sex==1
gen tb4m=112.684 if sex==1
gen tb5m=250.04 if sex==1

*set parameters for 50th percentile SBP females

replace t1m=106.7 if sex==2
replace t2m=140.7 if sex==2
replace t3m=154.0 if sex==2
replace t4m=160.5 if sex==2
replace t5m=168.9 if sex==2
replace ta1m=5.00 if sex==2
replace ta2m=10.70 if sex==2
replace ta3m=13.16 if sex==2
replace ta4m=14.51 if sex==2
replace ta5m=17.33 if sex==2
replace tb1m=6.701 if sex==2
replace tb2m=16.438 if sex==2
replace tb3m=46.80 if sex==2
replace tb4m=84.46 if sex==2
replace tb5m=203.608 if sex==2


*rename to merge on BP type

sort sex quantile
merge m:m sex quantile using "C:\Users\ggorman\Documents\Research\BP Norms do File\quantreg_coef.dta"
keep if (_merge==3 & type=="sys") | _merge==1

pause

*Create systolic BP
	
*was in SAS code: gen x=height

gen H=height

*was in SAS code: gen y=age

gen G=age

*replace w=(age-10)*(height-150) if sex==1
*replace w=(age-10)*(height-147) if sex==2

gen HG=(age-10)*(height-150) if sex==1
replace HG=(age-10)*(height-147) if sex==2


*do the regression


replace SBP95=b_0 + (x*bHt1 + x2s*bHt2 + x3s*bHt3 + x4s*bHt4) + (y*bage1 + y2s*bage2 + y3s*bage3 + y4s*bage4) + (w*bAgeHt1 + w2s*bAgeHt2 + w3s*bAgeHt3 + w4s*bAgeHt4)
  
	foreach num of numlist 1:99 {
		gen fxsys`num'= b0sys{`num'}+b1sys{i}*x1s+b2sys{i}*x2s+b3sys{i}*x3s+b4sys{i}*x4s
               +ba1sys{i}*y+ba2sys{i}*y2s+ba3sys{i}*y3s+ba4sys{i}*y4s
               +bb1sys{i}*w+bb2sys{i}*w2s+bb3sys{i}*w3s+bb4sys{i}*w4s;









pause
*may not need to do below since the file already calculates the coefficients for you; but where are the inputs

**********************************************************************
*generate values for the SBP coefficients (which are x) for Ht (which is H) using t1-5 
************************************************************************

* SAS Code: if x-t1m lt 0 then x2a=0;  else x2a=x-t1m;
gen x2a=H-t1m if type=="sys"
replace x2a=0 if H-t1m<0
	
*SAS Code: if x-t4m lt 0 then x2b=0;  else x2b=x-t4m;
gen x2b=H-t4m if type=="sys"
replace x2b=0 if H-t4m<0

*SAS Code: if x-t5m lt 0 then x2c=0;  else x2c=x-t5m;
gen x2c=H-t5m if type=="sys"
replace x2c=0 if (H-t5m)<0

*SAS Code: if x-t2m lt 0 then x3a=0;  else x3a=x-t2m;
gen x3a=H-t2m if type=="sys"
replace x3a=0 if H-t2m<0

*SAS Code: if x-t3m lt 0 then x4a=0;  else x4a=x-t3m;
gen x4a=H-t3m if type=="sys"
replace x4a=0 if H-t3m<0

*SAS Code: x2=x2a**3-x2b**3*(t5m-t1m)/(t5m-t4m)+x2c**3*(t4m-t1m)/(t5m-t4m);
gen x2=x2a^3-x2b^3*(t5m-t1m)/(t5m-t4m)+x2c^3*(t4m-t1m)/(t5m-t4m)
		
*SAS Code: x3=x3a**3-x2b**3*(t5m-t2m)/(t5m-t4m)+x2c**3*(t4m-t2m)/(t5m-t4m);
gen x3=x3a^3-x2b^3*(t5m-t2m)/(t5m-t4m)+x2c^3*(t4m-t2m)/(t5m-t4m)
	
*SAS Code: x4=x4a**3-x2b**3*(t5m-t3m)/(t5m-t4m)+x2c**3*(t4m-t3m)/(t5m-t4m);
gen x4=x4a^3-x2b^3*(t5m-t3m)/(t5m-t4m)+x2c^3*(t4m-t3m)/(t5m-t4m)
	
gen x1s=H	
replace x2s=x2/100
replace x3s=x3/100
replace x4s=x4/100

**********AGE************************************************************
*generate values for the SBP coefficients (which are y) for Age (which is G) using ta1-5 
************************************************************************

    
*SAS Code: if y-ta1m lt 0 then y2a=0;  else y2a=y-ta1m;
gen y2a=G-ta1m
replace y2a=0 if G-ta1m<0
	
	
*SAS Code: if y-ta4m lt 0 then y2b=0;  else y2b=y-ta4m;
gen y2b=G-ta4m
replace y2b=0 if G-ta4m<0
    
*SAS Code: if y-ta5m lt 0 then y2c=0;  else y2c=y-ta5m;
gen y2c=G-ta5m
replace y2c=0 if G-ta5m<0

*SAS Code: if y-ta2m lt 0 then y3a=0;  else y3a=y-ta2m;
gen y3a=G-ta2m
replace y3a=0 if G-ta2m<0

*SAS Code: if y-ta3m lt 0 then y4a=0;  else y4a=y-ta3m;
gen y4a=G-ta3m
replace y4a=0 if G-ta3m<0

*SAS Code: y2=y2a**3-y2b**3*(ta5m-ta1m)/(ta5m-ta4m)+y2c**3*(ta4m-ta1m)/(ta5m-ta4m);
gen y2=y2a^3-y2b^3*(ta5m-ta1m)/(ta5m-ta4m)+y2c^3*(ta4m-ta1m)/(ta5m-ta4m)
		
*SAS Code: y3=y3a**3-y2b**3*(ta5m-ta2m)/(ta5m-ta4m)+y2c**3*(ta4m-ta2m)/(ta5m-ta4m);
gen y3=y3a^3-y2b^3*(ta5m-ta2m)/(ta5m-ta4m)+y2c^3*(ta4m-ta2m)/(ta5m-ta4m)
		
*SAS Code: y4=y4a**3-y2b**3*(ta5m-ta3m)/(ta5m-ta4m)+y2c**3*(ta4m-ta3m)/(ta5m-ta4m);
gen y4=y4a^3-y2b^3*(ta5m-ta3m)/(ta5m-ta4m)+y2c^3*(ta4m-ta3m)/(ta5m-ta4m)

gen y1s=G
replace y2s=y2/100
replace y3s=y3/100
replace y4s=y4/100



**********AGE x Ht INTERACTION************************************************************
*generate values for the SBP coefficients (which are x) for Age-Ht interaction (which is HG) using tb1-5. Note final denominator is 100^2 
************************************************************************
    
*SAS Code: if w-tb1m lt 0 then w2a=0;  else w2a=w-tb1m;
gen w2a=HG-tb1m
replace w2a=0 if HG-tb1m<0
	
*SAS Code: if w-tb4m lt 0 then w2b=0;  else w2b=w-tb4m;
gen w2b=HG-tb4m
replace w2b=0 if HG-tb4m<0
	
*SAS Code: if w-tb5m lt 0 then w2c=0;  else w2c=w-tb5m;
gen w2c=HG-tb5m
replace w2c=0 if HG-tb5m<0

*SAS Code: if w-tb2m lt 0 then w3a=0;  else w3a=w-tb2m;
gen w3a=HG-tb2m
replace w3a=0 if HG-tb2m<0

*SAS Code: if w-tb3m lt 0 then w4a=0;  else w4a=w-tb3m;
gen w4a=HG-tb3m
replace w4a=0 if HG-tb3m<0


*SAS Code: w2=w2a**3-w2b**3*(tb5m-tb1m)/(tb5m-tb4m)+w2c**3*(tb4m-tb1m)/(tb5m-tb4m);
gen w2=w2a^3-w2b^3*(tb5m-tb1m)/(tb5m-tb4m)+w2c^3*(tb4m-tb1m)/(tb5m-tb4m)
	
*SAS Code: w3=w3a**3-w2b**3*(tb5m-tb2m)/(tb5m-tb4m)+w2c**3*(tb4m-tb2m)/(tb5m-tb4m);
gen w3=w3a^3-w2b^3*(tb5m-tb2m)/(tb5m-tb4m)+w2c^3*(tb4m-tb2m)/(tb5m-tb4m)
	
*SAS Code: w4=w4a**3-w2b**3*(tb5m-tb3m)/(tb5m-tb4m)+w2c**3*(tb4m-tb3m)/(tb5m-tb4m);
gen w4=w4a^3-w2b^3*(tb5m-tb3m)/(tb5m-tb4m)+w2c^3*(tb4m-tb3m)/(tb5m-tb4m)
	
gen w1s=HG	
replace w2s=w2/100^2
replace w3s=w3/100^2
replace w4s=w4/100^2

	
	
	
	

	

	
    
	
	
	
	foreach num of numlist 1:99 {
		gen fxsys`num'= b0sys{`num'}+b1sys{i}*x1s+b2sys{i}*x2s+b3sys{i}*x3s+b4sys{i}*x4s
               +ba1sys{i}*y+ba2sys{i}*y2s+ba3sys{i}*y3s+ba4sys{i}*y4s
               +bb1sys{i}*w+bb2sys{i}*w2s+bb3sys{i}*w3s+bb4sys{i}*w4s;
      difdias{i} =abs(diask5-fxdias{i})
	
	
	do i=1 to 99;
	fxdias{i}=b0sys{i}+b1sys{i}*x+b2sys{i}*x2s+b3sys{i}*x3s+b4sys{i}*x4s
               +ba1sys{i}*y+ba2sys{i}*y2s+ba3sys{i}*y3s+ba4sys{i}*y4s
               +bb1sys{i}*w+bb2sys{i}*w2s+bb3sys{i}*w3s+bb4sys{i}*w4s;
      difdias{i} =abs(diask5-fxdias{i});
    end;

    pdias=min(of difdias1-difdias99);
  run;

  data all2;
    set all2;

    array difdias{*} difdias1-difdias99;
    do i=1 to 99;
      if pdias=difdias{i} then diaspct=i;
    end;

    if diask5=. then diaspct=.; 

    if outlier=1 then diaspct=.;
	
	


	
	
	*this may not be right rename coefficients
	gen b1sys=x
	gen b2sys=x2s
	gen b3sys=x3s
	gen b4sys=x4s
	
	gen ba1sys=y
	gen ba2sys=y2s
	gen ba3sys=y3s
	gen ba4sys=y4s
	
	gen bb1sys=w
	gen bb2sys=w2s
	gen bb3sys=w3s
	gen bb4sys=w4s
	
	*his may not be right
	gen fxsys=Intercept+b1sys*x+b2sys*x2s+b3sys*x3s+b4sys*x4s+ba1sys*y+ba2sys*y2s+ba3sys*y3s+ba4sys*y4s+bb1sys*w+bb2sys*w2s+bb3sys*w3s+bb4sys*w4s if type=="sys"
	gen difsys=abs(sysbp-fxsys)


* Create diastolic BP

 if x-t1m lt 0 then x2a=0;  else x2a=x-t1m;
    if x-t4m lt 0 then x2b=0;  else x2b=x-t4m;
    if x-t5m lt 0 then x2c=0;  else x2c=x-t5m;
    x2=x2a**3-x2b**3*(t5m-t1m)/(t5m-t4m)+x2c**3*(t4m-t1m)/(t5m-t4m);
    if x-t2m lt 0 then x3a=0;  else x3a=x-t2m;
    x3=x3a**3-x2b**3*(t5m-t2m)/(t5m-t4m)+x2c**3*(t4m-t2m)/(t5m-t4m);
    if x-t3m lt 0 then x4a=0;  else x4a=x-t3m;
    x4=x4a**3-x2b**3*(t5m-t3m)/(t5m-t4m)+x2c**3*(t4m-t3m)/(t5m-t4m);
    x2s=x2/100;
    x3s=x3/100;
    x4s=x4/100;
    if y-ta1m lt 0 then y2a=0;  else y2a=y-ta1m;
    if y-ta4m lt 0 then y2b=0;  else y2b=y-ta4m;
    if y-ta5m lt 0 then y2c=0;  else y2c=y-ta5m;
    y2=y2a**3-y2b**3*(ta5m-ta1m)/(ta5m-ta4m)+y2c**3*(ta4m-ta1m)/(ta5m-ta4m);
    if y-ta2m lt 0 then y3a=0;  else y3a=y-ta2m;
    y3=y3a**3-y2b**3*(ta5m-ta2m)/(ta5m-ta4m)+y2c**3*(ta4m-ta2m)/(ta5m-ta4m);
    if y-ta3m lt 0 then y4a=0;  else y4a=y-ta3m;
    y4=y4a**3-y2b**3*(ta5m-ta3m)/(ta5m-ta4m)+y2c**3*(ta4m-ta3m)/(ta5m-ta4m);
    y2s=y2/100;
    y3s=y3/100;
    y4s=y4/100;
    if w-tb1m lt 0 then w2a=0;  else w2a=w-tb1m;
    if w-tb4m lt 0 then w2b=0;  else w2b=w-tb4m;
    if w-tb5m lt 0 then w2c=0;  else w2c=w-tb5m;
    w2=w2a**3-w2b**3*(tb5m-tb1m)/(tb5m-tb4m)+w2c**3*(tb4m-tb1m)/(tb5m-tb4m);
    if w-tb2m lt 0 then w3a=0;  else w3a=w-tb2m;
    w3=w3a**3-w2b**3*(tb5m-tb2m)/(tb5m-tb4m)+w2c**3*(tb4m-tb2m)/(tb5m-tb4m);
    if w-tb3m lt 0 then w4a=0;  else w4a=w-tb3m;
    w4=w4a**3-w2b**3*(tb5m-tb3m)/(tb5m-tb4m)+w2c**3*(tb4m-tb3m)/(tb5m-tb4m);
    w2s=w2/100**2;
    w3s=w3/100**2;
    w4s=w4/100**2;

	
	foreach num of numlist 1:4 {
		gen b`num'sys=0
		gen ba`num'sys=0
		gen bb`num'sys=0
	}
	
	
	
	do i=1 to 99;
      fxdias{i}=b0sys{i}+b1sys{i}*x+b2sys{i}*x2s+b3sys{i}*x3s+b4sys{i}*x4s
               +ba1sys{i}*y+ba2sys{i}*y2s+ba3sys{i}*y3s+ba4sys{i}*y4s
               +bb1sys{i}*w+bb2sys{i}*w2s+bb3sys{i}*w3s+bb4sys{i}*w4s;
      difdias{i} =abs(diask5-fxdias{i});
    end;

    pdias=min(of difdias1-difdias99);
  run;

  data all2;
    set all2;

    array difdias{*} difdias1-difdias99;
    do i=1 to 99;
      if pdias=difdias{i} then diaspct=i;
    end;

    if diask5=. then diaspct=.; 

    if outlier=1 then diaspct=.;
	
	


*merge BP coefficient file and calculate


use "C:\Users\ggorman\Documents\Research\BP Norms do File\quantreg_coef.dta", clear

gen x=b1sys 
gen x2s=b2sys 
gen x3s=b3sys 
gen x4s=b4sys
gen y=ba1sys 
gen y2s=ba2sys 
gen y3s=ba3sys 
gen y4s=ba4sys
gen w=bb1sys 
gen w2s=bb2sys 
gen w3s=bb3sys 
gen w4s=bb4sys


