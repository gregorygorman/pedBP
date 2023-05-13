*Make sure the open STATA file to be analyzed has the following variables 
*      id, sex, age, height, sysbp, diask5.

 *      sex with value: 1 for male; 2 for female.
 *     age is in years.
 *     height is in CENTIMETERS.

 *    Note: The records may (a) list the above variables in any order;
 *                      and (b) contain missing data for 
 *                               either sysbp or diask5;
                       and (c) contain other variables.

 *           The percentiles for systolic and diastolic(?) are 
 *           provided as integer in the new variable names:
 *           syspct and diaspct.

 * This program calculates BP percentiles for children ages 0-17.99 years of age.


 
ageflg: if age is an integer then                     */
/*               set ageflg=1: +6.5 months when the age is     */
/*               converted to months in the program;           */
/*               if age is with decimal points then            */
/*               set ageflg=2: +0.5 months when the age is     */
/*               converted to months in the program  
 
gen ageflg=0
replace ageflg=1 if age is integer
replace ageflg=2 if age is byte
 

gen agemon=round(age*12,1.0)+6.5 if ageflg==1
replace agemon=round(age*12,1.0)+0.5 if ageflg==2


*merge main file with height file and calculate

gen zschtnw=(((height/m)**l)-1)/(l*s)
gen probzht=probnorm(zschtnw);
gen ht1 =m*(-3.09*l*s+1)**(1/l)
gen ht2 =m*(3.09*l*s+1)**(1/l)

gen outlier=0
replace outlier=1 if height<ht1 or height>ht2


*merge BP coefficient file and cvalculate

gen x=b1sys x2s=b2sys x3s=b3sys x4s=b4sys
  	   y=ba1sys y2s=ba2sys y3s=ba3sys y4s=ba4sys
	   w=bb1sys w2s=bb2sys w3s=bb3sys w4s=bb4sys
	
    if sex='M' then sex1=1;
    else if sex='F' then sex1=2;
    drop sex;

