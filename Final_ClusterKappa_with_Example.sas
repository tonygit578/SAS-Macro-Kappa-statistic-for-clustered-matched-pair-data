/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
### SAS macro to calculate the Kappa statistic and its variance under clustered data
 		Code to create the results in Examples section of the manuscript;
### "Kappa statistic for clustered matched-pair data";
### for Statistics in Medicine (in press);
### by Zhao Yang, Ming Zhou;
### Version Date: 23May2013
** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;
/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
##### Revision History
On 30Jan2014
1. remove the macro parameter trt1V	and trt2V which are not really needed in the run
2. remove the restriction that the macro parameters clusV and unitV need to be numeric
3. Update the macro parameter description

On 04Feb2014
Change the variable name resp1 to resp1__, resp2 to resp2__ to avoid the name conflict

On 24Sep2014
Based on the comments from Dr Sophie Vanbelle, update the formula by changing
2 # PV_ * PV_` * BigOmega * P_V * P_V`
to
2 # P_V * P_V` * BigOmega * PV_ * PV_` 
in the component V22, however, there are no changes to the results

On 25Sep2014
More update on the the formula calculation
2 # P_V * P_V` * BigOmega * PV_ * PV_` 
to
PV_ * P_V` * BigOmega *  PV_ * P_V` + P_V * PV_` * BigOmega *  P_V * PV_`
in the component V22 
** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;

%macro ClusterKappa (dsin = , /* the input dataset */
					 clusV = , /* the cluster variable*/
					 unitV = , /* the unit variable*/
					 Resp1V =, /* Response variable for procedure 1, Row variable
					 			need to be numeric integer and > 0 */
					 Resp2V =  /* Response variable for procedure 2, Column Variable
								need to be numeric integer and > 0 */);
	proc sql noprint;
		create table clustDS as
			select distinct &clusV from &dsin
			order by &clusV; quit;
	data clustDS;
		set clustDS; clusterID + 1;
	proc sort data = &dsin;
		by &clusV &unitV;
	data resp_;
		merge &dsin clustDS;
		by &clusV;
	proc sort data = resp_;
		by clusterID &unitV;
	data F_Ana_DS (drop = &clusV &unitV rename = (clusterID = cluster));
		set resp_;
		by clusterID &unitV;
		if first.clusterID then subject = 0;
			subject + 1; 
	data anal_D;
		set F_Ana_DS (rename = (&Resp1V = resp1__ &Resp2V = resp2__ ));
		maxR = max(resp1__, resp2__);
	run;

	proc sql noprint;
		select max(maxR) into: N_Cat /*** to get the # of response category */ 
			from anal_D;
		select max(cluster) into: nclus /*** to get the # of clusters ***/
			from anal_D;
		create table count1 as  /** this is the calculate the marginal total for each cluster of treatment 1 */
			select distinct cluster, resp1__, count(*) as count1
				from anal_D
				group by cluster, resp1__;  
		create table count2 as  /** this is the calculate the marginal total for each cluster of treatment 2 */
			select distinct cluster, resp2__, count(*) as count2
				from anal_D
				group by cluster, resp2__; 
	quit; 

	*** this is to construct the dataset with 0 in each combination;
	data constru1; /*** generate a dataset with full combination of tabulation */
		do cluster = 1 to &nclus;
			do resp1__ = 1 to &N_Cat;
				count1 = 0; output;
			end;
		end;
	proc sort data = constru1;  
		by cluster resp1__;
	data count1;	/** merge the datasets to get the full combination of tabulation */
		merge constru1 count1 ;
		by cluster resp1__;

	data constru2;
		do cluster = 1 to &nclus;
			do resp2__ = 1 to &N_Cat;
				count2 = 0; output;
			end;
		end;
	proc sort data = constru2;
		by cluster resp2__;
	data count2; /** merge the datasets to get the full combination of tabulation */
		merge constru2 count2 ;
		by cluster resp2__;
	run; quit;

	proc summary data = count1; /*** contains the # of total subjects */
		class cluster;
		var count1;
		output out = total (drop = _freq_ where = (_TYPE_ = 0) ) sum = total ;
	data _null_;
		set total (keep = total);
		call symput("total", total);
	run;

	proc summary data = count1; /*** contains the # of total subjects for each cluster;*/
		class cluster;
		var count1;
		output out = n_k (drop = _freq_ where = (_TYPE_ = 1) ) sum = n_k ;
	proc sort data = n_k (drop = _type_) out = n_k;
		by cluster;
	run;  

	data concordA; 
		set anal_D;
		if resp1__ = resp2__ then freq = 1; else freq = 0;

	proc summary data = concordA; /** to calculate the concordant pair at the diagnoal of the k x k table;*/
		class cluster;
		var freq;
		output out = concord (drop = _freq_ where = (_TYPE_ = 1) ) sum = n_ii_Sk ;
	proc sort data = concord (drop = _type_) out = concord;
		by cluster;
	run;  

	%macro outputM(countds = );	
		%do i = 1 %to &N_Cat;
			%let seqD = %sysfunc( substr( &countds, 6, 1));
			data margin&seqD.&i (rename = (count&seqD = %if &seqD = 1 %then %do; n&i._k %end; %else %do; n_&i.k %end; ) drop = resp&seqD.__); 
				set &countds;
				if resp&seqD.__ = &i then output;
			proc sort data = margin&seqD.&i;
				by cluster;
			run;
		%end;
	%mend;
	%outputM(countds = Count1); /*** output the marginal total for treatment 1*/;
	%outputM(countds = Count2); /*** output the marginal total for treatment 2*/;

	data finalD;
		merge N_k 
			%do i = 1 %to &N_Cat;
				margin1&i margin2&i 
			%end;	concord;
		by cluster;
	data finalD;
		set finalD;
		BigN = &total;
	run;

	data finalD_A (drop = i);
		set finalD;
		array varlistA {*} %do i = 1 %to &N_Cat; n&i._k n_&i.k %end; n_ii_Sk;
		array varlistB {*} %do i = 1 %to &N_Cat; p&i._k p_&i.k %end; p_ii_Sk;
		do i = 1 to dim(varlistA);
			varlistB {i} = varlistA {i}/n_k;
		end;
		omega_k = n_k/BigN;
	run;
			
	*** calculation ussing the Kappa under the independence (non-clusterd situation);
	data fakestr;
		do resp1__ = 1 to &N_Cat;
			do resp2__ = 1 to &N_Cat;
				count = 0;	output;
			end;
		end;
	proc sort data = fakestr;
		by resp1__ resp2__;
	run;

	proc summary data = anal_D;
		class resp1__ resp2__;
		output out = interm_D (where = (_TYPE_ = 3) rename = (_FREQ_ = count) );
	proc sort data = interm_D (drop = _TYPE_) out = interm_D;
		by resp1__ resp2__;
	data A_Response;
		merge fakestr interm_D;
		by resp1__ resp2__;
	run;

	title2 "*** Part 1: Kappa statistics assuming independence *** ";	
	proc freq data = A_Response;
		tables resp1__ * resp2__/agree;
		weight count/zeros;
	run;

	proc iml;
		* Read data into IML ;
		use finalD_A;  read all ;
		kappOD = J( 1, 6, .);
		*** this is used to calculate the point estimation;
		%do i = 1 %to &N_Cat;
			  p&i._ = omega_k` * p&i._k;
			  p_&i = omega_k` * p_&i.k;
		%end;

		P_O = omega_k` * p_ii_Sk;
		P_E = P1_ * P_1  %do i = 2 %to &N_Cat; + P&i._ * P_&i %end; ;

		kappa = (P_O - P_E) / (1 - P_E); ** this is the calculated Kappa value;

		*** end of point estimation calculation;
		*** this part is for the calculation of the variance based on Delta method;
		PO_V = p_ii_Sk;

		P_V  = p1_k %do i = 2 %to &N_Cat; || p&i._k %end; ;
		PV_  = p_1k %do i = 2 %to &N_Cat; || p_&i.k %end; ;

		omega_k2 = omega_k # omega_k;
		BigOmega = ( I(&nclus) - omega_k * J(1, &nclus, 1) ) * diag(omega_k2) * ( I(&nclus) - J(&nclus, 1, 1) * t(omega_k) );
			  
		V11 = PO_V` * BigOmega * PO_V;
		V12 = PO_V` * BigOmega * (P_V * PV_` + PV_ * P_V`) * omega_k;
		V21 = V12;
		V22 = omega_k` * (PV_ * P_V` * BigOmega * P_V * PV_` + P_V * PV_` * BigOmega * PV_ * P_V` + 
							PV_ * P_V` * BigOmega *  PV_ * P_V` + P_V * PV_` * BigOmega *  P_V * PV_`) * omega_k;

		VMatrix = J( 2, 2, .);
		LVector = J( 1, 2, .);

		VMatrix[1,1] = V11;
		VMatrix[1,2] = V12;
		VMatrix[2,1] = V21;
		VMatrix[2,2] = V22;
		VMatrix = (&nclus # &nclus) # VMatrix/(&nclus - 1);
			  
		LVector[1] = 1 / (1 - P_E);
		LVector[2] = (P_O - 1) / ( (1 - P_E) # (1 - P_E) );

		ASE_Del = sqrt( LVector * VMatrix * LVector`/&nclus ); *** variance using Delta method;
			  
		lower = kappa - 1.96 # ASE_Del;
		upper = kappa + 1.96 # ASE_Del;
			 * print V11, V12, V21, V22 VMatrix Var_Del P_O kappa;

			  *** end of the calculation of the variance based on Delta method;

			  ******* this part is for pooling the final results information;
		
		title2 "*** Part 2: Kappa statistics for clustered matched-pair data *** ";	
		print P_O P_E kappa ASE_Del[label="Standard Error"],, 
			  lower[label="95% CI Lower limit"] upper[label="95% CI Upper limit"];
		quit;
%mend;

/***** Here is the example-1 ********/;
data response;
 	input Cluster Subject Trt2 Resp2 Trt1 Resp1;
cards;
1 1 1 0 2 0
1 2 1 0 2 1
1 3 1 0 2 1
2 1 1 1 2 1
2 2 1 1 2 1
2 3 1 0 2 1
3 1 1 1 2 1
3 2 1 1 2 1
3 3 1 1 2 1
4 1 1 1 2 1
5 1 1 1 2 1
5 2 1 1 2 1
5 3 1 0 2 1
6 1 1 1 2 1
6 2 1 1 2 1
6 3 1 1 2 1
6 4 1 1 2 1
7 1 1 1 2 1
7 2 1 1 2 1
7 3 1 1 2 1
8 1 1 1 2 1
8 2 1 1 2 1
9 1 1 1 2 1
9 2 1 1 2 0
10 1 1 1 2 1
11 1 1 1 2 1
11 2 1 1 2 1
11 3 1 0 2 0
12 1 1 1 2 1
12 2 1 1 2 1
13 1 1 1 2 1
13 2 1 1 2 1
13 3 1 1 2 1
14 1 1 1 2 1
14 2 1 1 2 1
15 1 1 0 2 1
15 2 1 0 2 1
16 1 1 1 2 1
16 2 1 1 2 1
16 3 1 0 2 0
17 1 1 1 2 1
17 2 1 1 2 1
17 3 1 0 2 0
18 1 1 1 2 1
18 2 1 1 2 1
18 3 1 0 2 1
19 1 1 1 2 1
19 2 1 1 2 1
20 1 1 1 2 1
21 1 1 1 2 1
21 2 1 1 2 1
;
run;

data response;
	set response;
	if Resp1 = 1 then Resp1 = 2; else Resp1 = 1;
	if Resp2 = 1 then Resp2 = 2; else Resp2 = 1;
run;

%ClusterKappa(dsin = response,  clusV = cluster,  unitV = subject,  Resp1V = Resp1, Resp2V = Resp2);
