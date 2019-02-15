options nocenter ls=85 ps=55;
data cstat;INFILE '\\Client\C$\Users\lennonj\Desktop\CCF\arima_stoich.csv' lrecl = 500 dlm=',' firstobs=2 missover;
INPUT sample time BNL2 PNL2 BNL3 PNL3 BNL5 PNL5 BPL2 PPL2 BPL4 PPL4 BPL5 PPL5;

/*Log10-transform abundance data */
BNL2=LOG10(BNL2);
PNL2=LOG10(PNL2);
BNL3=LOG10(BNL3);
PNL3=LOG10(PNL3);
BNL5=LOG10(BNL5);
PNL5=LOG10(PNL5);
BPL2=LOG10(BPL2);
PPL2=LOG10(PPL2);
BPL4=LOG10(BPL4);
PPL4=LOG10(PPL4);
BPL5=LOG10(BPL5);
PPL5=LOG10(PPL5);

/*PROC PRINT data=cstat;	
run;*/

/*Great description of PROC ARIMA here: 
http://support.sas.com/documentation/cdl/en/etsug/63939/HTML/default/viewer.htm#etsug_arima_sect004.htm

/*identify: produces series of diagnostics to test whether data are stationary
including autocorr plots, where you're looking for rapid decay in ACF, 
along with diagnostics, where larger p-values indicate stationarity

/*Also generates white noise test: is there data to model? Yes, if <0.05*/

/* `p` statement for auotregressive (AR) model (e.g., AR[1])
predicts the change in phage as an average change, plus some fraction of the previous change, 
plus a random error */

/* `q` statement for moving average (MA) model 

/* when and `p` and `q` are both used, you are creating a mixed autoregressive and moving-average model*/
/* An ARMA(1,1) model predicts the change in phage as an average change, plus some fraction of the previous change,
/* plus a random error, plus some fraction of the random error in the preceding period.*/

proc arima data=cstat; 
identify var=BPL5;
estimate p=1 q=1;

identify var=PPL5;
estimate p=1 q=1;

identify var=BPL5 crosscorr=PPL5 outcov=covEstimates nlag=8;
proc print data=covEstimates;
run;	
quit;
/*Clear results/output */
/* Highlight text, right click, and choose submit selection*/
dm 'odsresults; clear';

