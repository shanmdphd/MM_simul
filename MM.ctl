$PROB Breast Cancer Overall Survival Modeling with Right Censoring
;$DES Approach, Proportional Hazard Model, Exponential Distribution of Hazard (Constant Hazard over Time)
;OCS=0, death (event)
;OCS=1, alive (right censoring)
;TIME, year

$INPUT ID TIME DV REP

$DATA D:\AMC\Research\In_Vitro_DDI_Simulation\MM\MM_sim_data.csv IGNORE=@

$SUBR ADVAN=6 TOL=6

$MODEL
 COMP=(CONC)

$PK
 IF (NEWIND.LE.1) IBSL = DV  		        ; Concentration at TIME=0
 A_0(1) = IBSL

 VM   = THETA(1)*EXP(ETA(1))            ; Vmax
 KM   = THETA(2)                        ; Km
 

$DES
                
 DADT(1) = - VM*A(1)/(KM+A(1))           ; A(1), Regarded as concentration

$ERROR

  IPRED = A(1)
  Y = IPRED + EPS(1) + IPRED*EPS(2)

$THETA
 (0, 100)
 (0, 500)

$OMEGA 0 FIX

$SIGMA 0.1 0.1

$ESTIMATION SORT MAXEVAL=9999 PRINT=2 METHOD=COND SLOW INTER MSFO=MM.MSF
$COVARIANCE PRINT = E; MATRIX=S
$TABLE ID TIME IBSL MDV IPRED CWRES
       FILE=MM.FIT NOPRINT ONEHEADER
$TABLE ID ETA(1)
       FILE=MM.PAR NOPRINT ONEHEADER FIRSTONLY NOAPPEND
$TABLE ID VM KM
       FILE=MM.IPK NOPRINT ONEHEADER FIRSTONLY NOAPPEND

