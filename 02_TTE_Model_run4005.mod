;; 1. Based on: run4004
;; 2. Description: Time To First Recurrence, Gompertz Hazard
;; x1. Author: Och

;; 3. Label:

$PROBLEM TTE
$INPUT ID USUBJID=DROP TIME DV OBS CMT MDV EVID AMT RATE TYPE FE EVENT1 AGE SEX TIME_PRIOR_RAND ANOGEPSD SMOKER_NO CONTRACEPTIVE_DRUG ALBUMIN PROTEIN ALBUMIN_BASE PROTEIN_BASE TREATMENT_NO HEIGHTBL WEIGHTBL ET1 ET2 PK_UNKNOWN
$DATA     ..\DATASET\HDIT101_TTE_V01.csv IGNORE=@

;Sim_start : remove from simulation model   
IGNORE=(FE.EQ.2) 
;Sim_end

$SUBROUTINES ADVAN6 TOL=5

$MODEL
COMP = (CENTRAL1)
COMP = (PERI1)

COMP = (ABS2)
COMP = (CENTRAL2)

COMP = (dummy)
COMP = (RECURRENCE)

COMP = (TRANSIT1)
COMP = (TRANSIT2)
COMP = (TRANSIT3)
COMP = (TRANSIT4)
COMP = (TRANSIT5)
COMP = (TRANSIT6)
COMP = (TRANSIT7)
COMP = (TRANSIT8)
COMP = (TRANSIT9)
COMP = (TRANSIT10)
COMP = (TRANSIT11)
COMP = (TRANSIT12)

$PK
IF (NEWIND.NE.2) TP=0
;----HDIT PK----
STUDY=2
CL	= 0.01 *EXP(ET1)*(WEIGHTBL/72)**0.654
V	= 3.03*EXP(ET2)*((HEIGHTBL/174)**2.12)*(1.13**(STUDY-1))
Q   = 0.013
VP  = 2.71*(0.728**SEX)*(1+0.0178*(TIME_PRIOR_RAND-2.8598))
S1	= V 
dummy=ETA(1)

;----Valaciclovir PK Vezina et al,2013----
VAL_KA=0.28         ;{typical value of Ka}
VAL_CL=49.9         ;{typical value of CL/F}
VAL_VC=74.1 

VAL_KEL=VAL_CL/VAL_VC

;----Hazard PD----

MTT = THETA(6)   ;Transit compartment rate constant

ktr = 12/MTT         ;Baseline KTR
A_0(7)= 1           ;Set transit compartment1 initial value
A_0(8)= 1           ;Set transit compartment2 initial value
A_0(9)= 1           ;Set transit compartment3 initial value
A_0(10)= 1          ;Set transit compartment4 initial value
A_0(11)= 1          ;Set transit compartment5 initial value
A_0(12)= 1          ;Set transit compartment6 initial value
A_0(13)= 1          ;Set transit compartment7 initial value
A_0(14)= 1          ;Set transit compartment8 initial value
A_0(15)= 1          ;Set transit compartment9 initial value
A_0(16)= 1          ;Set transit compartment10 initial value
A_0(17)= 1          ;Set transit compartment11 initial value
A_0(18)= 1          ;Set transit compartment12 initial value

EMAX=THETA(4)
EC50=exp(THETA(5))

MTDIFF=1
MTIME(1)=THETA(3)
COVARIATE1=MPAST(1)

ANOGEPSD_COV=EXP(THETA(7)*ANOGEPSD)
T2R_COV=(TIME_PRIOR_RAND/5)**THETA(8)

LAM0=THETA(1)*COVARIATE1*ANOGEPSD_COV*T2R_COV
SHP0=exp(THETA(2))*(-1)

LAM		= LAM0/1000	;scale factor
SHP		= (SHP0/100000)	;shape factor

$DES
;----HDIT PK----
DADT(1) = - CL/V*A(1) - Q/V*A(1) + Q/VP*A(2)
DADT(2) =               Q/V*A(1) - Q/VP*A(2)

;----Valaciclovir PK----
DADT(3) = -VAL_KA * A(3) 
DADT(4) =  VAL_KA * A(3)  - VAL_KEL * A(4) 	             

;----dummy---
DADT(5) = 0

;----Hazard----

DELAY	=   1E-6
LAMI	= 	LAM
SHPI	=	SHP

DADT(6) = (LAMI*EXP(SHPI*(T+DELAY)))*A(18)

;----Transit----
DRUG_CONC=A(1)/V
DRUG  = 1-(EMAX*DRUG_CONC)/(EC50+DRUG_CONC)  
   
DADT(7) =  ktr*(DRUG)-ktr*A(7)
DADT(8) =  ktr*A(7)-ktr*A(8)
DADT(9) =  ktr*A(8)-ktr*A(9)
DADT(10) =  ktr*A(9)-ktr*A(10)
DADT(11) =  ktr*A(10)-ktr*A(11)
DADT(12) =  ktr*A(11)-ktr*A(12)
DADT(13) =  ktr*A(12)-ktr*A(13)
DADT(14) =  ktr*A(13)-ktr*A(14)
DADT(15) =  ktr*A(14)-ktr*A(15)
DADT(16) =  ktr*A(15)-ktr*A(16)
DADT(17) =  ktr*A(16)-ktr*A(17)
DADT(18) =  ktr*A(17)-ktr*A(18)

$ERROR
;----Hazard PD----

LAMII	= 	LAM
SHPII	=	SHP
HAZNOW= (LAMII*EXP(SHPII*(TIME+1E-6)))*A(18)

IF(NEWIND.NE.2) OLDCHZ=0                 ;Initialization of OLDCHZ for each ID
CHZ =A(6)-OLDCHZ                         ;Cumulative hazard since latest event
;Sim_start : add/remove for simulation
;OLDCHZ=A(6)                             ;Initialization of OLDCHZ for each ID
IF(DV.NE.0.AND.CMT.EQ.6) OLDCHZ=A(6)     ;Initialization of OLDCHZ for each ID
;Sim_end
SUR=EXP(-CHZ)                   	     ;Survival
IF(DV.EQ.0) Y=SUR                        ;Non-event DV=0
IF(DV.NE.0) Y=SUR*HAZNOW                 ;event DV=1

TP=TIME

IF(ICALL.EQ.4) THEN                      ; for simulation
   CALL RANDOM (2,R)
   DV=0
   RTTE = 0
   IF(TIME.EQ.4320) RTTE = 1             ; for the censored observation at days 180 (4320 hrs)
   IF(R.GT.SUR) RTTE = 1
   IF(R.GT.SUR) DV=1
 ENDIF


$THETA
(0, 0.414)     ; Lambda
(0, 2.54)      ; Shape
(10, 120)FIX   ; Mtime
(0, 1) FIX     ; Emax (HDIT101)
(0, 4.67)      ; Ec50 (HDIT101)
(0, 846)       ; MTT  (HDIT101)
(-2, 0.103)    ; Herpes episodes
(-2, -0.292)   ; Disease duration

$OMEGA  
 0 FIX  ; dummy

;Sim_start : add/remove for simulation
;$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION SUB=100

$ESTIM MAXEVAL=9990 METHOD=1 LAPLACE LIKE PRINT=1 MSFO=msfb4005SIMPRO SIGL=9 NSIG=3 NOABORT
$COV PRINT=E
;Sim_end

$TABLE ID TIME ETAS(1:LAST)CMT HAZNOW COVARIATE1 DRUG TREATMENT_NO AGE SEX HEIGHTBL WEIGHTBL TIME_PRIOR_RAND EVID NOPRINT ONEHEADER FILE=sdtab4005
