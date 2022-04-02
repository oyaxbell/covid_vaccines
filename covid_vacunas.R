## Effectiveness of a nation-wide COVID-19 vaccination program in Mexico against symptomatic COVID-19, hospitalizations, 
## and death: a retrospective analysis of national surveillance data
## Data Analysis: Omar Yaxmehen Bello-Chavolla (oyaxbell@yahoo.com.mx)
## Latest version of Analysis 01-Apr-2022
## Any question regarding analysis, please contact Omar Yaxmehen Bello-Chavolla

#### Dataset loading and cleaning ####
pacman::p_load(tidyverse, readr, lmerTest, lme4, lubridate, data.table, MatchIt, jtools,rms,bit64,iDOVE,stringi,optmatch,
               survival, parallel, MASS, flextable, officer, lmtest, sandwich, survminer, oddsratio, ggsci, ggpubr)
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/Vacunacion COVID-19/")
setwd("C:/Users/Investigacion/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/Vacunacion COVID-19")
base_final<-read_rds("base_final.rds")
efficacy<- function(x){
  ef<-(1-exp(coef(x)[1]))*100
  ul<-(1-exp(as.data.frame(confint(x)[1,1])))*100
  ll<-(1-exp(as.data.frame(confint(x)[1,2])))*100
  out<-cbind(round(ef,2), round(ll,2),round(ul,2))
  out1<-paste0(out[1]," (",out[2],"-",out[3],")")
  return(out1)
}

extract<-function(x){
  ef<-as.numeric(str_extract_all(x, "[[:digit:]]+\\.*[[:digit:]]*")[[1]][1])
  le<-as.numeric(str_extract_all(x, "[[:digit:]]+\\.*[[:digit:]]*")[[1]][2])
  ue<-as.numeric(str_extract_all(x, "[[:digit:]]+\\.*[[:digit:]]*")[[1]][3])
  df<-data.frame(ef,le,ue)
  return(df)
}

base_final$hosp_covid<-ifelse(base_final$covid==1 & base_final$TIPO_PACIENTE==2, 1, 0)
base_final$muerte_covid<-ifelse(base_final$covid==1 & base_final$muerte==1, 1, 0)
base_final$TIPO_VAC_COV[base_final$TIPO_VAC_COV=="Convidencia"]<-"CanSino"
base_final$DIABETES<-ifelse(base_final$DIABETES=="1" | base_final$DIABETES=="SI",1,0)
base_final$OBESIDAD<-ifelse(base_final$OBESIDAD=="1" | base_final$OBESIDAD=="SI",1,0)
base_final$months_from_vaccine<-round(base_final$months_from_vaccine)
base_final<- base_final %>%mutate(months=dplyr::recode(months_from_vaccine, "0"="1", "1"="1","2"="2", "3"="3", 
                                                       "4"="4", "5"="5", "6"=">6", "7"="7", "8"="8"))

base_final$edad_cat[base_final$EDAD<30]<-1
base_final$edad_cat[base_final$EDAD>=30 & base_final$EDAD<60]<-2
base_final$edad_cat[base_final$EDAD>=60]<-3
table(base_final$edad_cat)
base_final$edad_cat<-factor(base_final$edad_cat, labels = c("18-29 years", "30-59 years", ">60 years"))

#### Vaccine effectiveness CanSino ####
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) #%>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time),
                                       txtime = time_to_vaccine,
                                       fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                CanSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino,data= sdata)
cansino_e1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+EDAD+SEXO+comorb+strata(id),data= sdata)
cansino_e3<-efficacy(m1)

### Mortalidad ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time_death),
                                       txtime = time_to_vaccine,
                                       fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                CanSino = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(tstart, tstop, muerte)~CanSino, data=sdata%>% filter(covid==1))
cansino_d1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~CanSino+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1))
cansino_d3<-efficacy(m1)

### Hospitalization ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))

sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                CanSino = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(tstart, tstop, hosp_covid)~CanSino, data=sdata%>% filter(covid==1))
cansino_h1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~CanSino+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1))
cansino_h3<-efficacy(m1)

#### Vaccine effectiveness Sputnik ####
base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18)%>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time),
                                       txtime = time_to_vaccine,
                                       fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))
### Esquema incompleto ###
m1<-coxph(Surv(tstart, tstop, covid)~Sputnik, data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
sputnik_e1_p<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
sputnik_e3_p<-efficacy(m1)

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid)~Sputnik, data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
sputnik_e1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
sputnik_e3_c<-efficacy(m1)

### Mortalidad INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time_death),
                                       txtime = time_to_vaccine,
                                       fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, muerte)~Sputnik, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sputnik_d1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, muerte)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sputnik_d3_p<-efficacy(m1)

### Mortalidad COMPLETA ###
m1<-coxph(Surv(tstart, tstop, muerte)~Sputnik, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sputnik_d1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sputnik_d3_c<-efficacy(m1)

### Hospitalización INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, hosp_covid)~Sputnik, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sputnik_h1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, hosp_covid)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sputnik_h3_p<-efficacy(m1)

### Hospitalización COMPLETA ###
m1<-coxph(Surv(tstart, tstop, hosp_covid)~Sputnik, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sputnik_h1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~Sputnik+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sputnik_h3_c<-efficacy(m1)

#### Vaccine effectiveness AstraZeneca ####
base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")

tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time),
                                       txtime = time_to_vaccine,
                                       fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))
### Esquema incompleto ###
m1<-coxph(Surv(tstart, tstop, covid)~AstraZeneca, data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
astra_e1_p<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
astra_e3_p<-efficacy(m1)

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid)~AstraZeneca, data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
astra_e1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
astra_e3_c<-efficacy(m1)

### Mortalidad INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time_death),
                                       txtime = time_to_vaccine,
                                       fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, muerte)~AstraZeneca, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
astra_d1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, muerte)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
astra_d3_p<-efficacy(m1)

### Mortalidad COMPLETA ###
m1<-coxph(Surv(tstart, tstop, muerte)~AstraZeneca, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
astra_d1_c<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
astra_d3_c<-efficacy(m1)

### Hospitalización INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, hosp_covid)~AstraZeneca, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
astra_h1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, hosp_covid)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
astra_h3_p<-efficacy(m1)

### Hospitalización COMPLETA ###
m1<-coxph(Surv(tstart, tstop, hosp_covid)~AstraZeneca, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
astra_h1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~AstraZeneca+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
astra_h3_c<-efficacy(m1)

#### Vaccine effectiveness Pfizer ####
base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")

tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                     futime= pmax(.5, time),
                                     txtime = time_to_vaccine,
                                     fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))
### Esquema incompleto ###
m1<-coxph(Surv(tstart, tstop, covid)~Pfizer, data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
pfizer_e1_p<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
pfizer_e3_p<-efficacy(m1)

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid)~Pfizer, data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
pfizer_e1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
pfizer_e3_c<-efficacy(m1)

### Mortalidad INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                     futime= pmax(.5, time_death),
                                     txtime = time_to_vaccine,
                                     fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, muerte)~Pfizer, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
pfizer_d1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, muerte)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
pfizer_d3_p<-efficacy(m1)

### Mortalidad COMPLETA ###
m1<-coxph(Surv(tstart, tstop, muerte)~Pfizer, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
pfizer_d1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
pfizer_d3_c<-efficacy(m1)

### Hospitalización INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, hosp_covid)~Pfizer, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
pfizer_h1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, hosp_covid)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
pfizer_h3_p<-efficacy(m1)

### Hospitalización COMPLETA ###
m1<-coxph(Surv(tstart, tstop, hosp_covid)~Pfizer, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
pfizer_h1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~Pfizer+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
pfizer_h3_c<-efficacy(m1)

#### Vaccine effectiveness J&J ####
base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")

tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time),
                                       txtime = time_to_vaccine,
                                       fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ,data= sdata)
jj_e1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+EDAD+SEXO+comorb+strata(id),data= sdata)
jj_e3<-efficacy(m1)

### Mortalidad ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time_death),
                                       txtime = time_to_vaccine,
                                       fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(tstart, tstop, muerte)~JJ, data=sdata%>% filter(covid==1))
jj_d1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~JJ+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1))
jj_d3<-efficacy(m1)

### Hospitalization ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))

sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(tstart, tstop, hosp_covid)~JJ, data=sdata%>% filter(covid==1))
jj_h1<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~JJ+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1))
jj_h3<-efficacy(m1)

#### Vaccine effectiveness Sinovac ####
base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")

tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))
### Esquema incompleto ###
m1<-coxph(Surv(tstart, tstop, covid)~Sinovac, data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
sinovac_e1_p<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
sinovac_e3_p<-efficacy(m1)

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid)~Sinovac, data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
sinovac_e1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
sinovac_e3_c<-efficacy(m1)

### Mortalidad INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, muerte)~Sinovac, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sinovac_d1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, muerte)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sinovac_d3_p<-efficacy(m1)

### Mortalidad COMPLETA ###
m1<-coxph(Surv(tstart, tstop, muerte)~Sinovac, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sinovac_d1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sinovac_d3_c<-efficacy(m1)

### Hospitalización INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, hosp_covid)~Sinovac, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sinovac_h1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, hosp_covid)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
sinovac_h3_p<-efficacy(m1)

### Hospitalización COMPLETA ###
m1<-coxph(Surv(tstart, tstop, hosp_covid)~Sinovac, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sinovac_h1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~Sinovac+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
sinovac_h3_c<-efficacy(m1)

#### Vaccine effectiveness Moderna ####
base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")

tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time),
                                       txtime = time_to_vaccine,
                                       fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))
### Esquema incompleto ###
m1<-coxph(Surv(tstart, tstop, covid)~Moderna, data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
moderna_e1_p<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="COMPLETA"))
moderna_e3_p<-efficacy(m1)

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid)~Moderna, data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
moderna_e1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, covid)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata %>% filter(VACUNA_COV!="INCOMPLETA"))
moderna_e3_c<-efficacy(m1)

### Mortalidad INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                       futime= pmax(.5, time_death),
                                       txtime = time_to_vaccine,
                                       fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, muerte)~Moderna, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
moderna_d1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, muerte)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
moderna_d3_p<-efficacy(m1)

### Mortalidad COMPLETA ###
m1<-coxph(Surv(tstart, tstop, muerte)~Moderna, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
moderna_d1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, muerte)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
moderna_d3_c<-efficacy(m1)

### Hospitalización INCOMPLETA  ###
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = hosp_covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                hosp_covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

m1<-coxph(Surv(time_death, hosp_covid)~Moderna, data=sdata %>% filter(covid==1, VACUNA_COV!="COMPLETA"))
moderna_h1_p<-efficacy(m1)

m1<-coxph(Surv(time_death, hosp_covid)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="COMPLETA"))
moderna_h3_p<-efficacy(m1)

### Hospitalización COMPLETA ###
m1<-coxph(Surv(tstart, tstop, hosp_covid)~Moderna, data=sdata %>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
moderna_h1_c<-efficacy(m1)

m1<-coxph(Surv(tstart, tstop, hosp_covid)~Moderna+EDAD+SEXO+comorb+strata(id), data=sdata%>% filter(covid==1, VACUNA_COV!="INCOMPLETA"))
moderna_h3_c<-efficacy(m1)

#### Table 1: Vaccine effectiveness for all evaluated outcomes ####
t1<-as.numeric(table(base_vacuna$TIPO_VAC_COV))
table1<-data.frame(Vaccine=c("mRNA BNT162b2", "mRNA BNT162b2", "ChAdOx1", "ChAdOx1", "Gam-COVID-Vac", "Gam-COVID-Vac", 
                             "CoronaVac", "CoronaVac", "mRNA-1273", "mRNA-1273", "Ad5-nCoV","Ad26.COV2.S"),
                   Status=c("Partial", "Complete","Partial", "Complete","Partial", "Complete","Partial", "Complete", "Partial", "Complete","Single dose", "Single dose"),
                   Unadjusted=c(pfizer_e1_p, pfizer_e1_c, astra_e1_p, astra_e1_c, sputnik_e1_p, sputnik_e1_c,
                                sinovac_e1_p, sinovac_e1_c, moderna_e1_p, moderna_e1_c, cansino_e1, jj_e1),
                   Adj2=c(pfizer_e3_p, pfizer_e3_c, astra_e3_p, astra_e3_c, sputnik_e3_p, sputnik_e3_c,
                          sinovac_e3_p, sinovac_e3_c, moderna_e3_p, moderna_e3_c, cansino_e3, jj_e3))

names(table1)<-c("Vaccine", "Status", "Unadjusted efficacy", "Adjusted Efficacy")
tab1 <-align(flextable(table1,cwidth = c(1.5,1,2,2)),align = "center",part = "all")
save_as_docx(tab1,path="tabla1.docx")

table2<-data.frame(Vaccine=c("mRNA BNT162b2", "mRNA BNT162b2", "ChAdOx1", "ChAdOx1", "Gam-COVID-Vac", "Gam-COVID-Vac", 
                             "CoronaVac", "CoronaVac", "mRNA-1273", "mRNA-1273", "Ad5-nCoV","Ad26.COV2.S"),
                   Status=c("Partial", "Complete","Partial", "Complete","Partial", "Complete","Partial", "Complete", "Partial", "Complete","Single dose", "Single dose"),
                   Unadjusted=c(pfizer_d1_p, pfizer_d1_c, astra_d1_p, astra_d1_c, sputnik_d1_p, sputnik_d1_c,
                                sinovac_d1_p, sinovac_d1_c, moderna_d1_p, moderna_d1_c, cansino_d1, jj_d1),
                   Adj2=c(pfizer_d3_p, pfizer_d3_c, astra_d3_p, astra_d3_c, sputnik_d3_p, sputnik_d3_c,
                          sinovac_d3_p, sinovac_d3_c, moderna_d3_p, moderna_d3_c, cansino_d3, jj_d3))
names(table2)<-c("Vaccine","Status", "Unadjusted efficacy", "Adjusted Efficacy")
tab2 <-align(flextable(table2,cwidth = c(1.5,1,2,2)),align = "center",part = "all")
save_as_docx(tab2,path="tabla2.docx")

table3<-data.frame(Vaccine=c("mRNA BNT162b2", "mRNA BNT162b2", "ChAdOx1", "ChAdOx1", "Gam-COVID-Vac", "Gam-COVID-Vac", 
                             "CoronaVac", "CoronaVac", "mRNA-1273", "mRNA-1273", "Ad5-nCoV","Ad26.COV2.S"),
                   Status=c("Partial", "Complete","Partial", "Complete","Partial", "Complete","Partial", "Complete", "Partial", "Complete","Single dose", "Single dose"),
                   Unadjusted=c(pfizer_h1_p, pfizer_h1_c, astra_h1_p, astra_h1_c, sputnik_h1_p, sputnik_h1_c,
                                sinovac_h1_p, sinovac_h1_c, moderna_h1_p, moderna_h1_c, cansino_h1, jj_h1),
                   Adj2=c(pfizer_h3_p, pfizer_h3_c, astra_h3_p, astra_h3_c, sputnik_h3_p, sputnik_h3_c,
                          sinovac_h3_p, sinovac_h3_c, moderna_h3_p, moderna_h3_c, cansino_h3, jj_h3))

names(table3)<-c("Vaccine","Status", "Unadjusted efficacy", "Adjusted Efficacy")
tab3 <-align(flextable(table3,cwidth = c(1.5,1,2,2)),align = "center",part = "all")
save_as_docx(tab3,path="tabla3.docx")


#### Figure 2: Cumulative incidence graphs for all evaluated outcomes ####
nejm<-pal_nejm("default")(8)
mod2_kma<-survfit(Surv(time, covid)~TIPO_VAC_COV, data=base_final)
f3a<-ggsurvplot(mod2_kma, fun = "event", data = base_final, size = 1,palette = nejm,cconf.int = T, risk.table = F,pval = F, xlab="Days since December 24th, 2020",
                  ylab="Cumulative incidence (%)", title="SARS-CoV-2 infection", 
                legend.labs =c("ChAdOx1", "Ad5-nCoV","Gam-COVID-Vac","Ad26.COV2.S", "mRNA-1273","Unvaccinated","BNT162b2","CoronaVac"),
                ylim= c(0,0.6), xlim=c(0, 230), break.y.by= c(0.1), break.x.by= c(45), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ) )


## COVID-19 hospitalizations ###
fit <- survfit(Surv(time_death, hosp_covid)~TIPO_VAC_COV, data=base_final %>% filter(covid==1))
f3b<-ggsurvplot(fit, fun = "event", data = base_final, size = 1,palette = nejm,cconf.int = T, risk.table = F,pval = TRUE, xlab="Days since December 24th, 2020",
                ylab="Cumulative incidence (%)", title="COVID-19 related hospitalizations", pval.coord = c(0, 0.8),
                legend.labs =c("ChAdOx1", "Ad5-nCoV","Gam-COVID-Vac","Ad26.COV2.S", "mRNA-1273","Unvaccinated","BNT162b2","CoronaVac"),
                ylim= c(0,0.2), xlim=c(0, 230), break.y.by= c(0.04), break.x.by= c(45), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ) )

## COVID-19 deaths ###

fit1 <- survfit(Surv(time_death, muerte)~TIPO_VAC_COV, data=base_final %>% filter(covid==1))
f3c<-ggsurvplot(fit1, fun = "event", data = base_final, size = 1,palette = nejm,cconf.int = T, risk.table = F,pval = TRUE, xlab="Days since December 24th, 2020",
                ylab="Cumulative incidence (%)", title="COVID-19 related death", pval.coord = c(0, 0.8),
                legend.labs =c("ChAdOx1", "Ad5-nCoV","Gam-COVID-Vac","Ad26.COV2.S", "mRNA-1273","Unvaccinated","BNT162b2","CoronaVac"),
                ylim= c(0,0.2), xlim=c(0, 230), break.y.by= c(0.04), break.x.by= c(45), ggtheme =  ggpubr::theme_pubclean() + 
                  (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), text = element_text(hjust=0.5, family = "sans", size=15),
                         legend.text = element_text(hjust=0.5, size=12)) ) )

fig3<-ggarrange(f3a$plot, f3b$plot, f3c$plot, labels = c("A", "B", "C"), nrow=1, ncol=3)

ggsave(fig3,filename = "Figure2.jpg", 
       width = 55, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Delta vs. B.1.1.519 ####
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"))
cansino_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"))
cansino_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01" & VACUNA_COV!="INCOMPLETA"))
sputnik_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sputnik_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
astra_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
astra_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
pfizer_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
pfizer_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"))
jj_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"))
jj_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sinovac_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sinovac_delta<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
moderna_519<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
moderna_delta<-efficacy(m1)

vaccines<-c("ChadOx1", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_519), extract(cansino_519), extract(pfizer_519),
            extract(sputnik_519),extract(sinovac_519),extract(moderna_519), extract(jj_519))
delta<-rbind(extract(astra_delta), extract(cansino_delta), extract(pfizer_delta),extract(sputnik_delta),
             extract(sinovac_delta),extract(moderna_delta), extract(jj_delta))

variant<-rbind(b519, delta)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("B.1.1.519",7),rep("Delta",7))

v1<-variant %>%
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-3.2),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,110)+
  ylab("Effectiveness for symptomatic COVID-19")+labs(fill="Predominant SARS-CoV-2 Variant")+xlab("")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

### Death ###
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"))
cansino_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"))
cansino_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01" & VACUNA_COV!="INCOMPLETA"))
sputnik_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 &FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sputnik_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
astra_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
astra_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
pfizer_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
pfizer_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"))
jj_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"))
jj_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sinovac_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
sinovac_delta_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI<"2021-07-01"& VACUNA_COV!="INCOMPLETA"))
moderna_519_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+EDAD+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & FECINISI>="2021-07-01"& VACUNA_COV!="INCOMPLETA"))
moderna_delta_d<-efficacy(m1)

vaccines<-c("ChadOx1", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_519_d), extract(cansino_519_d), extract(pfizer_519_d),
            extract(sputnik_519_d),extract(sinovac_519_d),extract(moderna_519_d), extract(jj_519_d))
delta<-rbind(extract(astra_delta_d), extract(cansino_delta_d), extract(pfizer_delta_d),extract(sputnik_delta_d),
             extract(sinovac_delta_d),extract(moderna_delta_d), extract(jj_delta_d))

variant<-rbind(b519, delta)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("B.1.1.519",7),rep("Delta",7))

variant<-rbind(b519, delta)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("B.1.1.519",7),rep("Delta",7))
variant$le[7]<-14.5
variant$le[6]<-65.0
variant$ue[6]<-100.0

d1<-variant %>% 
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2.4),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,110)+
  ylab("Effectiveness for COVID-19 death")+labs(fill="Predominant SARS-CoV-2 Variant")+xlab("")+
  theme(axis.ticks = element_blank())

#### Age comparison ####
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60))
cansino_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60))
cansino_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60 & VACUNA_COV!="INCOMPLETA"))
sputnik_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60& VACUNA_COV!="INCOMPLETA"))
sputnik_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60& VACUNA_COV!="INCOMPLETA"))
astra_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60& VACUNA_COV!="INCOMPLETA"))
astra_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60& VACUNA_COV!="INCOMPLETA"))
pfizer_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60& VACUNA_COV!="INCOMPLETA"))
pfizer_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60))
jj_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60))
jj_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60& VACUNA_COV!="INCOMPLETA"))
sinovac_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60& VACUNA_COV!="INCOMPLETA"))
sinovac_u65<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD>=60& VACUNA_COV!="INCOMPLETA"))
moderna_o65<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(EDAD<60& VACUNA_COV!="INCOMPLETA"))
moderna_u65<-efficacy(m1)

vaccines<-c("ChadOx1", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_o65), extract(cansino_o65), extract(pfizer_o65),
            extract(sputnik_o65),extract(sinovac_o65),extract(moderna_o65), extract(jj_o65))
delta<-rbind(extract(astra_u65), extract(cansino_u65), extract(pfizer_u65),extract(sputnik_u65),
             extract(sinovac_u65),extract(moderna_u65), extract(jj_u65))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("Under 60 years",7),rep("60 years or older",7))
variant<-variant %>%mutate(across(variant, factor, levels=c("Under 60 years","60 years or older")))

v2<-variant %>%
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,105)+
  ylab("Effectiveness for symptomatic COVID-19")+labs(fill="Age group")+xlab("")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

### Death ###
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60))
cansino_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60))
cansino_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60 & VACUNA_COV!="INCOMPLETA"))
sputnik_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60 & VACUNA_COV!="INCOMPLETA"))
sputnik_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60& VACUNA_COV!="INCOMPLETA"))
astra_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60& VACUNA_COV!="INCOMPLETA"))
astra_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60& VACUNA_COV!="INCOMPLETA"))
pfizer_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60& VACUNA_COV!="INCOMPLETA"))
pfizer_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60))
jj_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60))
jj_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60& VACUNA_COV!="INCOMPLETA"))
sinovac_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60& VACUNA_COV!="INCOMPLETA"))
sinovac_u65_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD>=60& VACUNA_COV!="INCOMPLETA"))
moderna_o65_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & EDAD<60& VACUNA_COV!="INCOMPLETA"))
moderna_u65_d<-efficacy(m1)

vaccines<-c("ChadOx", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_o65_d), extract(cansino_o65_d), extract(pfizer_o65_d),
            extract(sputnik_o65_d),extract(sinovac_o65_d),extract(moderna_o65_d), extract(jj_o65_d))
delta<-rbind(extract(astra_u65_d), extract(cansino_u65_d), extract(pfizer_u65_d),extract(sputnik_u65_d),
             extract(sinovac_u65_d),extract(moderna_u65_d), extract(jj_u65_d))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("Under 60 years",7),rep("60 years or older",7))
variant<-variant %>%mutate(across(variant, factor, levels=c("Under 60 years","60 years or older")))

d2<-variant %>%
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2.4),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,110)+
  ylab("Effectiveness for COVID-19 death")+labs(fill="Age group")+xlab("")+
  theme(axis.ticks = element_blank())

#### Effectiveness in diabetes ####
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1))
cansino_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1))
cansino_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
sputnik_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
sputnik_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
astra_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
astra_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
pfizer_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
pfizer_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1))
jj_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1))
jj_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
sinovac_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
sinovac_nodiab<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
moderna_diab<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
moderna_nodiab<-efficacy(m1)

### Build Figure ###
vaccines<-c("ChadOx1", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_diab), extract(cansino_diab), extract(pfizer_diab),
            extract(sputnik_diab),extract(sinovac_diab),extract(moderna_diab), extract(jj_diab))
delta<-rbind(extract(astra_nodiab), extract(cansino_nodiab), extract(pfizer_nodiab),extract(sputnik_nodiab),
             extract(sinovac_nodiab),extract(moderna_nodiab), extract(jj_nodiab))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("No diabetes",7),rep("Diabetes",7))
variant<-variant %>%mutate(across(variant, factor, levels=c("No diabetes", "Diabetes")))
  

v3<-variant %>%
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,105)+
  ylab("Effectiveness for symptomatic COVID-19")+labs(fill="Status")+xlab("")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

### Death ###
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1))
cansino_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1))
cansino_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1 & VACUNA_COV!="INCOMPLETA"))
sputnik_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1 & VACUNA_COV!="INCOMPLETA"))
sputnik_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1& VACUNA_COV!="INCOMPLETA"))
astra_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1& VACUNA_COV!="INCOMPLETA"))
astra_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1& VACUNA_COV!="INCOMPLETA"))
pfizer_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1& VACUNA_COV!="INCOMPLETA"))
pfizer_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1))
jj_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1))
jj_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1& VACUNA_COV!="INCOMPLETA"))
sinovac_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1& VACUNA_COV!="INCOMPLETA"))
sinovac_nodiab_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES==1& VACUNA_COV!="INCOMPLETA"))
moderna_diab_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & DIABETES!=1& VACUNA_COV!="INCOMPLETA"))
moderna_nodiab_d<-efficacy(m1)

vaccines<-c("ChadOx", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_diab_d), extract(cansino_diab_d), extract(pfizer_diab_d),
            extract(sputnik_diab_d),extract(sinovac_diab_d),extract(moderna_diab_d), extract(jj_diab_d))
delta<-rbind(extract(astra_nodiab_d), extract(cansino_nodiab_d), extract(pfizer_nodiab_d),extract(sputnik_nodiab_d),
             extract(sinovac_nodiab_d),extract(moderna_nodiab_d), extract(jj_nodiab_d))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("No diabetes",7),rep("Diabetes",7))
variant<-variant %>%mutate(across(variant, factor, levels=c("No diabetes", "Diabetes")))

d3<-variant %>% 
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2.5),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,105)+
  ylab("Effectiveness for COVID-19 death")+labs(fill="Status")+xlab("")+
  theme(axis.ticks = element_blank())

#### Effectiveness in obesity ####
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1))
cansino_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1))
cansino_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
sputnik_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
sputnik_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
astra_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
astra_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
pfizer_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
pfizer_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1))
jj_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1))
jj_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
sinovac_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
sinovac_noob<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time),
                                      txtime = time_to_vaccine,
                                      fustat = covid))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                covid = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
moderna_ob<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, covid) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
moderna_noob<-efficacy(m1)

### Build Figure ###
vaccines<-c("ChadOx", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_ob), extract(cansino_ob), extract(pfizer_ob),
            extract(sputnik_ob),extract(sinovac_ob),extract(moderna_ob), extract(jj_ob))
delta<-rbind(extract(astra_noob), extract(cansino_noob), extract(pfizer_noob),extract(sputnik_noob),
             extract(sinovac_noob),extract(moderna_noob), extract(jj_noob))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("No obesity",7),rep("Obesity",7))
variant<-variant %>%mutate(across(variant, factor, levels=c("No obesity", "Obesity")))

v4<-variant %>%
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,105)+
  ylab("Effectiveness for symptomatic COVID-19")+labs(fill="Status")+xlab("")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

### Death ###
base_vacuna<-base_final %>% filter(CanSino==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                CaSino = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1))
cansino_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ CanSino+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1))
cansino_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sputnik==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sputnik = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1 & VACUNA_COV!="INCOMPLETA"))
sputnik_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sputnik+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1 & VACUNA_COV!="INCOMPLETA"))
sputnik_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(AstraZeneca==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                AstraZeneca = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1& VACUNA_COV!="INCOMPLETA"))
astra_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ AstraZeneca+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1& VACUNA_COV!="INCOMPLETA"))
astra_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Pfizer==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Pfizer = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1& VACUNA_COV!="INCOMPLETA"))
pfizer_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Pfizer+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1& VACUNA_COV!="INCOMPLETA"))
pfizer_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(JJ==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                JJ = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1))
jj_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ JJ+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1))
jj_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Sinovac==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Sinovac = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1& VACUNA_COV!="INCOMPLETA"))
sinovac_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Sinovac+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1& VACUNA_COV!="INCOMPLETA"))
sinovac_noob_d<-efficacy(m1)

base_vacuna<-base_final %>% filter(Moderna==1 | VACUNA_COV=="No vacunado") %>% filter(EDAD>=18) %>% filter(days_from_vaccine>=14 | VACUNA_COV=="No vacunado")
tdata <- with(base_vacuna, data.frame(ID_REGISTRO = ID_REGISTRO,
                                      futime= pmax(.5, time_death),
                                      txtime = time_to_vaccine,
                                      fustat = muerte))
sdata <- tmerge(base_vacuna, tdata, id=ID_REGISTRO,
                muerte = event(futime, fustat),
                Moderna = tdc(txtime),
                options= list(idname="subject"))

### Esquema completo ###
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD==1& VACUNA_COV!="INCOMPLETA"))
moderna_ob_d<-efficacy(m1)
m1<-coxph(Surv(tstart, tstop, muerte) ~ Moderna+SEXO+comorb+strata(id),data= sdata %>% filter(covid==1 & OBESIDAD!=1& VACUNA_COV!="INCOMPLETA"))
moderna_noob_d<-efficacy(m1)

vaccines<-c("ChadOx", "Ad5-nCoV", "BNT162b2","Gam-COVID-Vac","CoronaVac","mRNA-1273","Ad26.COV2.S")
b519<-rbind(extract(astra_ob_d), extract(cansino_ob_d), extract(pfizer_ob_d),
            extract(sputnik_ob_d),extract(sinovac_ob_d),extract(moderna_ob_d), extract(jj_ob_d))
delta<-rbind(extract(astra_noob_d), extract(cansino_noob_d), extract(pfizer_noob_d),extract(sputnik_noob_d),
             extract(sinovac_noob_d),extract(moderna_noob_d), extract(jj_noob_d))

variant<-rbind(delta,b519)
variant$vaccine<-c(vaccines, vaccines)
variant$variant<-c(rep("No obesity",7),rep("Obesity",7))

d4<-variant %>% 
  ggplot(aes(x=variant, y=ef, fill=variant)) + 
  geom_bar(stat="identity", position = "dodge", color="black", width=0.7, alpha=0.8) + 
  geom_errorbar(aes(ymin=le, ymax=ue), width=0.3) +
  geom_text(aes(label=ef, vjust=-2.0),col="black", size=4)+
  facet_wrap(~vaccine, ncol=7, nrow=1)+scale_fill_jama()+theme_pubclean()+ylim(0,110)+
  ylab("Effectiveness for COVID-19 death")+labs(fill="Status")+xlab("")+
  theme(axis.ticks = element_blank())

#### Figure 3: Stratified vaccine effectiveness per subgroups ####

fig3<-ggarrange(v1,v2,v3, nrow=3, ncol=1, labels = c("A", "B", "C"))

ggsave(fig3,filename = "Figure3.jpg", 
       width = 35, 
       height = 45,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

ggsave(v4,filename = "SuppFigure3.jpg", 
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

figs1<-ggarrange(d1,d2,d3,d4, nrow=4, ncol=1, labels = c("A", "B", "C", "D"))

ggsave(figs1,filename = "SuppFigure1.jpg", 
       width = 35, 
       height = 55,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

