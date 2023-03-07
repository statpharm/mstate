## https://www.jstatsoft.org/article/download/v038i07/444

library(mstate)
## Walkthrough the code from 

## mstate: An R Package for the Analysis of Competing Risks and Multi-State
#Models

# https://www.jstatsoft.org/article/download/v038i07/444

# multi-state analysis of data describing ...

## post-transplant events of patients

# with blood cancer. The data have been provided by the EBMT (the European Group
# for Blood and Marrow Transplantation).


## Data, questions, and model
#
# Three intermediate events are included in the model: Recovery (Rec), an
# Adverse Event (AE) and a combination of the two (AE and Rec). It is to be
# expected that recovery improves the prognosis and an adverse event
# deteriorates it.

# The model is suitable to show the size of these # effects, and to capture the
# influence of their timing and of the covariates on the prognosis.

#
# Four prognostic factors are known at baseline for all patients (see Table 1).
# They are: donor-recipient match (where gender mismatch is defined as female
# donor, male recipient), prophylaxis, year of transplant and age at transplant
# in years. All these covariates are treated as time-fixed categorical
# covariates.
#
# A multi-state approach is particularly appropriate for these data, since it
# can help to model both the disease-related and the treatment-related morbidity
# and mortality. These are mod- eled here by including the intermediate events
# recovery and adverse events

# We consider the following six-states model (see Figure 1):
# 
# 1. Alive and in remission, no recovery or adverse event;
# 2. Alive in remission, recovered from the treatment;
# 3. Alive in remission, occurrence of the adverse event;
# 4. Alive, both recovered and adverse event occurred;
# 5. Alive, in relapse (treatment failure);
# 6. Dead (treatment failure).
data("ebmt4")
ebmt <- ebmt4
head(ebmt)
# For instance, patient 1 had recovered after 22 days (transition from state 1
# to state 2) and was censored after 995 days without a further event. Patient 2
# experienced the adverse event after 12 days (transition from state 1 to state
# 3), then recovery after 29 days (transition from state 3 to state 4) and a
# relapse after 422 days (transition from state 4 to state 5). Finally, he/she
# died after 579 days, but this last event is not relevant to the model, because 
# the patient had already reached an absorbing state.

# Diverse clinical questions can be answered by fitting this model to the data, such as: 
# (1) how does age influence the different phases of the disease/recovery process? 
# (2) How does the occurrence of the adverse event after 2 months change the prognosis of 10-year survival for a patient?
# (3) Which risks should be monitored most carefully for a patient who has recovered after one month, taking into account certain covariate values?
tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                          c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", 
                                               "Rel", "Death"))
tmat
msebmt <- msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae",
                                                     "recae", "rel", "srv"), 
                 status = c(NA, "rec.s", "ae.s", "recae.s", 
                            "rel.s", "srv.s"), 
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 1, c(1:8, 10:12)]
events(msebmt)

covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
msebmt[msebmt$id == 1, -c(9, 10, 12:48, 61:84)]

msebmt[, c("Tstart", "Tstop", "time")] <- msebmt[, c("Tstart",
                                                     "Tstop", "time")]/365.25
c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msebmt,
           method = "breslow")
summary(c0)
c0
msf0 <- msfit(object = c0, vartype = "greenwood", trans = tmat)
msf0
summary(msf0)

plot(msf0, las = 1, lty = rep(1:2, c(8, 4)),
      xlab = "Years since transplantation")
msf0a <- msfit(object = c0, vartype = "aalen", trans = tmat)
head(msf0a$varHaz)


pt0 <- probtrans(msf0, predt = 0, method = "greenwood")
pt0a <- probtrans(msf0a, predt = 0, method = "aalen")
summary(pt0, from = 1)
plot(pt0)
plot(pt0a)


library("colorspace")
statecols <- heat_hcl(6, c = c(80, 30), l = c(30, 90),
  power = c(1/5, 2))[c(6, 5, 3, 4, 2, 1)]
 ord <- c(2, 4, 3, 5, 6, 1)
 plot(pt0, ord = ord, xlab = "Years since transplantation",
        las = 1, type = "filled", col = statecols[ord])

 
 
 
pt100 <- probtrans(msf0, predt = 100/365.25, method = "greenwood")
plot(pt100, ord = c(2, 4, 3, 5, 6, 1),
  xlab = "Years since transplantation", main = "Starting from transplant",
  xlim = c(0, 10), las = 1, type = "filled", col = statecols[ord])
plot(pt100, from = 3, ord = c(2, 4, 3, 5, 6, 1),
  xlab = "Years since transplantation", main = "Starting from AE",
  xlim = c(0, 10), las = 1, type = "filled", col = statecols[ord]) 



cfull <- coxph(Surv(Tstart, Tstop, status) ~ match.1 +
                  match.2 + match.3 + match.4 + match.5 + match.6 +
                  match.7 + match.8 + match.9 + match.10 + match.11 +
                  match.12 + proph.1 + proph.2 + proph.3 + proph.4 +
                  proph.5 + proph.6 + proph.7 + proph.8 + proph.9 +
                  proph.10 + proph.11 + proph.12 + year1.1 + year1.2 +
                  year1.3 + year1.4 + year1.5 + year1.6 + year1.7 +
                 year1.8 + year1.9 + year1.10 + year1.11 + year1.12 +
                 year2.1 + year2.2 + year2.3 + year2.4 + year2.5 +
                  year2.6 + year2.7 + year2.8 + year2.9 + year2.10 +
                  year2.11 + year2.12 + agecl1.1 + agecl1.2 + agecl1.3 +
                  agecl1.4 + agecl1.5 + agecl1.6 + agecl1.7 + agecl1.8 +
                  agecl1.9 + agecl1.10 + agecl1.11 + agecl1.12 + agecl2.1 +
                  agecl2.2 + agecl2.3 + agecl2.4 + agecl2.5 + agecl2.6 +
                  agecl2.7 + agecl2.8 + agecl2.9 + agecl2.10 + agecl2.11 +
                  agecl2.12 + strata(trans), data = msebmt, method = "breslow")


whA <- which(msebmt$proph == "yes" & msebmt$match == "no gender mismatch"
              & msebmt$year == "1995-1998" & msebmt$agecl == "<=20")
patA <- msebmt[rep(whA[1], 12), 9:12]
patA$trans <- 1:12
attr(patA, "trans") <- tmat
patA <- expand.covs(patA, covs, longnames = FALSE)
patA$strata <- patA$trans
msfA <- msfit(cfull, patA, trans = tmat)
ptA <- probtrans(msfA, predt = 0)

plot(ptA, ord = c(2, 4, 3, 5, 6, 1), main = "Patient A",
     las = 1, xlab = "Years since transplantation", xlim = c(0, 10),
     type = "filled", col = statecols[ord])
pt100A <- probtrans(msfA, predt = 100/365.25)

plot(pt100A, from = 3, ord = c(2, 4, 3, 5, 6, 1), main = "Patient A",
      las = 1, xlab = "Years since transplantation", xlim = c(0, 10),
      type = "filled", col = statecols[ord])
ptA10yrs <- probtrans(msfA, predt = 10, direction = "fixedhorizon")
head(ptA10yrs[[1]])
