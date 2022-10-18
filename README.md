COVID-SEIR-India
----------------
This repository stores the code, equations, and parameters associated with the Michael Group SEIR model.

# Project Guide
## Fitting and running the base scenario
The repository comes with India's daily case, death, and vaccination data up to
May 5th, which can be found in `India.csv`. The main project is run via the
script `Main.m`:

`matlab -nodisplay -nosplash < Main.m`

This will produce output `India.mat`, which contains the predictions of all
state functions until the end of the year. With this file loaded, you can plot
median proportion of people immune using the following code:

`plot(median(V+B+R2, 2));`

Other state functions can be visualized similarly.

## Running alternative social measure / vaccination scenarios

The impact of social measures on transmission is captured via a scaling factor,
`d`. To simulate 30 days of +25% increased social measures, add the following
lines to `diff_eqn1.m`, just before the definitions of the differential
equations:
```
if t > 427 && t < 427+30
     d = d*1.25;
end
 ```

where t = 427 is May 5th, the final data point.

To increase/decrease the vaccination rate, adjust line 444 in
`BM_SEIR_model.m`. For instance, to double the vaccination rate
going forward:

`totalv = 2*mean(Vaccinated(end-21:end));`

This would apply 2 times the average daily vaccination rate over the
last 21 days.

# Extended SEIR Model 
In this study, we simulated the ongoing SARS-CoV-2 outbreaks in India using a variation of an SEIR model described in detail in Newcomb and colleagues (1). The Ordinary Differential Equations (ODEs) describing the model are given fully below in Equations. Briefly, we assume each country is a closed population and ignore demographic changes such that the total population size remains constant. The population is divided into compartments representing various infection stages: susceptible (S), exposed (E), infectious asymptomatic (IA), infectious pre-symptomatic (IP), infectious with mild symptoms (IM), infectious with severe symptoms requiring hospitalization (IH), infectious with severe symptoms requiring intensive care including ventilation (IC), recovered and immune (R), first-dose vaccinated (V), completely vaccinated (B), and deceased (D). 

The specific transitions and rate parameters governing the evolution of the system, along with their prior and posterior fitted values, are described in the Table. The strength of social distancing measures as a result of public health policies to limit contacts is captured through the estimation of a scaling factor, d, which is in turn multiplied by the transmission rate, beta, to obtain the population-level transmission intensity operational at any given time in each population. This factor accounts for the transmission modifying effects of mask wearing, reductions in mobility and mixing, work from home, and any other deviations from the normal social behavior of each population prior to the epidemic. The vaccination data for India is directly applied by moving the proportion of the population that is vaccinated over a 10- day block from the S class to the V (1st dose) class. Individuals then move from the V to the B (2nd dose booster) class at a daily rate approximating a 21 day interval between vaccine doses. Average vaccination rates estimated from the last 3 weeks of vaccination data (April 15th -May 6th) were used to simulate into the future. The future impacts of changes in social mitigation interventions and vaccination rates are simulated by altering the values of d and the vaccination rate.

System of ODEs
-------------

![System of Equations](equations.png)

Table of Parameters/Priors
---------------------
Model parameter priors, along with best-fitting values.

| Parameter | Definition | **Prior range** | **Median Fit, India** | Units/notes |
| --- | --- | --- | --- | --- | --- |
| β | Infection transmission rate | **0.125 – 2.0** | **0.3254** | Estimated as R0\*gamma in SIR model |
| σ | Rate of moving from exposed class to infectious class | **0.16 – 0.5** | **0.3016** | 1/σ is the latent period; assumed 2-6 days |
| ⍴ | Proportion of exposed who become asymptomatic | **0.25 – 0.50** | **0.3785** |   |
| γA | Recovery rate of asymptomatic cases | **0.125 – 0.33** | **0.2314** | 1/γA is the infectious period; assumed 3-8 days |
| γM | Recovery rate of cases with mild symptoms | **0.125 – 0.33** | **0.2324** | 1/γM is the infectious period; assumed 3-8 days |
| γH | Recovery rate of cases with severe symptoms requiring hospitalization | **0.125 – 0.33** | **0.2203** | 1/γH is the infectious period of severe cases; assumed 3-8 days |
| γC | Recovery rate of cases with severe symptoms requiring intensive care | **0.125 – 0.33** | **0.2293** | 1/γC is the infectious period; assumed 3-8 days |
| δ1 | Rate of moving from presymptomatic class to mild symptomatic | **0.05 – 0.20** | **0.1600** | 1/time from start of infectious period to illness onset; assume 5-20 days |
| δ2 | Rate of moving from mild case to hospitalized class | **0.06 – 0.25** | **0.1474** | 1/time from illness onset to hospitalization; assume 4-15 days |
| δ3 | Rate of moving from hospitalized class to ICU | **0.09 – 1** | **0.4921** | 1/time from hospitalization to ICU; assume 1-11 days |
| m | Mortality rate of ICU class | **0.08 – 0.25** | **0.1519** | 1/time from ICU to death |
 | Proportion of cases detected by testing | **0.1 – 0.3** | **0.2054** |
| x1 | Proportion of mild cases that progress to hospital | **0.05 – 0.3** | **0.1573** | 5-30% of mild cases are hospitalized |
| x2 | Proportion of hospital cases that progress to ICU | **0.2 – 0.3** | **0.2497** | 20-30% of hospitalized cases require an ICU |
| x3 | Proportion of ICU cases that die | **0.2 – 0.8** | **0.4206** | Proportion of ICU cases that die |
| d | Reduction in transmission due to social distancing, face masks, etc. | **0.25 – 0.9** | **0.4626**  |
| 𝜀v | Vaccine Efficacy | **Fixed, 0.90** |   |
| 𝜀B | Booster Efficacy | **Fixed, 0.75** |   |
| ξv | Vaccination Rate | **Varies over time, according to vaccination data** |
| ξB | Booster Rate | **Fixed, 0.025** |
