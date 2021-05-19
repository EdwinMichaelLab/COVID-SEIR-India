COVID-SEIR-India
----------------
This repository stores the code, equations, and parameters associated with the Michael Group SEIR model.


System of ODEs
-------------

![System of Equations](equations.png)

Table of Parameters/Priors
---------------------
Model parameter priors, along with best-fitting values.

| Parameter | Definition | **Prior range** | **Median Fit, USA** | **Median Fit, India** | Units/notes |
| --- | --- | --- | --- | --- | --- |
| β | Infection transmission rate | **0.125 – 2.0** | **0.4188** | **0.3254** | Estimated as R0\*gamma in SIR model |
| σ | Rate of moving from exposed class to infectious class | **0.16 – 0.5** | **0.3083** | **0.3016** | 1/σ is the latent period; assumed 2-6 days |
| ⍴ | Proportion of exposed who become asymptomatic | **0.25 – 0.50** | **0.3887** | **0.3785** |   |
| γA | Recovery rate of asymptomatic cases | **0.125 – 0.33** | **0.2353** | **0.2314** | 1/γA is the infectious period; assumed 3-8 days |
| γM | Recovery rate of cases with mild symptoms | **0.125 – 0.33** | **0.2423** | **0.2324** | 1/γM is the infectious period; assumed 3-8 days |
| γH | Recovery rate of cases with severe symptoms requiring hospitalization | **0.125 – 0.33** | **0.2372** | **0.2203** | 1/γH is the infectious period of severe cases; assumed 3-8 days |
| γC | Recovery rate of cases with severe symptoms requiring intensive care | **0.125 – 0.33** | **0.2372** | **0.2293** | 1/γC is the infectious period; assumed 3-8 days |
| δ1 | Rate of moving from presymptomatic class to mild symptomatic | **0.05 – 0.20** | **0.1605** | **0.1600** | 1/time from start of infectious period to illness onset; assume 5-20 days |
| δ2 | Rate of moving from mild case to hospitalized class | **0.06 – 0.25** | **0.1493** | **0.1474** | 1/time from illness onset to hospitalization; assume 4-15 days |
| δ3 | Rate of moving from hospitalized class to ICU | **0.09 – 1** | **0.5098** | **0.4921** | 1/time from hospitalization to ICU; assume 1-11 days |
| m | Mortality rate of ICU class | **0.08 – 0.25** | **0.1497** | **0.1519** | 1/time from ICU to death |
 | Proportion of cases detected by testing | **0.1 – 0.3** | **0.2100** | **0.2054** |
| x1 | Proportion of mild cases that progress to hospital | **0.05 – 0.3** | **0.1657** | **0.1573** | 5-30% of mild cases are hospitalized |
| x2 | Proportion of hospital cases that progress to ICU | **0.2 – 0.3** | **0.2477** | **0.2497** | 20-30% of hospitalized cases require an ICU |
| x3 | Proportion of ICU cases that die | **0.2 – 0.8** | **0.4545** | **0.4206** | Proportion of ICU cases that die |
| d | Reduction in transmission due to social distancing, face masks, etc. | **0.25 – 0.9** | **0.4753** | **0.4626** |   |
| 𝜀v | Vaccine Efficacy | **Fixed, 0.90** |   |
| 𝜀B | Booster Efficacy | **Fixed, 0.75** |   |
| ξv | Vaccination Rate | **Varies over time, according to vaccination data** |
| ξB | Booster Rate | **Fixed, 0.025** |
