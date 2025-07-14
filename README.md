# superlearner_survival_tutorial

The super learner for time-to-event outcomes: A tutorial
Ruth Keogh, Karla Diaz-Ordaz, Jon Michael Gran, Nan van Geloven, possibly others

Code illustrating the use of the superlearner for time-to-event outcomes. Three methods are used: the discrete-time super learner of Polley and van der Laan (2011), the continuous-time super learner of Westling et al. (2023), and the continuous-time super learner of Munch and Gerds. (2024). The methods are applied to the rotterdam data set available in the 'survival' package in R. Methods are implemented using the packages 'superlearner', 'survSuperLearner' and 'statelearner'. By-hand implementations are also provided to clearly illustrate the steps involved. 

References

Munch, A. and Gerds, T.A. (2024). The state learner – a super learner for right-censored data. arXiv, arXiv:2405.17259.

Polley, E.C. and van der Laan, M.J. (2011). Super learning for right-censored data, 1st edition. New York: Springer., pp. 249–258.

Westling, T., Luedtke, A., Gilbert, P. B. and Carone, M. (2023). Inference for treatment-specific survival curves using machine learning. Journal of the American Statistical Association 119(546), 1541–1553.
