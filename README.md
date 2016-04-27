# Prj3-ZikaVirus
## CX4230 - Computer Simulation - Project 3 - Zika Virus Propagation

### By: Derrick Williams and Tilak Patel

#### Date: April 27, 2016

### Zika Virus Propagation Task:
Our task is to understand the potential propagation dynamics of the Zika virus via airline route network and a Reproduction Number of 4.

###Conceptual Model:
We used the Kermack-McKendrick SIR model with delay, incubation, and O’Leary’s add-on for vaccination.  We extended several existing models on domestic mosquito abundance and major airline route network to include the SIR model.

### How to Run zikaSim.py:
- The basic format if you are running from a command line prompt is "python zikaSim.py [-m] [-s] [-a] [-v] [--c (city to infect)] [--d (date to infect)] [--start (start date of simulation)] [--days (number of days to run simulation)] [-- tau (new tau for disease)] [--inc (new incubation days)] [--vac (new vaccination rate)] [--screen (new screening percent)] ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv
- If running in say pycharm, set edit configuration to "[-m] [-s] [-a] [-v] [--c (city to infect)] [--d (date to infect)] [--start (start date of simulation)] [--days (number of days to run simulation)] [-- tau (new tau for disease)] [--inc (new incubation days)] [--vac (new vaccination rate)] [--screen (new screening percent)] ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv
".
- "-m" is for producing a map of the network simulation for one city
- "-s" is for showing stats at starting infection city
- "-a" run every simulation with starting infection in every month and in every city
- "-v" vaccinate using default percent based on stopping infection with 1 - 1/TAU
- "--c (city to infect)" specify what city to infect
- "--d (date to infect)" specity the date the infection should start
- "--start (start date of simulation)" specify the start date of the simulation
- "--days (number of days to run simulation)" specify how many days to run the simulation after the start date
- "--tau (new tau for disease)" specify new tau for a different disease or assume less aggressive tau for the zika virus
- "--inc (new incubation days)" specify new incubination days for virus to start showing symptoms of disease and person can spread to others then
- "--vac (new vaccination rate)" specify new vaccination rate to work with O'Leary vaccination formula
- "--screen (new screening percent)" specify new screening percent to screen out passengers from airline travel
