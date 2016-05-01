# Prj3-ZikaVirus
## CX4230 - Computer Simulation - Project 3 - Zika Virus Propagation

### By: Derrick Williams and Tilak Patel

#### Date: April 27, 2016

### Zika Virus Propagation Task:
Our task is to understand the potential propagation dynamics of the Zika virus via airline route network and a Reproduction Number of 4.

###Conceptual Model:
We used the Kermack-McKendrick SIR model with delay, incubation, and O’Leary’s add-on for vaccination.  We extended several existing models on domestic mosquito abundance and major airline route network to include the SIR model.

### How to Run zikaSim.py:
- Main simulation file
- The basic format if you are running from a command line prompt is "python zikaSim.py [-m] [-s] [-a] [-v] [--c (city to infect)] [--d (date to infect)] [--start (start date of simulation)] [--days (number of days to run simulation)] [-- tau (new tau for disease)] [--inc (new incubation days)] [--vac (new vaccination rate)] [--screen (new screening percent)] ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv
- If running in say pycharm, set edit configuration to "[-m] [-s] [-a] [-v] [--c (city to infect)] [--d (date to infect)] [--start (start date of simulation)] [--days (number of days to run simulation)] [-- tau (new tau for disease)] [--inc (new incubation days)] [--vac (new vaccination rate)] [--screen (new screening percent)] ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv
".
- "-m" is for producing a map of the network
- "-s" is for showing stats at starting infection city
- "-a" run every simulation with starting infection in every month and in every city
- "-v" vaccinate using default percent based on stopping infection with 1 - 1/TAU
- "-z" show visualizations from one city for entire year
- "--c (city to infect)" specify what city to infect
- "--d (date to infect)" specity the date the infection should start
- "--start (start date of simulation)" specify the start date of the simulation
- "--days (number of days to run simulation)" specify how many days to run the simulation after the start date
- "--tau (new tau for disease)" specify new tau for a different disease or assume less aggressive tau for the zika virus
- "--inc (new incubation days)" specify new incubination days for virus to start showing symptoms of disease and person can spread to others then
- "--vac (new vaccination rate)" specify new vaccination rate to work with O'Leary vaccination formula
- "--screen (new screening percent)" specify new screening percent to screen out passengers from airline travel
- "--rec" change recovered time for humans

### How to Run convertAirports.py:
- Converts airports to appropriate data format for use in simulation file
- The basic format if you are running from a command line prompt is "python convertAirports.py ./Data/airports.dat"
- If running in say pycharm, set edit configuration to "./Data/airports.dat"

### How to Run convertAirlinePassenger.py:
- Converts passenger and airline data to appropriate data format for use in simulation file
- The basic format if you are running from a command line prompt is "python convertAirlinePassenger.py ./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv ./Data/airportsMin.csv"
- If running in say pycharm, set edit configuration to "./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv ./Data/airportsMin.csv"

### How to Run powerLawAirlineRoutes.py:
- Creates figure of flight data based on number of out-degrees vs count of the out-degree airports
- The basic format if you are running from a command line prompt is " python powerLawAirlineRoutes.py ./Data/airports.dat ./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv"
- If running in say pycharm, set edit configuration to "./Data/airports.dat ./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv"

### Data Files:
- mosCurves.csv - was created using Fig. 2 of this paper: 
Monaghan AJ, Morin CW, Steinhoff DF, Wilhelmi O, Hayden M, Quattrochi DA, Reiskind M, Lloyd AL, Smith K, Schmidt CA, Scalf PE, Ernst K. On the Seasonal Occurrence and Abundance of the Zika Virus Vector Mosquito Aedes Aegypti in the Contiguous United States. PLOS Currents Outbreaks. 2016 Mar 16 . Edition 1. doi: 10.1371/currents.outbreaks.50dfc7f46798675fc63e7d7da563da76.
- airportsMin.csv - was created using airports.dat file which came from http://openflights.org/data.html and running the file through convertAirports.py
- airlineRoutesPassengerData.csv - was created using travel data (182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv) from United States Department of Transportation - Bureau of Transportation Statistics for airline travel from all of 2015 and running the file through convertAirlinePassenger.py
- airports.dat - data from http://openflights.org/data.html
- 182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv - data from United States Department of Transportation - Bureau of Transportation Statistics
- airport_data.xlsx - data from wikipedia.org for determining which cities/hubs to consider for the simulation
- powerLawData.csv - output data from running powerLawAirlineRoutes.py; in format of [number of degrees, count]
