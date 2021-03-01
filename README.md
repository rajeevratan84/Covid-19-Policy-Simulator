# Covid-19-Policy-Simulator

A simple RShiny App that simulates changes in covid-19 infection rates when government restriction policies change. It can be adapted for any country, however this version is tailored specifically for Trinidad and Tobago. 

It uses a simple SIR model and changes in policy correspond to changes in Rnaught.

### Full Dashboard
![](https://raw.githubusercontent.com/rajeevratan84/Covid-19-Policy-Simulator/main/Full.png)

### Outputs
1. Cases over time, showing the periods of policy change (dotted line).

![](https://raw.githubusercontent.com/rajeevratan84/Covid-19-Policy-Simulator/main/Cases.png)

2. What happens if we didn't have any lockdowns?
3. 
![](https://raw.githubusercontent.com/rajeevratan84/Covid-19-Policy-Simulator/main/without_lockdown.png)

### Inputs (3 Phases in this version)
![](https://raw.githubusercontent.com/rajeevratan84/Covid-19-Policy-Simulator/main/policy.png)


**Parameters Used for SIR Model:**
Detected cases are usally ~10% of actual infections
R0 of 2.25
Gamma = 0.2, beta <- 4.5e-07
beta <- 4.5e-07
Population = 1,394,000
Initial Infect People = 2
Day 0 is 3rd March 2020
