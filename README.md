# Sensivitiy Analysis

inputLATTE.py - reads in yangtab (reference density functional theory table), runs LATTE calculation on it, and then returns DFT energy/atom, TB energy/atom, DFT force, TB forces

parseParam.py - reads in original parameter files (.ref files) and prints them out to the LATTE input files (ppots.nonortho, bondints.nonortho, electrons.dat)

updateParam.py - reads in "log" which consists of new updated parameters (from Markov Chain Monte Carlo) and prints them out to LATTE input files

analysis.py - sensitivity analysis on output information
