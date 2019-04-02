# Adaptive integrate-and-fire model
Supporting datasets and code for the manuscript "Adaptive integrate-and-fire model reproduces the dynamics of olfactory receptor neuron responses in moth"

========================
 Data files description
========================

The minimum duration of a pheromone puff/blank was chosen either as 50 ms or 100 ms
(folders "50ms" and "100ms")

4 doses of pheromone were applied together with each minimum puff/blank duration
(folders "1pg", "10pg", "100pg", "1ng")

Complete data about a recording of a neuron 'IDENTIFIER' is given in 2 separate files:

IDENTIFIER_valve_states.txt
 - information on the stimulus time course
 - first column gives times when the stimulus was switched on or switched off
 - second column gives the type of the switch (1 = ON, -1 = OFF)

IDENTIFIER_spikes_times.txt
 - one column with times of spikes
