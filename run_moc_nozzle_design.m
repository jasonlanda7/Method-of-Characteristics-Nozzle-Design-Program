clc; clear; close all;
disp('=== METHOD OF CHARACTERISTICS NOZZLE DESIGN ===');
P0 = input('Enter chamber pressure P0 [Pa]: ');
g  = input('Enter ratio of specific heats (gamma): ');
TR = input('Enter throat radius TR [mm]: ');
out = moc_bell_nozzle(P0, g, TR);
