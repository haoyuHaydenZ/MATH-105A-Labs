%% Lab 4 Test script 

close all
clc
clear

% A new .mat file will be used to check your implementation of Gauss
% Elminination. Use this script to test your formatting. 

load data.mat 

tic
xhat = your_gauss_elim(AA,bb);
toc

norm(xx - xhat)
norm(AA*xhat - bb)
