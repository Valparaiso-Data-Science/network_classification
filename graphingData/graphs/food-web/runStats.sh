#!/bin/bash
# Run network stats on edge lists
for file in edges/*
do
  print $file
  ../../graphpack.exe -f $file --all -a 0 > "output/${file:6:15}.txt";
done