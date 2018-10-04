#!/bin/bash
# Run network stats on edge lists
for file in edges/*
do
  echo "${file:11:4}"
  ../../../graphpack -f $file -a 0 > "output/${file:11:4}.txt";
done