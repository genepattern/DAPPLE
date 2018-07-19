#!/bin/bash

awk '$2!="NA"{print}' $1_seedScores | sort -g -k4 > tmp1
awk '$2=="NA"{print}' $1_seedScores > tmp2
cat tmp tmp2 > $1_seedScores
