#!/bin/bash

for i in {1..100}
do 
	bsub -q hour -W 4:0 -o bsuboutput -J 700runs python ../BuildNetwork_STRING.py heiko_input keyword=heiko_input_700_P$i database=700 permute=5
done

bsub -q hour -o bsuboutput -w 'ended(700runs)' python mergeDAPPLE.py


for i in {1..100}
do 
	bsub -q hour -W 4:0 -o bsuboutput -J 400runs python ../BuildNetwork_STRING.py heiko_input keyword=heiko_input_400_P$i database=400 permute=5
done

#for i in {1..500}
#do 
#	bsub -q week -o bsuboutput -J fullruns python ../BuildNetwork_STRING.py heiko_input keyword=heiko_input_full_P$i database=full permute=1
#done

bsub -q hour -o bsuboutput -w 'ended(400runs)' -J merge python mergeDAPPLE.py

bsub -q hour -o bsuboutput -w 'ended(merge)' bash cleanup.bash
