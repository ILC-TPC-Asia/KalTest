#!/bin/sh
nevt=5000
pt=-100.0

t0=000
./EXKalTest -b $nevt $pt $t0
mv h.root p${pt}t${t0}.root

t0=007
./EXKalTest -b $nevt $pt $t0
mv h.root p${pt}t${t0}.root

t0=014
./EXKalTest -b $nevt $pt $t0
mv h.root p${pt}t${t0}.root

t0=021
./EXKalTest -b $nevt $pt $t0
mv h.root p${pt}t${t0}.root
