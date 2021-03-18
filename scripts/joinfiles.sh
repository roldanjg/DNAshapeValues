#!/usr/bin/env bash

awk 'NF' tomato_SL4.0ch00.Roll.wig tomato_SL4.0ch01.Roll.wig tomato_SL4.0ch02.fa.Roll.wig tomato_SL4.0ch03.Roll.wig tomato_SL4.0ch04.Roll.wig tomato_SL4.0ch05.Roll.wig tomato_SL4.0ch06.Roll.wig tomato_SL4.0ch07.Roll.wig tomato_SL4.0ch08.Roll.wig tomato_SL4.0ch09.Roll.wig tomato_SL4.0ch10.Roll.wig tomato_SL4.0ch11.Roll.wig tomato_SL4.0ch12.Roll.wig > SL4.0.Roll.wig
../../programs/wigToBigWig SL4.0.Roll.wig genome.index.txt SL4.0.Roll.wig.bw
awk 'NF' tomato_SL4.0ch00.MGW.wig tomato_SL4.0ch01.MGW.wig tomato_SL4.0ch02.fa.MGW.wig tomato_SL4.0ch03.MGW.wig tomato_SL4.0ch04.MGW.wig tomato_SL4.0ch05.MGW.wig tomato_SL4.0ch06.MGW.wig tomato_SL4.0ch07.MGW.wig tomato_SL4.0ch08.MGW.wig tomato_SL4.0ch09.MGW.wig tomato_SL4.0ch10.MGW.wig tomato_SL4.0ch11.MGW.wig tomato_SL4.0ch12.MGW.wig > SL4.0.MGW.wig
../../programs/wigToBigWig SL4.0.MGW.wig genome.index.txt SL4.0.MGW.wig.bw
awk 'NF' tomato_SL4.0ch00.ProT.wig tomato_SL4.0ch01.ProT.wig tomato_SL4.0ch02.fa.ProT.wig tomato_SL4.0ch03.ProT.wig tomato_SL4.0ch04.ProT.wig tomato_SL4.0ch05.ProT.wig tomato_SL4.0ch06.ProT.wig tomato_SL4.0ch07.ProT.wig tomato_SL4.0ch08.ProT.wig tomato_SL4.0ch09.ProT.wig tomato_SL4.0ch10.ProT.wig tomato_SL4.0ch11.ProT.wig tomato_SL4.0ch12.ProT.wig > SL4.0.ProT.wig
../../programs/wigToBigWig SL4.0.ProT.wig genome.index.txt SL4.0.ProT.wig.bw
awk 'NF' tomato_SL4.0ch00.HelT.wig tomato_SL4.0ch01.HelT.wig tomato_SL4.0ch02.fa.HelT.wig tomato_SL4.0ch03.HelT.wig tomato_SL4.0ch04.HelT.wig tomato_SL4.0ch05.HelT.wig tomato_SL4.0ch06.HelT.wig tomato_SL4.0ch07.HelT.wig tomato_SL4.0ch08.HelT.wig tomato_SL4.0ch09.HelT.wig tomato_SL4.0ch10.HelT.wig tomato_SL4.0ch11.HelT.wig tomato_SL4.0ch12.HelT.wig > SL4.0.HelT.wig
../../programs/wigToBigWig SL4.0.HelT.wig genome.index.txt SL4.0.HelT.wig.bw