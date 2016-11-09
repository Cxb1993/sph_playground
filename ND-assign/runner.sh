#!/bin/bash


dim=1
tasktype='infslb'
spacing='1 0.5 0.25 0.125 0.06 0.03 0.01'
ktype='n2w fab'
storebase='/Users/sergeibiriukov/_MoCA/Data'
for k in $ktype;
do
  errfname=$k'c-'$tasktype'-1-'$dim'D-err'
  header='ARG: lx3y. '$k'. '$tasktype'. Errors. k1=k2=1. t_max=5. '
  header=$header$dim'D. {| dx | n | err(t/3) | err(2t/3) | err(t) | t/3 | 2t/3 | t |}'
  `echo $header > $errfname`
done

it=0
for i in $spacing;
do
  for k in $ktype;
  do
    errfname=$k'c-'$tasktype'-1-'$dim'D-err'
    runcmd="time ./execute $dim $tasktype $i $errfname $k 1&>/dev/null"
    echo $runcmd
    `$runcmd`
    runcmd="mkdir -p $storebase/$k-$tasktype-1/$it-$i"
    `$runcmd`
    runcmd="mv output/* $storebase/$k-$tasktype-1/$it-$i"
    `$runcmd`
    echo "Done $k-$tasktype-1 for $i"
  done
  it=$((it+1))
done
