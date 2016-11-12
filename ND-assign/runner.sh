#!/bin/bash

dim=1
tasktype='infslb'
spacing='1 0.5 0.25 0.125 0.06 0.03 0.01'
ktype='n2w fab'
storebase='/Users/sergeibiriukov/_MoCA/Data'
dtprefix=`date +%Y%m%d%H%M`
taskfull='c-'$tasktype'-10'

`mkdir -p output`
for k in $ktype;
do
  errfname=$dtprefix'-'$k$taskfull'-'$dim'D-err'
  header='ARG: lx3y3e3e. '$k'. '$tasktype'. Errors. k1=k2=1. t_max=5. '
  header=$header$dim'D. {| dx | n | err(t/3) | bias(t/3) | t/3 | err(2t/3) | bias(2t/3) | 2t/3 | err(t) | bias(t) | t |}'
  `echo $header > $errfname`
done
echo '' > runresult.info

it=0
for i in $spacing; do
  for k in $ktype; do
    errfname=$dtprefix'-'$k$taskfull'-'$dim'D-err'
    runcmd="time ./execute $dim $tasktype $i $errfname $k 1&>/dev/null"
    echo $runcmd
    runresult=`$runcmd`
    echo "$runresult" >> runresult.info
    runcmd="mkdir -p $storebase/$k$taskfull/$it-$i"
    `$runcmd`
    runcmd="mv output/* $storebase/$k$taskfull/$it-$i"
    `$runcmd`
    echo -e "\nDone $k$taskfull for $i\n"
  done
  it=$((it+1))
done
