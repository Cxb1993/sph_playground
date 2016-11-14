#!/bin/bash

dim=1
tasktype='infslb'
spacing='1 0.5 0.25 0.125 0.06 0.03 0.01'
ktype='n2w fab'
storebase='/Users/sergeibiriukov/_MoCA/Data'
dtprefix=`date +%Y%m%d%H%M`
taskdir="$2-infslb-1"
brdx1='-10'
brdx2='10'
tfinish='5'
execname=$1

`mkdir -p output`
for k in $ktype;
do
  errfname=$dtprefix'-'$k$taskdir'-'$dim'D-err'
  header='ARG: lx3y3e3e. '$k'. '$tasktype'. Errors. k1=k2=1. t_max='$tfinish
  header=$header$dim'D. {| dx | n | err(t/3) | bias(t/3) | t/3 | err(2t/3) | bias(2t/3) | 2t/3 | err(t) | bias(t) | t |}'
  `echo $header > $errfname`
done
echo '' > runresult.info

it=0
for i in $spacing; do
  for k in $ktype; do
    errfname=$dtprefix'-'$k$taskdir'-'$dim'D-err'
    runcmd="time ./$execname $dim $tasktype $i $errfname $k $brdx1 $brdx2 $tfinish &>/dev/null"
    echo $runcmd
    runresult=`$runcmd`
    echo "$runresult" >> runresult.info
    runcmd="mkdir -p $storebase/$k$taskdir/$it-$i"
    `$runcmd`
    runcmd="mv output/* $storebase/$k$taskdir/$it-$i"
    `$runcmd`
    echo -e "\nDone $k$taskdir for $i\n"
  done
  it=$((it+1))
done
