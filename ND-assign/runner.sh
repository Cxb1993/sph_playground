#!/bin/bash

dim=1
tasktype='hc-sinx'
# spacing=`get_spacing -1. 1. 0.001 10 20`
spacing='0.1 0.05 0.025 0.0125 0.006 0.003 0.001'
ktype='n2w fab'
kbase='c q'
# storebase='/Users/sergeibiriukov/_MoCA/Data'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`

brdx1='-1'
brdx2='1'
tfinish='.5'

for kb in $kbase; do
  execname='execute-'$kb
  `mkdir -p output`
  echo '' > runresult.info

  it=0
  for i in $spacing; do
    for k in $ktype; do
      fullkernel=$k$kb
      errfname=$dtprefix'-'$tasktype'-'$fullkernel'-'$dim'D'

      if [ "$it" = "0" ]; then
        header='ARG: lx3y3e3e. '$fullkernel'. '$tasktype'. ['$brdx1';'$brdx2']. t_max='$tfinish'. '
        header=$header$dim'D. {| dx | n | err(t/3) | bias(t/3) | t/3 | err(2t/3) | bias(2t/3) | 2t/3 | err(t) | bias(t) | t |}'
        `echo $header > $errfname`
      fi

      runcmd="time ./$execname $dim $tasktype $i $errfname $k $brdx1 $brdx2 $tfinish &>/dev/null"
      echo $runcmd
      runresult=`$runcmd`
      echo "$runresult" >> runresult.info
      runcmd="mkdir -p $storebase/$fullkernel/$it-$i"
      `$runcmd`
      runcmd="mv output/* $storebase/$fullkernel/$it-$i"
      `$runcmd`
      echo -e "\nDone $tasktype $fullkernel $i\n"
    done
    it=$((it+1))
  done
done
