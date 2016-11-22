#!/bin/bash

dim=1
tasktype='hc-sinx'
spacing='0.1 0.05 0.025 0.0125 0.006 0.003'
ktype='n2w fab'
kbase='c q'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`

brdx1='-1'
brdx2='1'
tfinish='.5'

spstart=0.1
spend=0.
spstep=0.01
tstep=$spstart
flag=1
spacing=""
while [[ $flag -eq "1" ]]; do
  spacing=$spacing" "$tstep
  tstep=`echo "$tstep - $spstep" | bc`
  flag=`echo "$tstep > $spend" | bc`
done
echo "Spacings: $spacing"

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

      errfname=$dtprefix'-'$tasktype'-'$fullkernel'-'$dim'D'
      runcmd="time ./$execname $dim $tasktype $i $errfname $k $brdx1 $brdx2 $tfinish &>/dev/null"
      echo $runcmd
      runresult=`echo '\n' | $runcmd`
      echo "$runresult" >> runresult.info

      itsize=${#it}
      itspac=`tail -1 $errfname | awk '{print$1}'`
      itspac=`printf %f $itspac`
      if [ $itsize -eq '1' ]; then
        iti="0"$it
      else
        iti=$it
      fi
      runcmd="mkdir -p $storebase/$fullkernel/$iti-$itspac"
      `$runcmd`
      runcmd="mv output/* $storebase/$fullkernel/$iti-$itspac"
      `$runcmd`
      echo -e "\nDone $tasktype $fullkernel $i\n"
    done
    it=$((it+1))
  done
done

`rm -rf output`
