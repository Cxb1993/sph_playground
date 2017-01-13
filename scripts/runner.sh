#!/bin/bash

dimlist='1 2 3'
tasktype='hc-sinx'
ktype='n2w fab'
execnamelist='execute'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`

tfinish='.5'
spstart='.2'
spend='.01'
spstep='.01'
tstep=$spstart
flag='1'
spacing=""
while [[ $flag -eq "1" ]]; do
  spacing=$spacing" "$tstep
  tstep=`echo "$tstep - $spstep" | bc`
  flag=`echo "$tstep > $spend" | bc`
done
# spacing='0.2 0.19 0.2 0.19'
echo "Spacings: $spacing"

echo '' > runresult.info
`mkdir -p output`
for dim in $dimlist; do
  it=0
  for psp in $spacing; do
    for execname in $execnamelist; do
      for k in $ktype; do
        fullkernel=$k
        errfname=$dtprefix'-'$tasktype'-'$fullkernel'-'$dim'D'

        if [ "$it" = "0" ]; then
          header='ARG: lx3y3e3e. '$fullkernel'. '$tasktype'. t_max='$tfinish'. '
          header=$header$dim'D. {| dx | n | err(t/3) | bias(t/3) | t/3 | err(2t/3) | bias(2t/3) | 2t/3 | err(t) | bias(t) | t |}'
          `echo $header > $errfname`
        fi

        runcmd="time ./$execname $dim $tasktype $psp $errfname $k $tfinish &>/dev/null"
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
        echo -e "\nDone $tasktype $fullkernel $psp\n"
      done
    done
    it=$((it+1))
  done
done
`rm -rf output`
