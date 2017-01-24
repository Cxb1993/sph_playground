#!/bin/bash

dimlist='1 2 3'
tasktype='diff-laplace'
ktype='n2w fab'
execnamelist='execute'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`
kernelPrefix='quintic'

tfinish='100'
spstart='.2'
spend='.00'
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
        fullkernel=$kernelPrefix' '$k
        errfname=$dtprefix'-'$tasktype'-'$dim'D-'$k

        if [ "$it" = "0" ]; then
          header='ARG: xey. '$fullkernel'. '$tasktype'. '
          header=$header$dim'D. {| dx | partN | err l2 | 2nd err term | hfac |}'
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
        moveto="$storebase/$dim""D-""$k/$iti-$itspac"
        runcmd="mkdir -p $moveto"
        `$runcmd`
        runcmd="mv output/* $moveto"
        `$runcmd`
        echo -e "\nDone $tasktype $fullkernel $psp\n"
      done
    done
    it=$((it+1))
  done
done
`rm -rf output`
