#!/bin/bash

dimlist='1 2 3'
# dimlist='1 2'
# dimlist='2'
# dimlist='3'
tasktype='chi-laplace'
ktype='n2w fab'
# ktype='n2w'
# ktype='fab'
# kbase='q c'
kbase='cubic'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`
kernelprefix=$kbase

# This is spacing now
tfinish='-1'
realspacing='.002'
spstart='1'
spend='2'
spstep='.01'
tstep=$spstart
flag='1'
# in this case it is "sk" multiplier of number of neighbours
spacing=""
while [[ $flag -eq "1" ]]; do
  spacing=$spacing" "$tstep
  tstep=`echo "$tstep + $spstep" | bc`
  flag=`echo "$tstep < $spend" | bc`
done
echo "Spacings: $spacing"

echo '' > result.info
`mkdir -p output`
for dim in $dimlist; do
  it=0
  for psp in $spacing; do
    for kb in $kbase; do
      execname='execute'
      for k in $ktype; do
        fullkernel=$kernelprefix' '$k
        errfname=$dtprefix'-'$tasktype'-'$dim'D-'$k

        if [ "$it" = "0" ]; then
          header='ARG: lx3y3e3e. '$fullkernel'. '$tasktype'. t_max='$tfinish'. '
          header=$header$dim'D. {| dx | n | 0 | 0 | h | chi11 | chi12 | chi13 | chi21 | chi22 | chi23 | chi31 | chi32 | chi33 |}'
          `echo $header > $errfname`
        fi

        runcmd="time ./$execname $dim $tasktype $realspacing $errfname $k $tfinish $psp &>/dev/null"
        echo $runcmd
        runresult=`echo '\n' | $runcmd`
        echo "$runresult" >> result.info

        itsize=${#it}
        itspac=`tail -1 $errfname | awk '{print$NF}'`
        itspac=`printf %f $itspac`
        if [ $itsize -eq '1' ]; then
          iti="00"$it
        elif [ $itsize -eq '2' ]; then
          iti="0"$it
        else
          iti=$it
        fi
        # `mkdir -p $storebase/$dim""D-""$k/`
        # moveto="$storebase/$dim""D-""$k/$iti-$itspac.zip"
        # runcmd="zip -9 $moveto ./output/*"
        # runresult=`$runcmd`
        # echo "$runresult" >> result.info
        # `rm -rf output/*`
        # echo -e "\nDone $tasktype $fullkernel $psp\n"
      done
    done
    it=$((it+1))
  done
done

`rm -rf output`
