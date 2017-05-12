#!/bin/bash

env | grep OMP_NUM_THREADS
dimlist='1 2 3'
# dimlist=$1
# dimlist='1 2'
# tasktype='diff-laplace'
tasklist='chi-graddiv chi-laplace'
# tasktype='chi-laplace'
# tasktype=$2
# ktype='n2w'
# ktype='fab'
# ktype='2nw'
kernellist='n2w fab'
execname='execute'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`
# kernelPrefix='quintic'
kernelPrefix='cubic'
# kernelPrefix='mgauss'
# kernelPrefix='sinc'
difftype='diff'
# difftype='symm'
#----------------------------
calcinfluence='no'
tfinish='100'
spstart='1.'
spend='2.'
# spend='3.'
spstep='.01'
tstep=$spstart
flag='1'
spacing=""
while [[ $flag -eq "1" ]]; do
  spacing=$spacing" "$tstep
  tstep=`echo "$tstep + $spstep" | bc`
  flag=`echo "$tstep < $spend" | bc`
done
# spacing='0.2'
echo "Spacings: $spacing"

echo '' > result.info
`mkdir -p output/1 output/2 output/3`
for dim in $dimlist; do
  for curtask in $tasklist; do
    for curkern in $kernellist; do
      fullkernel=$kernelPrefix' '$curkern
      errfname=$dtprefix'-'$curtask'-'$difftype'-'$dim'D-'$curkern
      it=0

      for psp in $spacing; do
        if [ "$it" = "0" ]; then
          header='ARG. '$fullkernel'. '$curtask'. '$difftype'. '
          header=$header$dim'D. {| dx | partN | err l2 | 2nd err term | hfac |}'
          `echo $header > $errfname`
        fi

        itsize=${#it}
        # itspac=`tail -1 $errfname | awk '{print$1}'`
        itspac=`printf %f $psp`
        if [ $itsize -eq '1' ]; then
          iti="00"$it
        else
          if [ $itsize -eq '2' ]; then
            iti="0"$it
          else
            iti=$it
          fi
        fi

        runcmd="time ./$execname --dim $dim --tasktype $curtask --spacing 0.06 \
                      --errfilename $errfname --kerneltype $curkern --tfinish $tfinish \
                      --hfac $psp --difftype $difftype --silent yes"

        if [ "$calcinfluence" = "yes" ]; then
          `echo "hfac=$psp\n" > output/$dim/influence$iti.info`
          runcmd="$runcmd --kerninfluencefile output/$dim/influence$iti.info"
        fi

        runcmd="$runcmd &>/dev/null"
        echo $runcmd
        runresult=`echo '\n' | $runcmd`
        echo "$runresult" >> result.info
        # `mkdir -p $storebase/$dim""D-""$k/`
        # moveto="$storebase/$dim""D-""$k/$iti-$itspac.zip"
        # runcmd="zip -9 $moveto ./output/*"
        # runresult=`$runcmd`
        # echo "$runresult" >> result.info
        # `rm -rf output/*`
        # echo -e "\nDone $tasktype $fullkernel $psp\n"
        it=$((it+1))
      done
    done
  done
done
runcmd="zip -9 -r output.zip ./output/*"
cmdres=`$runcmd`
# `rm -rf output`
echo "$errfname"
