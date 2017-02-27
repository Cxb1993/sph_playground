#!/bin/bash

env | grep OMP_NUM_THREADS
# dimlist='1 2 3'
dimlist=$1
# tasktype='diff-graddiv'
tasktype=$2
# ktype='n2w fab'
ktype='n2w'
execnamelist='execute'
storebase=`pwd`
dtprefix=`date +%Y%m%d%H%M`
# kernelPrefix='quintic'
kernelPrefix='cubic'
# kernelPrefix='mgauss'
# kernelPrefix='sinc'

tfinish='100'
spstart='1.'
spend='2.'
# spend='3.'
spstep='.02'
tstep=$spstart
flag='1'
spacing=""
while [[ $flag -eq "1" ]]; do
  spacing=$spacing" "$tstep
  tstep=`echo "$tstep + $spstep" | bc`
  flag=`echo "$tstep < $spend" | bc`
done
# spacing='0.2 0.19 0.2 0.19'
echo "Spacings: $spacing"

echo '' > result.info
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

        runcmd="time ./$execname $dim $tasktype 0.06 $errfname $k $tfinish $psp &>/dev/null"
        echo $runcmd
        runresult=`echo '\n' | $runcmd`
        echo "$runresult" >> result.info

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
runcmd="zip -9 output.zip ./output/*"
cmdres = `$runcmd`
`rm -rf output`
echo "$errfname"
