#! /bin/bash

for file in mlfsom*predin.txt; do
  tempnum=`echo $file | sed 's/[^0-9]*//g'`
  #tempcell=`cat $file | awk '/CELL/{print $2}'`
  mos=`cat $file | awk '/MOSAIC/{print $2}'`
  #echo $tempcell
  #dose=`echo $tempcell | awk -v cell_k=0.3 '{print ($1-77.25)/cell_k}'`
  #echo $dose
  rename "s/mlfsom_tempfile${tempnum}/mlfsom_tempfile_mos=${mos}/" *
done
