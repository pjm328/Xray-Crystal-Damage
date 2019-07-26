#! /bin/bash

for file in mlfsom*predin.txt; do
  tempnum=`echo $file | sed 's/[^0-9]*//g'`
  tempmos=`cat $file | awk '/MOSAIC/{print $2}'`
  rename "s/mlfsom_tempfile${tempnum}/mlfsom_tempfile_mos=${tempmos}/" *
done
