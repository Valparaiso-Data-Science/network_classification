#!/bin/bash
while read LINE
	do wget -nc "http://nrvis.com/download/data/$LINE.zip"
done < list
unzip -o \*.zip
rm -rf ./*.zip

for i in *.mtx; do
	../graphpack -f $i -a 0 -t 4
done
