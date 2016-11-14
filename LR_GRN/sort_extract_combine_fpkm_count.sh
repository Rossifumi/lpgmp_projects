#!/bin/bash

# sort files by column 1
for file in *.txt;
    do 
        sort -k 1 "$file" > "${file/%txt/sort}";
    done

# extract certain column    
for file in *.sort;
    do
        awk '{print $6}' "$file" >"${file/%sort/sort.fpkm}";
        awk '{print $5}' "$file" >"${file/%sort/sort.count}";
    done

awk '{print $1}' 0h_hisat_sort.count.sort > sort.id;
rm *.sort;

# construct file names and combine columns
time_array=(0h 0.5h 1h 2h 3h 4h 5h 6h 8h 10h 12h 14h 16h 20h 24h 28h 32h 36h 40h 44h 48h 56h 64h 72h);

my_appendix_1='_hisat_sort.count.sort.fpkm';
my_appendix_2='_hisat_sort.count.sort.count';
my_files_1='';
my_files_2='';

for ((i=0; i< ${#time_array[@]}; i++));
    do
        my_files_1="$my_files_1 ${time_array[i]}$my_appendix_1";
        my_files_2="$my_files_2 ${time_array[i]}$my_appendix_2";
    done

paste sort.id $my_files_1 > Bna_exp_fpkm.tab;
paste sort.id $my_files_2 > Bna_exp_count.tab;

rm *sort.fpkm;
rm *sort.count;
rm sort.id;

