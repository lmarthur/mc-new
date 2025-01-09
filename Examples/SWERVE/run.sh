#!/bin/bash

while getopts i:d:p:s:e:n: flag
do
    case "${flag}" in
	i)input_file=${OPTARG};;
	d)run_dir=${OPTARG};;
	p)param=${OPTARG};;
	s)start=${OPTARG};;
	e)end=${OPTARG};;
	n)num=${OPTARG};;
    esac
done

echo "$input_file"
echo "$run_dir"
echo "$param"
echo "$start"
echo "$end"
echo "$num"

cd "$run_dir"
rm -r *
cd ..

for ((i=1; i<="$num"; i++))
do
    dir_name=input"$i"
    mkdir "$run_dir"/"$dir_name"
    cp "$input_file" "$run_dir"/"$dir_name"/input"$i".flap
done

python3 modify_input.py --dir="$run_dir" --param="$param" --start="$start" --end="$end" --num="$num"

for dir in "$run_dir"/*/ ; do
    echo "$dir"
    cp MC-NEW "$dir"/
    cd "$dir"
    ./MC-NEW < $(basename "$dir").flap > out
    cd ..
    cd ..
done

cd "$run_dir"
touch all_aero_coefs.dat
cd ..
cp cmy_alpha.plot "$run_dir"/
cp general.plot "$run_dir"/

python3 combine_data.py "$run_dir" "$param"
