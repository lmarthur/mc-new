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

echo "Input file: $input_file"
echo "Current directory: $(pwd)"
echo "Run directory: $run_dir"
echo "Independent variable: $param, iterated from $start to $end with $num steps"

# Search for run directory and if it doesn't exist, create it
if [ ! -d "$run_dir" ]; then
    mkdir "$run_dir"
fi

# Create input files in the run directory
for ((i=1; i<="$num"; i++))
do
    dir_name=input"$i"
    mkdir "$run_dir"/"$dir_name"
    cp "$input_file" "$run_dir"/"$dir_name"/input"$i".flap
done

# Run python script to modify the input files
python modify_input.py --dir="$run_dir" --param="$param" --start="$start" --end="$end" --num="$num"

# Run the MC-NEW executable on each input file
for dir in "$run_dir"/*/ ; do
    echo "$dir"

    ./MC-NEW < "$dir"/$(basename "$dir").flap > "$dir"/out

done

# Combine the outputs into a single data file

cd "$run_dir"
echo "Current directory: $(pwd)"
touch aero_coefs.dat
cd ..
cp cmy_alpha.plot "$run_dir"/
cp general.plot "$run_dir"/

# Uncomment the following when it is fully implemented
echo $(pwd)
python combine_data.py "$(pwd)" "$param"

# Run the plotting and analysis scripts

