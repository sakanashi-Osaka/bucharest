#!/bin/bash

# Check if the script received one argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <4-digit-run-number>"
  exit 1
fi

# Assign the first argument to the run variable
run=$1

# Validate that the run number is a 4-digit number
if [[ ! $run =~ ^[0-9]{4}$ ]]; then
  echo "Error: The run number must be a 4-digit number."
  exit 1
fi

# Loop for i from 0 to 15
for i in {0..15}; do
  # Loop for j from 0 to 15
  for j in {0..15}; do
    # Execute the command
    echo "Running: ../calib $run $i $j"
    ./calib $run $i $j
  done
done
