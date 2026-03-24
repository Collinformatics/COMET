#!/bin/bash

# Parameters
colorR='\e[91m'
color='\e[34m'
rst='\e[0m'

# Define port
port="$1"
if [ -z $1 ]; then
  port="9090"
fi
printf "Checking processes at port: %b%b%b\\n" "$colorR" "$port" "$rst"

# List processes
killProc=false
procs=$(lsof -i :$port)
if [[ -n "$procs" ]]; then
  echo -e "$procs"
  echo "$procs" | awk 'NR > 1 {print $1, $2}' > proc.txt
  p=()
  while IFS= read -r line; do
      read -ra parts <<< "$line"  # Split line into words
      p+=("${parts[@]}")
  done < proc.txt
  rm proc.txt

  # Kill python processes
  for (( i=0; i<${#p[@]}; i+=2 )); do
    x=${p[i]}
    if [[ $x == *python* ]]; then
      killProc=true
      pid=${p[i+1]}
      printf "Kill PID: %b%b%b\\n" "$colorR" "$pid" "$rst"
      kill $pid
    fi
  done
fi
if [[ $killProc != true ]]; then
  echo "No processes were terminated"
fi
