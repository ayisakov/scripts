#!/bin/bash

# This script will set the first non-LVDS output (monitor) as the
# primary monitor if multiple outputs exist

MONITORS=$(xrandr -q | grep " connected" | awk '{ print $1 }')
if [ -z "$MONITORS" ]; then
    echo "No monitors detected" >2
    exit 2
fi

# set the first output that is not the LVDS (built-in display) output
# as primary if more than one output is connected
if [ $(echo "$MONITORS" | wc -l) -eq 1 ]; then
    exit 0 # nothing to do if only one output is connected
fi

while IFS= read -r line; do
    if [ "$line" != "LVDS" ]; then
        xrandr --output $line --primary
        break
    fi
done < <(printf '%s\n' "$MONITORS")

