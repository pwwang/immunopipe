#!/usr/bin/env bash

# run the command and capture the stdout
out=$(command $@)

echo "$out"

if [[ $(echo "$out" | grep -c reduce) -gt 0 ]]; then
    echo ""
    echo "Original Rscript failed, try reduce = 0 ..."
    echo "Modifying vj_pairing_plot.r ..."
    sed -i 's/mat, annotationTrack/mat, reduce = 0, annotationTrack/' vj_pairing_plot.r
    cmd=$(echo "$out" | grep "\\[RUtil\\] Executing" | sed 's/\[RUtil\] Executing //')
    echo "Run again: $cmd"
    command $cmd
else
    echo ""
    echo "Not a 'reduce' failure for original script, nothing to do."
fi
