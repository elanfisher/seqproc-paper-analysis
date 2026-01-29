#!/bin/bash
# Debug 10x extraction
SPLITCODE="/home/ubuntu/combine-lab/splitcode/build/src/splitcode"
CONFIG="configs/splitcode/10x_v2.config"
INPUT="head_10x.fq"

echo "--- Testing with current config: ---"
cat $CONFIG
echo "------------------------------------"

$SPLITCODE --x-only -c $CONFIG $INPUT > debug_10x_out.txt 2> debug_10x.log

echo "Exit code: $?"
echo "Log content:"
cat debug_10x.log
echo "Output content (head):"
head debug_10x_out.txt
wc -l debug_10x_out.txt
