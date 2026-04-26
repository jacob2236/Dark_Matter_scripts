#!/bin/bash

BASE_DIR="scan_results"

for dir in "$BASE_DIR"/*/; do
   TARGET="${dir%/}/Events/run_01"
   
   if [ -d "$TARGET" ]; then
      FILE=$(find "$TARGET" -maxdepth 1 -name "*banner.txt" | head -n 1)
      
      if [ -f "$FILE" ]; then
         number_of_events=$(awk '$2 =="Number" && $3=="of" && $4=="Events" {print $NF; exit}' "$FILE")
         cross_section_value=$(awk '$2=="Integrated" && $3=="weight" {print $NF; exit}' "$FILE")
         echo "File: $FILE"
         echo "Number of events: $number_of_events"
         echo "Cross Section (σ): $cross_section_value"
      else
         echo "Banner file not found in $TARGET"
      fi
   else
      echo "Directory $TARGET not found"
   fi
   echo
done