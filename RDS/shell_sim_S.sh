#!/bin/bash
# Define the directory to run the simulations:
RUN_DIR="/home/seanleishman/Personal/claire/Simulations"
# Define the main directory containing all the case folders:
TARGET_DIR="/home/seanleishman/Personal/claire/RDS/Simulations"
# List of known case numbers
CASE_NUMBERS="283"
#2025-04-13 269 270 271
#2025-03-20 247 248 249 250 252 253 256 257 258 261 262 264 265 266 254 259 267 255 260 263 268 248 251
##2025-03-20 215 216 217 218 219 220 221 222 224 223 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 (aborted during 241)
#2025-03-20 212 213 214
#2025-02-228 194
#2025-02-18 207
#2025-02-03 204 205 206
#2024-01-26 1641 1932 1654
#2025-01-23 1501 1641 1932 1654 1931 203
#202 203 1881 1882 1931 1932 204 205 206 207
#184 IS ACTUALLY 172 REMEMBER TO CHANGE; I HAVE REDEFINED PARAMS IN 184 CAUSE THIS FILE WILL BE WAY TOO MASSIVE
# 2_test 3_test

LOG_FILE="/home/seanleishman/Personal/claire/simulation_log.txt"
exec > >(tee -a "$LOG_FILE") 2>&1

for CASE_NUM in $CASE_NUMBERS; do
    # Check and cleanup RUN_DIR to ensure only ADEv8.cc is in the RUN_DIR ***CHANGE
    echo "Cleaning up $RUN_DIR..."
    find "$RUN_DIR" -type f ! -name 'ADEv8.cc' -exec rm {} +
    find "$RUN_DIR" -type d ! -path "$RUN_DIR" -exec rm -r {} +
    #
    # Specify the case directory, from which parameter file is extracted and output files are uploaded after simulation.
    CASE_DIR="$TARGET_DIR/Case_$CASE_NUM"
    #cd "$CASE_DIR"
    #
    #Print line to indicate the case that is processed
    echo "Processing Case $CASE_NUM"
    #
    # Check for the existence of 'data.dat' in the case directory
    if [ -f "$CASE_DIR/Data.dat" ]; then
        echo "Found 'Data.dat' in $CASE_DIR. Skipping Case $CASE_NUM."
        continue  # Skip to the next iteration of the loop
    fi
    #
    #Check and transfer the parameters file
    PARAM_FILE="$CASE_DIR/parameters.dat"
    if [ -f "$PARAM_FILE" ]; then
        mv "$PARAM_FILE" "$RUN_DIR/parameters.dat"
        echo "parameters.dat moved to $RUN_DIR for case $CASE_NUM."
    else
        echo "parameters.dat not found in $CASE_DIR. Skipping Case $CASE_NUM"
        continue  # Skip to the next case
    fi
    #
    #if [ ! -f "$PARAM_FILE" ]; then
     #   echo "parameters.dat not found in $CASE_DIR"
     #   continue  # Skip this case
    #fi
  #  mv "$PARAM_FILE" "$RUN_DIR/parameters.dat"
#Navigate to the run directory to run the case
    cd "$RUN_DIR"
    echo "Currently in $(pwd)"
#Compile the simulation program REMEMBER TO CHANGE V888
    g++ ADEv8.cc -O3 -o binary
    if [ $? -ne 0 ]; then
        echo "Compilation failed for case $CASE_NUM."
        echo "parameters.dat moved to $CASE_DIR"
        mv "$RUN_DIR/parameters.dat" "$CASE_DIR"
        continue  # Skip to the next case
    fi
#Run the simulation
    ./binary
    if [ $? -ne 0 ]; then
        echo "Simulation failed in $CASE_DIR"
        mv "$RUN_DIR/parameters.dat" "$CASE_DIR"
        continue  # Skip to the next case
    fi
#Organise output files into case specific directory
    mv vtk_fluid "$CASE_DIR"
    mv Data.dat "$CASE_DIR"
    mv parameters.dat "$CASE_DIR"
    mv ReactionProfile.dat "$CASE_DIR"
    mv binary "$CASE_DIR"
    echo "Case $CASE_NUM processed successfully."

done
