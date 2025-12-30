# UAV Base Station Positioning in Urban Scenarios

This project efficiently and accurately solves the problem of UAV base station positioning and site selection in complex and dense urban scenarios by constructing a Realistic Coverage Model and the GeoAIS optimization algorithm.

## Implementation Steps

### Step 1: Generate UAV Base Station Candidate Set
Run the `1_ABSCandidateSet.py` file:
- **Input**: `ABS_selection.3idx` index file (refer to the Q-view algorithm for index, or you can download the preprocessed data from this link:https://pan.baidu.com/s/1jm8V6GVaoDZZUKsiDlvhoA, Extraction code：2025).
- **Output**: Candidate information for 3498 UAV base stations (including base station ID, XYZ coordinates, LOS radius, and NLOS radius), stored in `/data/coverage_results.csv`.

### Step 2: Build Realistic Coverage Model
Run the `2_realistic_coverage_model.py` file:
- **Input**: `/data/coverage_results.csv` (generated in Step 1).
- **Output**: Realistic coverage data of the base stations (stored in binary format), saved to `/data/coverage.bin`.

### Step 3: Execute GeoAIS Optimization Algorithm
Run the `GeoAIS.py` file:
- **Input**: `/data/coverage.bin` and `/data/coverage_results.csv`.
- **Output**: Actual coverage, overlap, corresponding base station IDs and running time for the required number of base stations (e.g., 20 ABSs).
