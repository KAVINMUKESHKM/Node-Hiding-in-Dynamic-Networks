# BIH Dynamic Hiding Pipeline

This repository contains an implementation of the Based Importance Hiding (BIH) pipeline
that runs ArchAngel community detection on a dynamic network, applies rewiring to hide
overlapping nodes, and compares results.

## Quick Start

1. **Clone the repositories:**
   ```bash
   # Clone this repository
   git clone https://github.com/KAVINMUKESHKM/Node-Hiding-in-Dynamic-Networks.git
   cd Node-Hiding-in-Dynamic-Networks
   
   # Clone the ArchAngel library (required dependency)
   git clone https://github.com/GiulioRossetti/ANGEL.git angel
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Prepare your input data:**
   - Place your dynamic network file (3-column `.ncol` format: `nodeA\tnodeB\tsnapshot_id`) in the repo root.
   - Example dataset provided in `Datasets/Synthetic.ncol`

4. **Run the BIH pipeline:**
   ```bash
   python BihHiding.py Datasets/Synthetic.ncol
   ```

   Or with custom parameters:
   ```bash
   python BihHiding.py Datasets/Synthetic.ncol --max-candidate-sims 4 --hillclimb-max-iter 3
   ```

5. **Check results:**
   - Output directory: `BIH_Results/`
   - Modified network: `BIH_Results/networks/merged_snapshots_modified.ncol`
   - Comparison CSV: `BIH_Results/comparison/hiding_comparison.csv`

## Main Files

- `BihHiding.py` - Main BIH pipeline script
- `Datasets/` - Sample datasets for testing
- `requirements.txt` - Python dependencies
- `BIH_Results/` - Generated output directory (created on first run)

## Dependencies

This project requires the **ArchAngel** community detection library, which must be cloned separately:

```bash
git clone https://github.com/GiulioRossetti/ANGEL.git angel
```

The `angel` folder should be placed in the project root directory.

## CLI Usage

### Basic command:
```bash
python BihHiding.py <input_file> [options]
```

### Arguments:
- `input_file` (required) - Path to your .ncol file

### Options:
- `--base-path PATH` - Base directory for outputs (default: current directory)
- `--no-simulate-targets` - Disable simulation-based target selection
- `--debug-sim` - Write debug JSON files showing candidate scores
- `--max-candidate-sims N` - Limit candidate communities to simulate (default: all)
- `--no-hillclimb` - Disable hill-climb local search
- `--hillclimb-max-iter N` - Max hill-climb iterations (default: 3)

### Example commands:

**Run with default settings:**
```bash
python BihHiding.py Datasets/Synthetic.ncol
```

**Run with limited simulation and more hill-climbing:**
```bash
python BihHiding.py Datasets/Synthetic.ncol --max-candidate-sims 5 --hillclimb-max-iter 5
```

**Debug mode with detailed logs:**
```bash
python BihHiding.py Datasets/Synthetic.ncol --debug-sim
```

**Fast mode (no simulation, no hill-climb):**
```bash
python BihHiding.py Datasets/Synthetic.ncol --no-simulate-targets --no-hillclimb
```

## Output Format

After running the pipeline, `BIH_Results/` contains:

- **Per-snapshot community files:**
  - `original_detectionArchAngel_coms_<snapshot>.txt`
  - `modified_detectionArchAngel_coms_<snapshot>.txt`
  
- **Modified network:**
  - `networks/merged_snapshots_modified.ncol`

- **Comparison results:**
  - `comparison/hiding_comparison.csv` (per-node hiding success)
  - `comparison/sim_debug_snapshot_<id>.json` (if `--debug-sim` used)

- **Visualization:**
  - `hiding_cost_per_snapshot.png` (operations per snapshot chart)

## Notes

- All paths are relative to the repository root — no hardcoded paths
- Output directories are created automatically on first run
- **Important:** The ArchAngel library must be cloned separately into an `angel/` folder in the project root
- Simulation-based target selection and hill-climbing are enabled by default for better hiding success rates

