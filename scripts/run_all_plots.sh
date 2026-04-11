#!/usr/bin/env bash
# Batch plotting script - runs all plotting scripts for each dataset

# Configuration
OUT_DIR_PLOTS_SMALL="${OUT_DIR_PLOTS_SMALL:-plots/small}"
OUT_DIR_PLOTS_SYNTHETIC="${OUT_DIR_PLOTS_SYNTHETIC:-plots/synthetic}"
IND_DIR_RESULTS_SMALL="${IND_DIR_RESULTS_SMALL:-results/small}"
IND_DIR_RESULTS_SYNTHETIC="${IND_DIR_RESULTS_SYNTHETIC:-results/synth}"

# Plotting scripts
PLOT_SCRIPTS=(
    "scripts/plot_all.py"
    "scripts/plot_cycle_removal.py"
    "scripts/plot_hst.py"
    "scripts/plot_mendel.py"
)

# Dry-run mode (set to true to just see what would be processed)
DRY_RUN="${DRY_RUN:-false}"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Hardcoded datasets
SMALL_DATASETS=("BackBone" "Rocketfuel" "SNDLib")
SYNTHETIC_DATASETS=("fatclique" "2dgrids" "expander" "fattree")

# Helper function to run plots for a single CSV with all plot scripts
run_plots_for_csv() {
    local csv_file="$1"
    local output_dir="$2"
    local dataset_name="$3"

    if [[ ! -f "$csv_file" ]]; then
        echo -e "${RED}✗ CSV file not found: $csv_file${NC}"
        return 1
    fi

    echo -e "${BLUE}Processing: $csv_file${NC}"
    echo -e "${BLUE}  → Output: $output_dir${NC}"

    # In dry-run mode, just skip actual execution
    if [[ "$DRY_RUN" == "true" ]]; then
        echo -e "${YELLOW}(DRY-RUN: would process this dataset with all plot scripts)${NC}"
        return 0
    fi

    # Create output directory
    mkdir -p "$output_dir"

    local script_count=0
    local script_failed=0

    # Run each plotting script
    for plot_script in "${PLOT_SCRIPTS[@]}"; do
        if [[ ! -f "$plot_script" ]]; then
            echo -e "${YELLOW}  ⚠ Script not found: $plot_script${NC}"
            continue
        fi

        echo -e "${BLUE}  Running: $plot_script${NC}"

        if python3 "$plot_script" --input "$csv_file" --output "$output_dir"; then
            echo -e "${GREEN}  ✓ $plot_script completed${NC}"
            ((script_count++)) || true
        else
            echo -e "${RED}  ✗ $plot_script failed${NC}"
            ((script_failed++)) || true
        fi
    done

    if [[ $script_failed -eq 0 ]]; then
        echo -e "${GREEN}✓ Successfully generated all plots for $dataset_name${NC}"
        return 0
    else
        echo -e "${YELLOW}⚠ Some plots failed for $dataset_name ($script_failed/$((script_count + script_failed)))${NC}"
        return 1
    fi
}

# Main function to process datasets
process_datasets() {
    local dataset_type="$1"
    local input_dir="$2"
    local output_base_dir="$3"
    shift 3
    local datasets=("$@")

    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Processing $dataset_type datasets${NC}"
    echo -e "${BLUE}========================================${NC}"

    local processed_count=0
    local failed_count=0

    # Process each dataset
    for dataset_name in "${datasets[@]}"; do
        echo -e "${BLUE}Checking dataset: $dataset_name${NC}"

        # Determine CSV file path
        local dataset_dir="$input_dir/$dataset_name"
        local csv_file=""

        if [[ -f "$dataset_dir/combined.csv" ]]; then
            csv_file="$dataset_dir/combined.csv"
        elif [[ -f "$dataset_dir/${dataset_name}.csv" ]]; then
            csv_file="$dataset_dir/${dataset_name}.csv"
        else
            echo -e "${YELLOW}⚠ No CSV file found for $dataset_name${NC}"
            ((failed_count++)) || true
            continue
        fi

        local output_dir="$output_base_dir/$dataset_name/"

        # Run all plotting scripts
        if run_plots_for_csv "$csv_file" "$output_dir" "$dataset_name"; then
            ((processed_count++)) || true
        else
            ((failed_count++)) || true
        fi
    done

    echo ""
    echo -e "${BLUE}Summary for $dataset_type:${NC}"
    echo -e "  ${GREEN}Processed: $processed_count${NC}"
    if [[ $failed_count -gt 0 ]]; then
        echo -e "  ${RED}Failed: $failed_count${NC}"
    fi
}

# Main execution
echo ""
echo -e "${YELLOW}Starting plot generation for all datasets...${NC}"
echo ""

# Process small datasets
process_datasets "Small" "$IND_DIR_RESULTS_SMALL" "$OUT_DIR_PLOTS_SMALL" "${SMALL_DATASETS[@]}"

# Process synthetic datasets
process_datasets "Synthetic" "$IND_DIR_RESULTS_SYNTHETIC" "$OUT_DIR_PLOTS_SYNTHETIC" "${SYNTHETIC_DATASETS[@]}"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}✓ All plots generated successfully!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Output directories:"
echo "  Small datasets:     $OUT_DIR_PLOTS_SMALL"
echo "  Synthetic datasets: $OUT_DIR_PLOTS_SYNTHETIC"
echo ""
