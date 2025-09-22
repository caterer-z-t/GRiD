#!/bin/bash
# common.sh — shared functions and validation

# -------------------------------
# USAGE function (generic template)
# -------------------------------
usage_common() {
    echo "Usage: $0 [-c CRAM_DIR] [-o OUTPUT_FILE] [-r REF_FASTA] [-C CHR] [-s START] [-e END] [-w WORK_DIR] [-h]"
    echo ""
    echo "  -c   Directory containing CRAM files"
    echo "  -o   Output file path"
    echo "  -r   Reference FASTA file"
    echo "  -C   Chromosome (integer or string: 6, X, Y, MT, with/without chr)"
    echo "  -s   Start position of the VNTR region"
    echo "  -e   End position of the VNTR region"
    echo "  -w   Working directory"
    echo "  -h   Show this help message"
    exit 1
}

# -------------------------------
# Normalize chromosome input
# -------------------------------
normalize_chr() {
    local input="$1"
    if [[ "$input" =~ ^[0-9]+$ ]]; then
        echo "chr${input}"
    elif [[ "$input" =~ ^chr[0-9]+$ ]] || [[ "$input" =~ ^chr(X|Y|MT)$ ]]; then
        echo "$input"
    elif [[ "$input" =~ ^(X|Y|MT)$ ]]; then
        echo "chr${input}"
    else
        echo "Error: Invalid chromosome format '$input'." >&2
        exit 1
    fi
}

# -------------------------------
# Ensure WORK_DIR matches output file directory
# -------------------------------
sync_workdir_with_output() {
    local output_file="$1"
    local work_dir="$2"

    local output_dir
    output_dir=$(dirname "$output_file")

    if [ "$output_dir" != "$work_dir" ]; then
        echo "$output_dir"
    else
        echo "$work_dir"
    fi
}

# -------------------------------
# Argument & Input Validation Helpers
# -------------------------------

# Ensure required arguments are provided
# Usage: require_args VAR1 VAR2 VAR3 ...
require_args() {
    local missing=false
    for var in "$@"; do
        if [ -z "${!var}" ]; then
            echo "Error: Missing required argument: $var"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        usage
        exit 1
    fi
}

# Verify that a directory exists
# Usage: require_dir "/path/to/dir"
require_dir() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' does not exist."
        exit 1
    fi
}

# Verify that a file exists
# Usage: require_file "/path/to/file"
require_file() {
    local file=$1
    if [ ! -f "$file" ]; then
        echo "Error: File '$file' not found."
        exit 1
    fi
}

# -------------------------------
# Find CRAM files in a directory
# -------------------------------
# Usage: find_cram_files "/path/to/cram_dir"
# Returns: array of cram files
find_cram_files() {
    local dir=$1
    shopt -s nullglob
    local files=("$dir"/*.cram)
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
        echo "Error: No CRAM files found in $dir" >&2
        exit 1
    fi

    echo "${files[@]}"
}
