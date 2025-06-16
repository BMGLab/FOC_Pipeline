#!/bin/bash

set -e

# === CONFIGURATION ===
SPECIES_NAME="fol"
GENOME="../ref/GCF_000149955.1_ASM14995v2_genomic.fna"
ANNOTATION="../ref/GCF_000149955.1_ASM14995v2_genomic.gff"
PROTEINS="GCF_000149955.1_ASM14995v2_protein.faa.gz"
WORKDIR="augustus_training_fol"
FLANKING=1000

# === ENVIRONMENT SETUP ===
echo "[*] Activating augustus environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate augustus

# Copy AUGUSTUS config to writable local directory if needed
if [ -z "$AUGUSTUS_CONFIG_PATH" ]; then
    export AUGUSTUS_CONFIG_PATH="$HOME/augustus_config"
    if [ ! -d "$AUGUSTUS_CONFIG_PATH" ]; then
        echo "[*] Creating writable AUGUSTUS config at $AUGUSTUS_CONFIG_PATH"
        cp -r $CONDA_PREFIX/config "$AUGUSTUS_CONFIG_PATH"
    fi
fi

# === DIRECTORY SETUP ===
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# === STEP 1: Create Species Model ===
if [ ! -d "$AUGUSTUS_CONFIG_PATH/species/$SPECIES_NAME" ]; then
    echo "[*] Creating AUGUSTUS species: $SPECIES_NAME"
    new_species.pl --species=$SPECIES_NAME
fi

# === STEP 2: Convert GFF + FASTA to GenBank Format ===
echo "[*] Converting GFF and FASTA to GenBank format..."
gff2gbSmallDNA.pl ../$ANNOTATION ../$GENOME $FLANKING fol.gb

# === STEP 3: Create Training and Test Sets ===
echo "[*] Splitting GenBank file into training and test sets..."
randomSplit.pl fol.gb 100

# === STEP 4: Training AUGUSTUS ===
echo "[*] Training AUGUSTUS model for $SPECIES_NAME..."
etraining --species=$SPECIES_NAME fol.gb.train

# === STEP 5: Testing Model ===
echo "[*] Running test prediction to validate model..."
augustus --species=$SPECIES_NAME fol.gb.test > test_output.gff

# === STEP 6: Optional - Accuracy Report ===
if command -v computeAccuracy.pl &> /dev/null; then
    echo "[*] Computing prediction accuracy..."
    computeAccuracy.pl --pred=test_output.gff --ref=fol.gb.test
fi

echo "[✓] Training complete! Model stored in: $AUGUSTUS_CONFIG_PATH/species/$SPECIES_NAME"