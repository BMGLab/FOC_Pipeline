clear && nextflow \
    -log logs/.nextflow.log \
    run main.nf \
    --skip_deeptmhmm true \
    --samplesheet_path 'samplesheet copy.csv' \
    # -resume