clear && nextflow \
    -log logs/.nextflow.log \
    run main.nf \
    --skip_deeptmhmm true \
    -resume
