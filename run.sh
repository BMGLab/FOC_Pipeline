export NXF_VER=24.10.5
export SLACK_WEBHOOK="https://hooks.slack.com/services/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"
clear && nextflow \
    -log logs/.nextflow.log \
    run main.nf \
    --skip_deeptmhmm true \
    --samplesheet_path 'samplesheet.csv' \