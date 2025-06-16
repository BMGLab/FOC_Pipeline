// 
process FILTER_CORE_GENES {
    tag "Filtering core genes from: ${gffs.size()} samples"

    input:
    path gffs
    path faas

    output:
    path "core_filtered/*", emit: out_core_filtered

    script:
    def faa_flag = faas ? "--faa_dir faas" : ""
    """
    mkdir -p core_filtered
    python filter_core_genes_augustus_liftoff.py \
    --gff_dir $gffs \
    $faa_flag \
    --output_dir core_filtered
    """
}