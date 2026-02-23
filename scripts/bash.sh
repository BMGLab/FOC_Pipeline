conda activate fox-phylo

cd ~/Projects/FoxGenome_wd

mkdir -p beyza_article_analysis/ncbi_genome_comparison/phylogenomics_filtered

python scripts/fox_busco_phylogeny2.py \
  --keep-samples beyza_article_analysis/ncbi_genome_comparison/mash_cluster/keep.txt \
  --out-dir beyza_article_analysis/ncbi_genome_comparison/phylogenomics_filtered \
  --tree-programs iqtree3 raxml-ng \
  --raxml-model partition \
  --tree-threads 16 \
  --iqtree-bootstrap 1000


