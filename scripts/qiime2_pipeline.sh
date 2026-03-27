pwd

conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
conda activate qiime2-amplicon-2024.10
qiime --version
qiime --help
conda deactivate

python qiime_manifest.py
head -5 qiime_manifest_v2.tsv
cat -vet qiime_manifest_v2.tsv | head

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest_v2.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

qiime tools view demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 0 \
  --p-n-reads-learn 250000 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime tools view stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv

qiime feature-table filter-samples \
  --i-table table.qza \
  --p-min-frequency 1000 \
  --o-filtered-table table-filtered.qza

qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv

qiime tools view table-filtered.qzv

qiime tools export \
  --input-path table-filtered.qza \
  --output-path exported-table

biom convert \
  -i exported-table/feature-table.biom \
  -o exported-table/feature-table.tsv \
  --to-tsv

head -1 exported-table/feature-table.tsv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime metadata tabulate \
  --m-input-file metadata_clean.tsv \
  --o-visualization metadata_check.qzv

qiime tools view metadata_check.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 18000 \
  --m-metadata-file metadata_clean.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata_clean.tsv \
  --o-visualization shannon_groups.qzv

qiime tools view shannon_groups.qzv

qiime tools export \
  --input-path core-metrics-results/bray_curtis_distance_matrix.qza \
  --output-path exported-bray

conda env create \
  -n qiime2-2023.9 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2023.9-py38-osx-conda.yml
conda activate qiime2-2023.9

qiime --version

wget https://data.qiime2.org/2023.9/common/silva-138-99-515-806-nb-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 2

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata_clean.tsv \
  --o-visualization taxa-barplot.qzv

qiime tools view taxa-barplot.qzv

qiime tools export \
  --input-path taxonomy.qza \
  --output-path taxonomy_export

conda deactivate
conda activate qiime2-amplicon-2024.10

qiime tools export \
  --input-path table.qza \
  --output-path table_export

biom convert \
  -i table_export/feature-table.biom \
  -o table_export/feature-table.tsv \
  --to-tsv

