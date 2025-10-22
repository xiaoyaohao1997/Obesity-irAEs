source ~/.bashrc
conda activate py38
cd ~/EGAS00001006762
echo "EGAS00001006762"

input_loom=BALF_Myeloid.loom
n_workers=16

tfs=~/scenic/allTFs_hg38.txt
feather=~/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=~/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# grn
pyscenic grn \
--num_workers $n_workers \
--output grn.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
pyscenic ctx \
grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
pyscenic aucell \
$input_loom \
ctx.csv \
--output aucell.loom \
--num_workers $n_workers
