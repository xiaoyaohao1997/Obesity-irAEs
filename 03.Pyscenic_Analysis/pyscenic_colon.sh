source ~/.bashrc
conda activate py38
cd ~/colon/pyscenic

input_loom=merged_seurat_int.loom
n_workers=24

tfs=~/scenic/allTFs_mm.txt
feather=~/scenic/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=~/scenic/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

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
