import pandas as pd

expr_data_file = "/Volumes/bioinformatics/projects/fei2021/results/merged/genes.fpkm.cufflinks.allSamples.xls"
expr_data_df = pd.read_csv(expr_data_file,sep="\t", index_col=0, header=0)

cell_days = [0,11,30]
EV_days = [0,3,5,7,9,11,12,14,16,18,20,22,24,26,28,30]


