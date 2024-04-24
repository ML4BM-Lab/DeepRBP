curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_isoform_tpm.gz -o ./TcgaTargetGtex_rsem_isoform_tpm.gz
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz -o ./TcgaTargetGtex_rsem_gene_tpm.gz
curl https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz -o ./TcgaTargetGTEX_phenotype.txt.gz 
gzip -d ./TcgaTargetGTEX_phenotype.txt


