wget -O GSE51032_RAW.tar "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51032&format=file"
tar xvf GSE51032_RAW.tar
for f in *.idat.gz; do gunzip $f; echo $f; done

Rscript minfi-norm.R

# you can move the file where you like::
#mv norm.beta.txt /drive/450k/
#gzip -9 /drive/450k/norm.beta.txt
