# Create database

```
# Initial database setup
cat build_webstr_db.sql | sqlite3 webstr2.db
```

```
# Download gene annotations
curl -o gene_annots/hg19.gtf.gz http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
curl -o gene_annots/hg38.gtf.gz http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
curl -o gene_annots/mm10.gtf.gz http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
curl -o gene_annots/rn7.gtf.gz http://ftp.ensembl.org/pub/release-107/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.107.gtf.gz

# Add gene annotations
for build in hg19 hg38 mm10 rn7
do
	echo "start $build"
	./add_annot.py $build gene_annots/${build}.gtf.gz webstr2.db
	echo "Completed $build"
done
```