--Create GENOMES table
CREATE TABLE GENOMES (
	genome_id int,
    genome_build varchar(10)
);
INSERT INTO GENOMES VALUES (1, 'hg19');
INSERT INTO GENOMES VALUES (2, 'hg38');
INSERT INTO GENOMES VALUES (3, 'mm10');
INSERT INTO GENOMES VALUES (4, 'rn7');

--Create table with gene annotations
CREATE TABLE GENEANNOTATIONS (
	feature_id int,
	genome_build varchar(10),
	feature_chrom varchar(10),
	feature_type varchar(20),
	feature_start int,
	feature_end int,
	feature_strand int,
	gene_name varchar(20),
	gene_id varchar(20),
	gene_type varchar(20),
	FOREIGN KEY (genome_build)
		REFERENCES GENOMES (genome_build)
);