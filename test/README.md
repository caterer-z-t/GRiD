# Testing

given the large size of the cram and crai files we will not upload these to github, however if you wish to test things please feel free to use your `.bam`/`.cram` files. Or feel free to download these `.cram` and `.crai` files from [International Genome](https://www.internationalgenome.org/) and their [Data Portal](https://www.internationalgenome.org/data-portal). For the purpose of testing I predominatly tested on [Sample HG00152](https://www.internationalgenome.org/data-portal/sample/HG00152) the 1KG_ONT_VIENNA hg38 `.cram` and `.crai` files. to download the files I used for testing please run the following commands. 

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00152.hg38.cram
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00152.hg38.cram.crai
```

Additionally ensure you have the correct FASTA file for hg38. 
```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

to test pleast run 
```bash
pytest -v test/
```