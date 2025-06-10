#!/bin/bash

# Code for database curation


# Gaby and Buckley Database ###########################

# Rename FASTA headers to only include NCBI accessions 
awk '{if(substr($0,1,1)==">") {split($0,a," "); print "" a[1]} else print $0}' nifH_database_2012.fasta > nifH_seqs_clean.fasta

# Obtain a list of NCBI accession IDs from the database file
grep "^>" nifH_seqs_clean.fasta | sed 's/^>//' | awk '{print $1}' | sort -u > accession_ids

# Access the Enterez NCBI Database to obtain the taxonomic lineage of each corresponding accession ID 
for i in `cat accession_ids`; do printf ${i}"\t"; esearch -db nuccore -query ${i} | elink -target taxonomy | efetch -format xml | \
xtract -pattern Taxon -element ScientificName,Lineage; done > accession_taxonomy

# Reformat output file with the appropriate taxonomic syntax
cat accession_taxonomy | awk -F'\t' '{
    split($3, lineage, "; ")
    printf "%s\t", $1
    if (length(lineage) >= 2) printf "d__%s;", lineage[2]
    if (length(lineage) >= 3) printf "k__%s;", lineage[3]
    if (length(lineage) >= 4) printf "p__%s;", lineage[4]
    if (length(lineage) >= 5) printf "c__%s;", lineage[5]
    if (length(lineage) >= 6) printf "o__%s;", lineage[6]
    if (length(lineage) >= 7) printf "f__%s;", lineage[7]
    if (length(lineage) >= 8) printf "g__%s;", lineage[8]
    printf "s__%s\n", $2
}' > taxonomy_final

# Process and reformat taxonomy
#!/bin/bash

# Set the replacement lineage
replacement="d__Bacteria;k__unclassified Bacteria;p__unclassified Bacteria;c__unclassified Bacteria;o__unclassified Bacteria;f__unclassified Bacteria;g__unclassified Bacteria;s__unclassified Bacteria"
unclassified="d__unclassified;k__unclassified;p__unclassified;c__unclassified;o__unclassified;f__unclassified;g__unclassified;s__unclassified"



# Process input file
awk -v replacement="$replacement" -v unclassified="$unclassified" -F'\t' -v OFS='\t' '
{
    if ($2 ~ /k__environmental samples/) {
        $2 = replacement
    }
    else if ($2 ~ /k__unclassified Bacteria/) {
        $2 = replacement
    }
    else if ($2 !~ /d__Bacteria/) {
        $2 = unclassified
    }
    print
}
' "$1"

#!/bin/bash

awk -F'\t' -v OFS='\t' '
{
    lineage = $2
    n = split(lineage, fields, ";")

    # Initialize taxonomy map
    delete tax
    for (i = 1; i <= n; i++) {
        split(fields[i], parts, "__")
        prefix = parts[1]
        value = parts[2]
        tax[prefix] = value
    }

    # Define taxonomic levels
    level_order[1] = "d"
    level_order[2] = "k"
    level_order[3] = "p"
    level_order[4] = "c"
    level_order[5] = "o"
    level_order[6] = "f"
    level_order[7] = "g"
    level_order[8] = "s"

    found = 0
    last_taxon = ""

    new_lineage = ""

    for (i = 1; i <= 8; i++) {
        prefix = level_order[i]
        value = tax[prefix]

        if (!found) {
            if (value ~ /environmental samples/) {
                found = 1
                # start replacing from this level onward
            } else if (value != "" && value !~ /uncultured/ && value !~ /unclassified/) {
                last_taxon = value
            }
        }

        if (found) {
            value = "unclassified " last_taxon
        } else if (value == "") {
            value = ""
        }

        # Append to lineage string
        if (i == 1) {
            new_lineage = prefix "__" value
        } else {
            new_lineage = new_lineage ";" prefix "__" value
        }
    }

    $2 = new_lineage
    print
}
' "$1"

# NFIX DB ##########
####################
# Combine all .fna files 
cat *.fna > NFixDB_nifH.fasta

# Fix headers in the new fasta file
sed -E 's/(>[^# ]*).*/\1/' NFixDB_nifH.fasta > NFixDB_nifH_headers.fasta

# Filter taxonomy file to only include relevant genes
awk '/^>/ {sub(/^>/,""); split($0,a," "); print a[1]}' NFixDB_nifH_headers.fasta > headers.txt

awk '
BEGIN {
    while (getline line < "headers.txt") {
        fasta[line] = 1
    }
}
NR == 1 { next }  # Skip header row

{
    match_found = ""
    for (h in fasta) {
        if (index($0, h) > 0) {
            match_found = h
            break
        }
    }
    if (match_found != "") {
        print match_found "\t" $3
    }
}
' filteredhits_i2-1.tsv > NFixDB_nifH_taxonomy.tsv

# Remove any sequences NOT in the taxonomy file. 
awk -v tax_file="NFixDB_nifH_taxonomy.tsv" '
BEGIN {
    FS = "\t"
    while ((getline line < tax_file) > 0) {
        split(line, fields, FS)
        id = fields[1]
        gsub(/^ +| +$/, "", id)  # trim spaces
        tax_ids[id] = 1
    }
    close(tax_file)
}
{
    if ($0 ~ /^>/) {
        header = substr($0, 2)
        gsub(/^ +| +$/, "", header)
        keep = (header in tax_ids)
    }
    if (keep) print
}
' NFixDB_nifH_headers.fasta > NfixDB_nifH_seqs_qiime2format.fasta

{ echo -e "Feature ID\tTaxon"; cat NFixDB_nifH_taxonomy.tsv; } > NFixDB_nifH_taxonomy_qiime2format.tsv

# Import into QIIME2
qiime tools import --input-path NFixDB_nifH_taxonomy_qiime2format.tsv --type 'FeatureData[Taxonomy]' --output-path NfixDB_nifH_taxonomy.qza
qiime tools import --input-path NfixDB_nifH_seqs_qiime2format.fasta --type 'FeatureData[Sequence]' --output-path NfixDB_nifH_seqs.qza

# Train the feature classifier 
qiime feature-classifier fit-classifier-naive-bayes  --i-reference-reads NfixDB_nifH_seqs.qza --i-reference-taxonomy NfixDB_nifH_taxonomy.qza --o-classifier NFixDB_classifier.qza

# Test classifier
qiime feature-classifier classify-sklearn --i-classifier NFixDB_classifier.qza --i-reads /mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq/nifH_sequences.qza --o-classification NFixDB_nifH_taxonomy.qza
qiime metadata tabulate --m-input-file NFixDB_nifH_taxonomy.qza  --o-visualization NFixDB_nifH_taxonomy.qzv
