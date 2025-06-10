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
