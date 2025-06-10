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

