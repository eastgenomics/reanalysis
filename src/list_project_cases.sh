#!/usr/bin/bash

# Identify cases in the 003_230124_caerus project with both a JSON file
# and at least one VCF. Print these case IDs to lines of a text file.

dx select 003_230124_caerus

cases=()

for folder in $(dx ls); do
    has_json=False
    has_vcf=False

    for file in $(dx ls "$folder"); do
        if [[ "$file" == *.json ]]; then
            has_json=True
        fi
        if [[ "$file" == *.vcf.gz ]]; then
            has_vcf=True
        fi
    done

    if $has_json and $has_vcf; then
        cases+=("$folder")
    fi
done

printf "%s\n" "${cases[@]}" > usable_project_cases.txt
