#!/bin/bash

# Given a case id for an EGLH 100k rare disease case, download all files
# from that case's DNAnexus folder (if they aren't an SV VCF or already
# downloaded.).

dx select 003_230124_caerus

if [[ $# == 1 ]]; then

    case_id="$1"
    echo "Retrieving files for ${case_id}"

    for file in $(dx ls "$case_id"); do
        if [[ "$file" != *.SV.* ]]; then
            if [[ ! -f "$file" ]]; then
                dx download "${case_id}/${file}" -f
            fi
        fi
    done

else
info_1="file_downloader.sh requires exactly one case ID as an argument. "
info_2="Example: 'caerus.sh SAP-35035-1'. "
info_3="Arguments passed: $*"
echo "${info_1}${info_2}${info_3}"
fi

echo "File download completed."
