#!/bin/bash

# Given a case id for an EGLH 100k rare disease case, download all files
# from that case's DNAnexus folder (if they aren't an SV VCF or already
# downloaded.).

dx select 003_230124_caerus

# check only 1 argument was given

if [[ $# == 1 ]]; then

    case_id="$1"

    # identify & download json file for case

    json_file=$(dx ls jsons/*"${case_id}"*)

    if [[ ! -f "$json_file" ]]; then
        echo "Downloading: ${json_file}"
        dx download "jsons/${json_file}" -f --no-progress

    else
        echo "Already downloaded: ${json_file}"

    fi

    # identify and download SNV VCFs for case

    for file in $(dx ls "$case_id"); do
        if [[ "$file" != *.SV.* ]]; then

            if [[ ! -f "$file" ]]; then
                echo "Downloading: ${file}"
                dx download "${case_id}/${file}" -f --no-progress

            else
                echo "Already downloaded: ${file}"

            fi
        fi
    done

else
    info_1="file_downloader.sh requires exactly one argument. "
    info_2="Example: 'file_downloader.sh SAP-35035-1'. "
    info_3="Arguments passed: $*"
    echo "${info_1}${info_2}${info_3}"

fi
