#!/bin/bash

# Opens a DNAnexus cloud workstation, and processes the specified 100K EGLH
# cases using the caerus.py pipeline.

# Example usage:
#   Process all cases:      caerus.sh all
#   Process single case:    caerus.sh <case_id>
#   Process multiple cases: caerus.sh <case_id_1> <case_id_2> ... <case_id_3>

# DNAnexus folder structure in 003_230124_caerus:
#   caerus.sh
# 	caerus.py
# 	jsons
#       <CASE_ID>.json
#       ...
#   vcfs
#       <CASE_ID>_proband.vcf
#       <CASE_ID>_<relation to proband>.vcf
#       ...

dx select 003_230124_caerus
dx run --ssh app-cloud_workstation
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

if [[ $# == 1 ]]; then  # if there's only one arg,
    if [[ "$1" == "all" ]]; then  # process all cases we have a json for
        dx download jsons/*.json -f
        for json in *.json; do
            case_id="${json}%.json"
            dx download vcfs/"$case_id"*.vcf*
            vcfs=""
            for vcf in "$case_id"*.vcf*; do
                vcfs="${vcfs} ${vcf}"
            done
            python caerus.py "$case_id" "$json" "$vcfs"
        done

    else  # assume single arg is a case id
        case_id="$1"
        json="${case_id}.json"
        dx download jsons/"$json" -f
        dx download vcfs/"$case_id"*.vcf*
        vcfs=""
        for vcf in "$case_id"*.vcf*; do
            vcfs="${vcfs} ${vcf}"
        done
        python caerus.py "$case_id" "$json" "$vcfs"
    fi

elif [[ $# -gt 1 ]]; then  # if there are multiple args, process multiple cases
    for case in "$@"; do
        case_id="$case"
        json="${case_id}.json"
        dx download vcfs/"$case_id"*.vcf*
        vcfs=""
        for vcf in "$case_id"*.vcf*; do
            vcfs="${vcfs} ${vcf}"
        python caerus.py "$case_id" "$json" "$vcfs"
        done
    done

else  # no args, no service
info_1="Specify cases to process as a list of space-separated case IDs."
info_2+="Example: caerus.sh CSA-XXXX SAP-YYYY"
info_3="Alternatively, to process all available cases, use 'caerus.sh all'."
echo "$info_1"
echo "$info_2"
echo "$info_3"

fi


# add a sanity check: are the static files present before we start the pipeline
# 	biomart api call template
# 	chromosome notation map
# 	genome file
# 	chain file
