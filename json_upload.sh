#!/bin/bash

# # script to upload jsons from cyto drive to 003_230124_caerus/jsons via putty

source /mnt/storage/apps/software/dnanexus/0.306.0/dx-toolkit/environment
cd /mnt/storage/samba/samba.ctrulab.uk/cytogenetics/staging_area/jay/jsons

dx login
dx select 003_230124_caerus
dx cd jsons

for old_json in *.json; do
    rd_check=$(cat "$old_json" | grep '"program": "rare_disease"')
    if [ ${rd_check:+1} ]; then
        case_id=${old_json%__irId=*}
        case_id=${case_id##*InterpretationDetail_caseID_Page}
        case_id=${case_id##*_}
        new_json="${case_id}.json"
        dx upload --brief --path "${new_json}" "$old_json"
    fi
done
