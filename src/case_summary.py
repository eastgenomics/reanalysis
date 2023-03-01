#!usr/bin/env python

# dx login
# dx select 003_230124_caerus
# dx run --ssh app-cloud_workstation
# unset DX_WORKSPACE_ID
# dx cd $DX_PROJECT_CONTEXT_ID:
# dx download case_summary.py -f
# dx cd jsons
# dx download Int*.json -f

import os
import json

output = {
    'assembly': {
        'b37': 0,
        'b38': 0,
        'other': []},
    'panels': [],
    'solved': {
        'not_p_or_f': [],
        'p_only': [],
        'f_only': [],
        'p_and_f': []},
    'variants': {
        'all_cases': [],
        'not_p_or_f': [],
        'p_only': [],
        'f_only': [],
        'p_and_f': []}}

for file in os.listdir():
    if file.startswith('Int') and file.endswith('.json'):
        with open(file, 'r') as reader:
            contents = json.load(reader)

        # get assembly details

        if contents['assembly'] == 'GRCh38':
            output['assembly']['b38'] += 1

        elif contents['assembly'] == 'GRCh37':
            output['assembly']['b37'] += 1

        else:
            output['assembly']['other'].append(contents['assembly'])

        # get panel details

        try:
            for panel in contents['interpretation_request_data']['json_request'][
            'pedigree']['analysisPanels']:

                output['panels'].append({
                    'name': panel['specificDisease'],
                    'id': panel['panelName'],
                    'version': panel['panelVersion']})

        except KeyError as e:
            print(f"\nError with panels for case {contents['case_id']}: {e}\n")

        # get solved info

        var_count = 0
        solved_fully = False
        solved_partially = False

        var_types = [
            'variants',
            'structuralVariants',
            'shortTandemRepeats',
            'chromosomalRearrangements']

        for report in contents['clinical_report']:
            if report['valid']:
                try:
                    solved = report['exit_questionnaire'][
                        'exit_questionnaire_data'][
                            'familyLevelQuestions']['caseSolvedFamily']

                    if solved == 'yes':
                        solved_fully = True
                        solved_partially = True

                    elif solved == 'unknown':
                        solved_partially = True

                    # get number of pcvs

                    for var_type in var_types:
                        try:
                            for var in report[
                                'clinical_report_data'][var_type]:

                                var_count += 1

                        except Exception:
                            pass
                except Exception:
                    pass

        # collect info

        output['variants']['all_cases'].append(var_count)

        if (not solved_fully) and (not solved_partially):
            output['solved']['not_p_or_f'].append(contents['case_id'])
            output['variants']['not_p_or_f'].append(var_count)

        elif (not solved_fully) and solved_partially:
            output['solved']['p_only'].append(contents['case_id'])
            output['variants']['p_only'].append(var_count)

        elif solved_fully and (not solved_partially):
            output['solved']['f_only'].append(contents['case_id'])
            output['variants']['f_only'].append(var_count)

        elif solved_fully and solved_partially:
            output['solved']['p_and_f'].append(contents['case_id'])
            output['variants']['p_and_f'].append(var_count)

# create json file

output_file = "summary_data_all_cases.json"

with open(output_file, 'w') as writer:
    json.dump(output, writer)

# dx upload --brief summary_data_all_cases.json
