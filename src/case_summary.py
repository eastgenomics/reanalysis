#!usr/bin/env python

# dx login
# dx select 003_230124_caerus
# dx run --ssh app-cloud_workstation
# unset DX_WORKSPACE_ID
# dx cd $DX_PROJECT_CONTEXT_ID:
# dx download case_summary.py -f
# dx download jsons/Int*.json -f

import os
import json

output = {
    'assembly': {
        'all_cases': {'37': [], '38': []},
        'not_p_or_f': {'37': [], '38': []},
        'p_only': {'37': [], '38': []},
        'f_only': {'37': [], '38': []},
        'p_and_f': {'37': [], '38': []}},
    'panels': {
        'all_cases': [],
        'not_p_or_f': [],
        'p_only': [],
        'f_only': [],
        'p_and_f': []},
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
        'p_and_f': []},
    'family': {
        'all_cases': [],
        'not_p_or_f': [],
        'p_only': [],
        'f_only': [],
        'p_and_f': []}}

for file in os.listdir():
    if file.startswith('Int') and file.endswith('.json'):
        with open(file, 'r') as reader:
            contents = json.load(reader)

        case_panels = []
        vcf_count = 0
        mother_vcf = False
        father_vcf = False
        both_parents = False

        req_data = contents.get('interpretation_request_data')
        if req_data:
            request = req_data.get('json_request')
            if request:
                pedigree = request.get('pedigree')
                if pedigree:

                    # get panel details

                    panels = pedigree.get('analysisPanels')

                    for panel in panels:
                        case_panels.append({
                            'name': panel.get('specificDisease'),
                            'id': panel.get('panelName'),
                            'version': panel.get('panelVersion')})

                    # get family vcf info

                    family = pedigree.get('members')

                    for member in family:
                        relation = None
                        has_vcf = False

                        if member.get('isProband'):
                            relation = 'proband'
                        else:
                            more_info = member.get('additionalInformation')
                            relation = more_info.get('relation_to_proband')

                        if relation:
                            relation = relation.lower()

                        samples = [s.get('sampleId') \
                            for s in member.get('samples') \
                            if s.get('product') == 'DNA']

                        if samples:
                            for vcf in request.get('vcfs'):
                                for sample in samples:
                                    if sample in vcf.get('sampleId'):
                                        has_vcf = True

                        if has_vcf:
                            vcf_count += 1

                            if relation == 'mother':
                                mother_vcf = True
                            elif relation == 'father':
                                father_vcf = True

        if mother_vcf and father_vcf:
            both_parents = True

        fam_dict = {'vcf_count': vcf_count, 'both_parents': both_parents}

        # get solved & variant info

        vars = []
        solved_fully = False
        solved_partially = False

        var_types = [
            'variants',
            'structuralVariants',
            'shortTandemRepeats',
            'chromosomalRearrangements']

        for report in contents.get('clinical_report'):
            if report.get('valid'):

                # identify whether case is fully or partially solved

                ques = report.get('exit_questionnaire')
                if ques:
                    ques_data = ques.get('exit_questionnaire_data')
                    if ques_data:
                        qs = ques_data.get('familyLevelQuestions')
                        if qs:
                            solved = qs.get('caseSolvedFamily')

                            if solved == 'yes':
                                solved_fully = True
                                solved_partially = True

                            elif solved == 'unknown':
                                solved_partially = True

                # get info on PCVs

                rep_data = report.get('clinical_report_data')

                for var_type in var_types:
                    if rep_data.get(var_type):
                        for var in rep_data[var_type]:

                                var_id = None
                                c = var.get('variantCoordinates')

                                if c:
                                    var_id = f"{c['assembly']}:{c['chromosome']}:{c['position']}:{c['reference']}:{c['alternate']}"

                                if var.get('reportEvents'):
                                    for event in var['reportEvents']:
                                        var_info = {'id': var_id,
                                            'tier': event.get('tier'),
                                            'moi': event.get('modeOfInheritance'),
                                            'penetrance': event.get('penetrance'),
                                            'segregation': event.get('segregationPattern'),
                                            'fully_explains': event.get('fullyExplainsPhenotype')}

                                        vars.append(var_info)

        # collect info

        output['variants']['all_cases'].append(vars)
        output['panels']['all_cases'].append(case_panels)
        output['family']['all_cases'].append(fam_dict)

        if (not solved_fully) and (not solved_partially):
            outcome = 'not_p_or_f'
        elif (not solved_fully) and solved_partially:
            outcome = 'p_only'
        elif solved_fully and (not solved_partially):
            outcome = 'f_only'
        elif solved_fully and solved_partially:
            outcome = 'p_and_f'

        output['solved'][outcome].append(contents['case_id'])
        output['variants'][outcome].append(vars)
        output['panels'][outcome].append(case_panels)
        output['family'][outcome].append(fam_dict)

        if contents['assembly'] == 'GRCh38':
            output['assembly']['all_cases']['38'].append(contents['case_id'])
            output['assembly'][outcome]['38'].append(contents['case_id'])
        elif contents['assembly'] == 'GRCh37':
            output['assembly']['all_cases']['37'].append(contents['case_id'])
            output['assembly'][outcome]['37'].append(contents['case_id'])

# create json file

output_file = "summary_data_all_cases.json"

with open(output_file, 'w') as writer:
    json.dump(output, writer)

# dx upload --brief summary_data_all_cases.json
