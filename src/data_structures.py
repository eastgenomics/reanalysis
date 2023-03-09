#!usr/bin/env python

import os
import re
import gzip
import json
import statistics
import subprocess

from datetime import datetime as dt
from panelapp import api, Panelapp
from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpencgaClient
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient


def generate_date():
    """ Returns current date as string in the form YYYYMMDD """

    current_date = dt.today()
    date = str(dt.strftime(current_date, "%Y%m%d"))

    return date


def session_setup(user, host):
    """ Define username and host server to use for OpenCGA login. Log in
    generates password prompt in terminal. Session token then generated and
    used to define OpenCGA client instance for the program.

    args:
        user [string]: username for OpenCGA login
        host [string]: host server URL

    returns:
        oc: OpenCGA client instance based on user credentials and token
    """

    print(f'{dt.now()} Setting up OpenCGA client')

    config_dict = {'rest': {'host': host}}
    config = ClientConfiguration(config_dict)

    oc = OpencgaClient(config)
    oc.login(user)

    token = oc.token
    oc = OpencgaClient(configuration=config, token=token)

    return oc


def get_all_project_ids(oc):
    """ Print and return the ids of all OpenCGA projects to which user
    has access.

    args:
        oc: OpenCGA client

    returns:
        project_ids [list]
    """

    all_projects = oc.projects.search().get_results()

    project_ids = [project['id'] for project in all_projects]
    project_ids.sort()

    # for id in project_ids:
    #     print(id)

    return project_ids


def get_project_study_ids(oc, project_id):
    """ Print and return the ids of all studies in a single project.

    args:
        oc: OpenCGA client
        project_id [str]

    returns:
        study_ids [list]
        """

    project_studies = oc.studies.search(project=project_id).get_results()

    study_ids = [study['id'] for study in project_studies]
    study_ids.sort()

    # for id in study_ids:
    #     print(id)

    return study_ids


def get_study_case_ids(oc, study_id, output_fp):
    """ Return and write out the ids of all cases in a single study.

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        case_ids [list]: list of case ID strings
    """

    cases = oc.clinical.search(study=study_id, include='id').get_results()

    case_ids = [case['id'] for case in cases]
    case_ids.sort()

    # print(f'Study {study_id} contains {len(case_ids)} cases')

    with open(output_fp, 'w') as writer:
        for id in case_ids:

            writer.write(f'{id}\n')

    return case_ids


def list_all_cases(oc, folder):
    """ Identify:

    - all projects to which user has access
    - all studies within those projects
    - all cases within those studies """

    all_cases = []

    project_ids = get_all_project_ids(oc)

    for project in project_ids:

        studies = get_project_study_ids(oc, project)

        for study in studies:

            study_id = f'emee-glh@{project}:{study}'
            output_fp = f'{folder}{project}_{study}_all_case_ids.txt'

            cases = get_study_case_ids(oc, study_id, output_fp)

            if study in ['rd37', 'rd38']:
                all_cases += [case for case in cases if case not in all_cases]

    all_cases.sort()

    with open(f'{folder}all_b37_b38_cases.txt', 'w') as writer:
        for case in all_cases:
            writer.write(f'{case}\n')


def get_single_case_data(oc, study_id, case_id):
    """ Retrieve all data from specified case in specified study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
        case_id [str]: id specifying case within study

    returns:
        case_data [dict]: all data associated with that case
    """

    print(f'{dt.now()} Retrieving case data')

    case_data = oc.clinical.info(
        study=study_id,
        clinical_analysis=case_id
        ).get_results()[0]

    return case_data


def look_at_json(fp):
    """ Examine the structure and contents of a JSON file.

    args:
        fp [str]: path to json file
    """

    types = [str, int, float, bool]

    with open(fp, 'r') as reader:
        data = json.load(reader)

    info = data

    # print(type(info))
    # print(len(info))
    # print(type(info[0]))
    # print(f'{key}: {value} [{type(value)}]')
    # print(f'{key} ({len(value)})')
    # print(f'{key} [{type(value)}]')

    # # IF YOU JUST WANT TO PRINT THE WHOLE THING
    # print(info)

    # IF YOU'RE LOOKING AT A DICT
    for key, value in info.items():

    # # IF YOU'RE LOOKING AT A LIST OF DICTS
    # for dict in info:
    #     for key, value in dict.items():

        # indent for lists of dicts, deindent for dicts

        if type(value) in types:
            print(f'{key}: {value} [{type(value)}]')

        elif type(value) == dict:
            if len(value) == 0:
                print(f'{key} [empty dict]')
            else:
                print(f'{key} [{type(value)}]')

        elif type(value) == list:
            if len(value) == 0:
                print(f'{key} [empty list]')
            else:
                print(f'{key} [list of {len(value)} {type(value[0])}]')

        else:
            print(f'{key} [{type(value)}]')


def look_at_100k_cases(cases):

    relationships = []
    affected_statuses = []
    adopted_statuses = []

    for case in cases:

        with open(case, 'r') as reader:
            data = json.load(reader)

        for person in data['interpretationRequestData']['pedigree']['members']:

            relation = person['additionalInformation']['relation_to_proband']
            affected = person['affectionStatus']
            adopted = person['adoptedStatus']

            if relation not in relationships:
                relationships.append(relation)

            if affected not in affected_statuses:
                affected_statuses.append(affected)

            if adopted not in adopted_statuses:
                adopted_statuses.append(adopted)

    print(relationships)
    print(affected_statuses)
    print(adopted_statuses)


def look_at_vcf(vcf):
    """ Get the number of lines in the header and variant sections of a
    VCF. Identify all unique values in the CHROM field.

    args:
        vcf [str]

    returns:
        head_count [int]: number of header lines
        body_count [int]: number of variant (body) lines
        chroms [list]: unique values in the CHROM field
    """

    head_count = count_vcf_section_lines(vcf, 'header')
    body_count = count_vcf_section_lines(vcf, 'variant')

    print(f"{vcf} has {head_count} header lines and {body_count} body lines")

    # head_lines = get_vcf_section_lines(vcf, 'header')
    body_lines = get_vcf_section_lines(vcf, 'variant')

    chroms = []

    for line in body_lines:
        chrom = re.search(r'^(.*?)\t', line)

        if chrom not in chroms:
            chroms.append(chrom)

    chroms.sort()

    print(f"{vcf} has {len(chroms)} unique CHROM values: {chroms}")

    return head_count, body_count, chroms


def get_vcf_section_lines(vcf, section):
    """ Given the path to a VCF file, return the header or variant lines
    as a list.

    args:
        vcf [str]
        section [str]: 'header' or 'variant'

    returns:
        lines [list]
    """

    if section == 'header':

        content = subprocess.run(f"bcftools view -h {vcf}",
            shell=True, capture_output=True, text=True).stdout

    elif section == 'variant':

        content = subprocess.run(f"bcftools view -H {vcf}",
            shell=True, capture_output=True, text=True).stdout

    lines = [line.strip() for line in content.split('\n') if line.strip()]

    return lines


def count_vcf_section_lines(vcf, section):
    """ Given the path to a VCF file, return the number of header or
    body lines.

    args:
        vcf [str]
        section [str]: 'header' or 'variant'

    returns:
        count [int]
    """

    if section == 'header':

        content = subprocess.run(f"bcftools view -h {vcf} | wc -l",
            shell=True, capture_output=True, text=True).stdout

    elif section == 'variant':

        content = subprocess.run(f"bcftools view -H {vcf} | wc -l",
            shell=True, capture_output=True, text=True).stdout

    count = int(content)
    print(f"{vcf} has {count} {section} lines")

    return count


def cb_client_annotation(input, output_json):
    """ Use the CB client to retrieve variant information """

    variants = input  # call using chrom:pos:ref:alt for single variant
    # variants = create_variant_list(input)

    cc = ConfigClient("resources/cb_config.json")
    cbc = CellBaseClient(cc)
    var_client = cbc.get_variant_client()

    var_info = var_client.get_annotation(variants, include=[
        'conservation', 'geneConstraints', 'functionalScore',
        'consequenceType', 'populationFrequencies'])

        # options for 'include':

        # 'hgvs', 'cytoband', 'conservation', 'geneConstraints',
        # 'functionalScore', 'variation' (='id'),
        # 'mirnaTargets' (='geneMirnaTargets'),
        # 'geneDisease' (='geneTraitAssociation')
        # 'drugInteraction' (='geneDrugInteraction'),
        # 'populationFrequencies' (= 'populationFrequencies' + 'id')
        # 'consequenceType' (= 'consequenceTypes' + 'displayConsequenceType')

    with open(output_json, 'w') as writer:
        json.dump(var_info, writer)

    return var_info


def create_variant_list(input_vcf):
    """ Given the path to a VCF file, read in the file and parse out a
    list of SNVs in the form 'chrom:pos:ref:alt'

    args:
        fp [str]

    returns:
        variants [list]
    """

    with open(input_vcf, 'r') as reader:
        lines = reader.readlines()

    variants = []

    for line in lines:
        if line[0] != '#':  # ignore header lines

            line_data = [value.strip() for value in line.split('\t')]

            chrom = line_data[0]
            pos = line_data[1]
            ref = line_data[3]
            alt = line_data[4]

            if (len(ref) == 1) and (len(alt) == 1):  # only look at SNVs

                string = f'{chrom}:{pos}:{ref}:{alt}'
                variants.append(string)

    return variants


def condense_annotation(variants, output_fp=None):
    """ Extract specific data values from an annotated list of variants.

    args:
        variants [list]: dicts of variant annotations
        output_fp [fp or None]: if supplied, dump output to json
    """

    data = []

    for variant in variants:

        var_dict = {
            'variant': None,
            'consequence': None,
            'allele_frequency': None,
            'scaled_cadd': None,
            'affected_transcripts': {},
            'constraint': {},
            'VDA': []}

        results = variant['results'][0]

        # get the variant id
        var_dict['variant'] = variant['id']

        # get the general variant consequence
        try:
            var_dict['consequence'] = results['displayConsequenceType']

        except KeyError:
            pass

        # get allele frequency
        try:
            for dct in results['populationFrequencies']:
                if dct['population'] == 'ALL':
                    var_dict['allele_frequency'] = dct['altAlleleFreq']

        except KeyError:
            pass

        # get the scaled CADD score
        try:
            for dct in results['functionalScore']:
                if dct['source'] == 'cadd_scaled':
                    var_dict['scaled_cadd'] = dct['score']

        except KeyError:
            pass

        # get gene constraint data
        try:
            for dct in results['geneConstraints']:
                if dct['name'] != 'exac_oe_lof':
                    var_dict['constraint'][dct['name']] = dct['value']

        except KeyError:
            pass

        # get affected genes and transcripts
        try:

            genes = {}

            for dct in results['consequenceTypes']:

                gene = dct['geneName']

                if gene not in genes.keys():
                    genes[gene] = {}

                transcript = dct['transcriptId']

                if transcript not in genes[gene].keys():
                    genes[gene][transcript] = []

                for effect in dct['sequenceOntologyTerms']:
                    if effect['name'] not in genes[gene][transcript]:
                        genes[gene][transcript].append(effect['name'])

            var_dict['affected_transcripts'] = genes

        except KeyError:
            pass

        # get DisGeNET variant-disease associations
        try:
            for ele in variant['geneTraitAssociation']:

                var_dict['VDA'].append({
                    'Trait': ele['name'],
                    'Score': ele['score']})

        except KeyError:
            pass

        data.append(var_dict)

    if output_fp:
        with open(output_fp, 'w') as writer:
            json.dump(data, writer)

    else:
        print(data)

    return var_dict


def annotation_filters(variants, parameters):
    """ Given annotated variant data and a set of filtering parameters,
    identify which variants meet filtering thresholds.

    args:
        variants [list]: 1 dict per variant, with selected annotation

        var_dict = {
            'variant': None,
            'consequence': None,
            'allele_frequency': None,
            'scaled_cadd': None,
            'affected_transcripts': {},
            'constraint': {},
            'VDA': []}

        parameters = {
            max_af,
            min_cadd,
            min_pli,
            max_oe_lof,
            max_oe_mis,
            max_oe_syn}

    returns:
        output [list]: 1 dict per variant
    """

    lof_vars = ['stop_gained', 'start_lost']

    for variant in variants:

        variant['filters_passed'] = 0

        filter_dct = {
            'af': False,
            'cadd': False,
            'constraint': False}

        # compare allele frequency to maximum

        if variant['allele_frequency'] <= parameters['max_af']:

            filter_dct['af'] = True
            variant['filters_passed'] += 1

        # compare CADD score to minimum

        if variant['scaled_cadd'] >= parameters['min_cadd']:

            filter_dct['cadd'] = True
            variant['filters_passed'] += 1

        # compare constraint to relevant values

        if (variant['consequence'] == 'synonymous_variant') and \
            (variant['constraint']['oe_syn'] <= parameters['max_oe_syn']):

            filter_dct['constraint'] = True
            variant['filters_passed'] += 1

        elif (variant['consequence'] == 'missense_variant') and \
            (variant['constraint']['oe_mis'] <= parameters['max_oe_mis']):

            filter_dct['constraint'] = True
            variant['filters_passed'] += 1

        elif (variant['consequence'] in lof_vars) and \
            ((variant['constraint']['oe_lof'] <= parameters['max_oe_lof']) or
            (variant['constraint']['exac_pLI'] >= parameters['min_pli'])):

            filter_dct['constraint'] = True
            variant['filters_passed'] += 1

    return variants


def get_cases_summary():

    # source /mnt/storage/apps/software/dnanexus/0.306.0/dx-toolkit/environment
    # cd /mnt/storage/samba/samba.ctrulab.uk/cytogenetics/staging_area/jay
    # ml python3.6.10

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

    path = "/mnt/storage/samba/samba.ctrulab.uk/cytogenetics/CIPAPI_JSON_FILES/C/"

    for file in os.listdir(path):
        if file.endswith('.json'):
            with open(f"{path}{file}", 'r') as reader:
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
                print(f"Error with panels for case {contents['case_id']}: {e}")

            # get solved info

            var_count = 0

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
                    except KeyError:
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


def process_cases_summary(file):

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

    with open(file, 'r') as reader:
        dict = json.load(reader)

    # print(f"There are {dict['assembly']['b37']} b37 and {dict['assembly']['b38']} b38 cases")
    # print(f"Overall, the mean number of potential causal variants per case is {statistics.mean(dict['variants']['all_cases'])}")

    # for key, case_list in dict['solved'].items():

    #     cases = dict['solved'][key]
    #     vars = dict['variants'][key]
    #     filename = f"case_list_{key}.txt"

    #     if vars:
    #         print(f"There are {len(cases)} {key} cases. {len(vars)} have potential causal variants. The average number of variants is {statistics.mean(vars)}.")
    #     else:
    #         print(f"There are {len(cases)} {key} cases. {len(vars)} have potential causal variants.")

    #     with open(filename, 'w') as writer:
    #         for case in case_list:
    #             writer.write(f"{case}\n")

    panels_count = {}

    for panel in dict['panels']:
        id = panel['id']
        name = panel['name']
        version = panel['version']

        # panel id in output dict
        if id in panels_count.keys():
            panels_count[id]['count'] += 1

            # panel version in id dict
            if version in panels_count[id].keys():
                panels_count[id][version]['count'] += 1

                # panel name in version dict
                if name in panels_count[id][version]['conditions'].keys():
                    panels_count[id][version]['conditions'][name] += 1

                # panel name not in version dict
                elif name not in panels_count[id][version]['conditions'].keys():
                    panels_count[id][version]['conditions'][name] = 1

            # panel version not in id dict
            elif version not in panels_count[id].keys():
                panels_count[id][version] = {'count': 1, 'conditions': {name: 1}}

        # panel id not in output dict
        elif id not in panels_count.keys():
            panels_count[id] = {'count': 1, version: {'count': 1, 'conditions': {name: 1}}}

    with open('panels_summary.txt', 'w') as writer:
        for panel_id, id_value in panels_count.items():
            writer.write(f"Panel {panel_id} ({id_value['count']} uses)\n")
            for key, value in id_value.items():
                if key != 'count':
                    writer.write(f"\tVersion {key} ({value['count']} uses)\n")
                    for cond, count in value['conditions'].items():
                        writer.write(f"\t\t{cond} ({count} uses)\n")
            writer.write("\n")


def get_panelapp_panel(panel_id, panel_version=None):
    """ Query the PanelApp API to retrieve data for the specified
    version of the specified panel. If version is not specified, returns
    the panel's current version.

    args:
        panel_id [str]
        panel_version [str/None]: retrieves current version if None

    returns:
        result [dict]: all panel data
    """

    try:
        result = Panelapp.Panel(str(panel_id), panel_version).get_data()

    # doesn't work for some older panel versions, for various reasons

    except Exception:
        path = ["panels", str(panel_id)]
        param = {"version": panel_version}
        url = api.build_url(path, param)
        result = api.get_panelapp_response(url)

    return result


def get_case_gene_list(case_id, json_fp):

    panels = []
    genes = []

    with open(json_fp, 'r') as reader:
        json_data = json.load(reader)

    json_panels = json_data['interpretation_request_data']['json_request']['pedigree']['analysisPanels']

    for panel in json_panels:
        panels.append({
            'name': panel['specificDisease'],
            'id': panel['panelName'],
            'version': panel['panelVersion']})

    assert len(panels) >= 1

    for panel in panels:
        panel_data = get_panelapp_panel(panel['id'])

        if panel_data:

            for gene in panel_data['genes']:
                if gene['confidence_level'] == '3' and \
                    gene['gene_data']['hgnc_id'] and \
                    (gene['gene_data']['hgnc_id'] not in genes):
                    genes.append(gene['gene_data']['hgnc_id'])

    print(f'{case_id} covers {len(genes)} genes')

    return genes


def test_vcf_regex(vcf):

    print('Getting variant line contents with bcftools...')

    var_list = []
    lines = get_vcf_section_lines(vcf, 'variant')

    print('Processing variant coordinates using regex...')

    for line in lines:

        c_p = re.search(r'^(.*?)\t(.*?)\t', line).groups()  # chrom & pos
        r_a = re.search(r'\t([ACGT]*)\t([ACGT]*)\t', line).groups()  # ref & alt

        for strings in c_p, r_a:
            assert len(strings) == 2, \
                f"variant regex: should be (field1, field2) but got {strings}"

        variant = f"{c_p[0]}:{c_p[1]}:{r_a[0]}:{r_a[1]}"

        var_list.append(variant)
        print(variant)

    print(f'regex found {len(var_list)} variants in {vcf}')


def get_all_refalt_chars(vcf):
    """ Identify all unique characters used in the REF or ALT fields of
    a VCF """

    chars = []

    var_lines = get_vcf_section_lines(vcf, 'variant')

    for line in var_lines:
        split = line.split('\t')
        ref = split[3]
        alt = split[4]

        for value in ref, alt:
            for char in value:
                if char not in chars:
                    chars.append(char)

    print(f'Unique REF/ALT characters in {vcf}: {chars}')


def main():
    """ Reference data """

    date = generate_date()

    res_dir = "resources/home/dnanexus/"
    bed_template = f"{res_dir}bed_template.txt"
    chrom_map = f"{res_dir}chrom_map.txt"
    genome_file = f"{res_dir}GCF_000001405.40_GRCh38.p14_genomic.fna"
    chain_file = f"{res_dir}hg19ToHg38.over.chain.gz"

    """ define filtering parameters """

    parameters = {
        'qual': 20,    # for initial QUAL value filtering
        'max_af': 0.05,  # seen in 5% of the population at most
        'min_cadd': 10,  # CADD >= 10 implies p(var not observed) >= 0.9
        'max_pli': 0.9,  # higher pLI implies less tolerant to truncation
        'max_oe_lof': 0.35,  # lower o/e implies gene less tolerant to variants
        'max_oe_mis': 0.35,
        'max_oe_syn': 0.35}

    """ function calls """

    json_input_1 = 'InterpretationDetail_caseID_Page9_SAP-32023-1__irId=32023__irVersion1_.json'
    json_output_1 = 'SAP-32023-1_temp_output.json'
    json_anno_1 = 'SAP-32023-1_annotation_group_4.json'
    vcf_1a = 'SAP-32023-1_1_sorted_proband.vcf.gz'
    vcf_1b = 'SAP-32023-1_2_normalised_proband.vcf.gz'
    vcf_1c = 'SAP-32023-1_3_renamed_proband.vcf.gz'
    vcf_1d = 'SAP-32023-1_4_tagged_proband.vcf.gz'
    vcf_1e = 'SAP-32023-1_6_filtered_proband.recode.vcf'
    vcf_1f = 'SAP-32023-1_7_sorted_proband.vcf.gz'

    json_input_2 = 'InterpretationDetail_caseID_Page9_SAP-33169-1__irId=33169__irVersion1_.json'
    json_output_2 = 'SAP-33169-1_temp_output.json'

    vcfs = [vcf_1f, vcf_1e, vcf_1d, vcf_1c, vcf_1b, vcf_1a]

    # get_variant_lines()
    # look_at_vcf(vcf)

    look_at_json(json_anno_1)

    # look_at_100k_cases(case_list)
    # process_cases_summary('summary_data_all_cases.json')

    # cb_client_annotation('7:117199533:G:A', output_json)

    # genes = get_case_gene_list('SAP-32023-1', json_1)

    # for vcf in vcfs:
    #     print(f"\nProcessing {vcf}\n")
    #     print('Counting variant lines with bcftools...')
    #     count_vcf_section_lines(vcf, 'variant')
    #     test_vcf_regex(vcf)

    # test_vcf_regex(vcf_1f)
    # get_all_refalt_chars(vcf_1f)


if __name__ == '__main__':
    main()
