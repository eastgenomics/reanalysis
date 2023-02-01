#!usr/bin/env python

import gzip
import json
from datetime import datetime as dt

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpencgaClient
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient


def gzip_file(fp, content):
    with gzip.open(fp, 'wt') as writer:
        writer.write(content)


def ungzip_file(fp):
    with gzip.open(fp, 'rt') as reader:
        data = reader.read()


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
    """ Examine the structure of a JSON file.

    args:
        fp [str]: path to json file
    """

    types = [str, int, float, bool]

    with open(fp, 'r') as reader:
        data = json.load(reader)

    info = data['family']['members'][0]

    # print(info)
    # print(type(info))
    # print(len(info))
    # print(type(info[0]))

    # print(f'{key}: {value} [{type(value)}]')
    # print(f'{key} ({len(value)})')
    # print(f'{key} [{type(value)}]')

    # IF YOU'RE LOOKING AT A DICT

    for key, value in info.items():

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

        # if type(value) == dict:
        #     print(f"Parent: {key}")

        #     for key1, value1 in value.items():
        #         if type(value1) in types:
        #             print(f'{key1}: {value1} [{type(value1)}]')

        #         elif type(value1) == list:
        #             if len(value1) == 0:
        #                 print(f'{key1} [empty list]')
        #             else:
        #                 print(f'{key1} [list of {len(value1)} {type(value1[0])}]')

        #         else:
        #             print(f'{key1} [{type(value1)}]')

    return fp


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
    """  """

    if vcf.endswith('.gz'):
        with gzip.open(vcf, "rt") as reader:
            lines = reader.readlines()

    else:
        with open(vcf, 'r') as reader:
            lines = reader.readlines()

    header = []
    body = []
    header_lines = 0
    body_lines = 0

    chroms = []

    for line in lines:
        if line.startswith('#'):
            header.append(line)
            header_lines += 1

        else:
            body.append(line)
            body_lines += 1

            split = line.split('\t')
            chrom = split[0].strip()

            if chrom not in chroms:
                print(chrom)
                chroms.append(chrom)

    chroms.sort()
    print(f"{vcf} has {len(chroms)} unique CHROM values")

    filename = f'{vcf}_stupid_chromosomes.txt'

    with open(filename, 'w') as writer:
        writer.write(f"{vcf} has {len(chroms)} unique values in CHROM field\n")
        for chrom in chroms:
            writer.write(f"{line}\n")


def cb_client_annotation(output, input):
    """ Use the CB client to retrieve variant information """

    t_start = dt.now()
    print(f'{t_start} Retrieving variant annotation using CB client')

    variants = create_variant_list(input)

    cc = ConfigClient("config.json")
    cbc = CellBaseClient(cc)
    var_client = cbc.get_variant_client()

    # var_info = var_client.get_annotation(single_var_1)

    var_info = var_client.get_annotation(variants, include=[
        # 'conservation',
        # 'cytoband',
        'functionalScore',
        'geneConstraints',
        # 'hgvs',
        # 'variation',  # 'id'
        'populationFrequencies',  # 'populationFrequencies' and 'id'
        'consequenceType',  # 'consequenceTypes' and 'displayConsequenceType'
        # 'mirnaTargets',  # 'geneMirnaTargets'
        'geneDisease', ])  # 'geneTraitAssociation'
        # 'drugInteraction'])  # 'geneDrugInteraction'

    with open(output, 'w') as writer:
        json.dump(var_info, writer)

    t_end = dt.now()
    t_diff = t_end - t_start
    print(f'Annotation data retrieval took {t_diff}')

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


def main():
    """ Reference data """

    date = generate_date()

    top_dir = 'do_not_upload/'
    json_dir = f'{top_dir}jsons/'
    vcf_dir = f'{top_dir}vcfs/'
    gen_dir = f'{top_dir}genomes/'
    bed_dir = f'{top_dir}bed_files/'

    case_list = [
        f'{json_dir}SAP-48034-1.json',
        f'{json_dir}SAP-55997-1.json',
        f'{json_dir}SAP-56069-1.json',
        f'{json_dir}SAP-56130-1.json',
        f'{json_dir}SAP-56172-1.json',
        f'{json_dir}SAP-56251-1.json']

    all_case_ids_file = 'all_case_ids.txt'
    case_id = 'SAP-55997-1'
    case_json = f'{json_dir}{case_id}.json'
    vcf_prefix = f'{vcf_dir}{case_id}'

    genome_file = f'{gen_dir}b38_genome.fa.bgz'
    chain = f'{gen_dir}hg19ToHg38.over.chain.gz'

    bed_file = f'{bed_dir}{case_id}.bed'

    hgnc_dump = '20220817_hgnc_dump.txt'
    # hgnc_df = import_hgnc_dump(hgnc_dump)

    """ define filtering parameters """

    qual = '20'  # for initial QUAL value filtering

    parameters = {
        'max_af': 0.05,  # seen in 5% of the population at most
        'min_cadd': 10,  # CADD >= 10 implies p(var not observed) >= 0.9
        'max_pli': 0.9,  # higher pLI implies less tolerant to truncation
        'max_oe_lof': 0.35,  # lower o/e implies gene less tolerant to variants
        'max_oe_mis': 0.35,
        'max_oe_syn': 0.35}

    """ define some semi-random variants for testing """

    test_vars = [
        '11:2159793:T:C',  # pathogenic INS var (synonymous)
        '19:17843086:G:T',  # pathogenic JAK3 var (missense)
        '13:32398437:C:G',  # pathogenic BRCA2 var (LOF, stop gained)
        '11:34452212:G:A',  # pathogenic CAT var (splice region)
        '17:43051062:C:T',  # pathogenic BRCA1 var (splice donor)
        '17:43049196:T:G',  # pathogenic BRCA1 var (splice acceptor)
        '17:43104967:A:C',  # pathogenic BRCA1 var (intronic)
        '14:21525060:G:A',  # VOUS in SALL2
        '7:117559479:G:A']  # most common CFTR variant in gnomad]

    """ function calls """

    vcf = f'{vcf_prefix}_1_proband.vcf.gz'
    b37_vcf = f'{vcf_dir}SAP-48034-1_1_proband.vcf'

    json_file = 'do_not_upload/data_structure/opencga_rd37_CSA-1000-1.json'

    look_at_json(json_file)
    # look_at_100k_cases(case_list)
    # look_at_vcf(vcf)
    # look_at_vcf(b37_vcf)


if __name__ == '__main__':
    main()
