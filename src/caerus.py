#!usr/bin/env python

"""
Program to process a single EGLH 100k rare disease case, filter and
prioritise variants from the original VCF output, and return potential
causal variants.

Inputs
    case_json [fp]: 100k JSON containing case data
    case_vcfs [fp array]: SNV vcf.gz files available for this case

Intermediate files


Outputs
    <case id>_report.xlsx [fp]: excel workbook of variant data
"""


import gzip
import json
import os
import pandas as pd
import re
import subprocess
import sys

from datetime import datetime as dt
from functools import wraps
from panelapp import api, Panelapp
from pycellbase.cbclient import CellBaseClient
from pycellbase.cbconfig import ConfigClient


""" utility functions """


def time_execution(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        t_start = dt.now()
        print(f'{t_start} Starting function: {func.__name__}')
        output = func(*args, **kwargs)
        t_end = dt.now()
        t_diff = t_end - t_start
        print(f'Duration: {t_diff}\n')
        return output
    return wrapper


def generate_date():
    """ Returns current date as string in the form YYYYMMDD """

    current_date = dt.today()
    date = str(dt.strftime(current_date, "%Y%m%d"))

    return date


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


""" process json """


def download_files(case_id):
    """ Download the JSON payload and any SNV VCFs associated with the
    specified case from the 003_230124_caerus DNAnexus project.

    args:
        case_id [str]

    returns:
        json_fp [str]
    """

    subprocess.run(["src/file_downloader.sh", case_id])

    for file in os.listdir():
        if file.endswith('.json') and (case_id in file):
            json_fp = file

    assert json_fp, f"Error: JSON file for {case_id} not downloaded"

    return json_fp


def check_rd_case(json_data):
    """ confirm rare disease case """

    assert json_data['program'].strip().lower() == 'rare_disease', \
        f"{json_data['case_id']} is not a rare disease case."


def setup_case_dict(json_data):
    """ Initialise the dictionary which will hold all of the case
    information.

    args:
        json_data [dict]

    returns:
        case [dict]
    """

    case = {
        'id': json_data['case_id'],
        'family': {
            'affected': {},
            'adopted': {},
            'vcfs': {
                'original': {}}},
        'initial_analysis': {
            'panels': [],
            'variants': {
                'reported': [],
                'b38': []},
            'solved': {
                'fully': False,
                'partially': False}},
        'reanalysis': {
            'panels': [],
            'variants': {}}}

    return case


def get_family_info(json_data, case):
    """ Get information on whether each family member is affected or
    adopted, and the IDs of relevant DNA samples and VCFs.

    args:
        json_data [dict]
        case [dict]

    returns:
        case [dict]
    """

    family = json_data['interpretation_request_data'][
        'json_request']['pedigree']['members']

    # identify whether the proband is adopted

    adopted = True

    for person in family:
        if person['isProband']:
            if person['adoptedStatus'].lower() == 'notadopted':
                adopted = False

    # get info for proband, or all family if proband is not adopted

    for person in family:
        if person['isProband'] or not adopted:

            # get relationship to proband

            if person['isProband']:
                relation = 'proband'
            else:
                relation = person['additionalInformation'][
                    'relation_to_proband'].lower()

            if relation == 'unrelated':
                continue  # don't use data for unrelated individuals

            # get adopted status

            case['family']['adopted'][relation] = True

            if person['adoptedStatus'] == 'notadopted':
                case['family']['adopted'][relation] = False

            if ('sibling' in relation) and case['family']['adopted'][relation]:
                continue  # omit data for any adopted siblings

            # get affection status

            affect_map = {
                'AFFECTED': 'True',
                'UNCERTAIN': 'Unknown',
                'UNAFFECTED': 'False'}

            case['family']['affected'][relation] = affect_map[
                person['affectionStatus']]

            # get IDs of all DNA samples, use to get paths to vcfs

            samples = [s['sampleId'] for s in person['samples']
                if s['product'] == 'DNA']

            if samples:
                vcfs = []

                for vcf in json_data['interpretation_request_data'][
                    'json_request']['vcfs']:

                    for sample in samples:

                        if (sample in vcf['sampleId']) and \
                            (vcf['fileType'] == 'VCF_small'):

                            vcfs.append(vcf['uriFile'])

                if vcfs:
                    assert len(vcfs) == 1, \
                        f"{case['id']} {relation} has multiple VCF paths"

                    pattern = r'(?<=Variations\/)(.*?)(?=$)'
                    vcf_ids = re.findall(pattern, vcfs[0])

                    case['family']['vcfs']['original'][relation] = vcf_ids[0]

    return case


def get_original_panels(json_data, case):
    """ Update a case dict with information on which panels were used
    for the original analysis, and their versions.

    args:
        json_data [dict]
        case [dict]

    returns:
        case [dict]
    """

    for panel in json_data['interpretation_request_data']['json_request'][
        'pedigree']['analysisPanels']:

        case['initial_analysis']['panels'].append({
            'name': panel['specificDisease'],
            'id': panel['panelName'],
            'version': panel['panelVersion']})

    assert len(case['initial_analysis']['panels']) >= 1, "Case has no panels"

    return case


def get_current_panels(case):
    """ Retrieve data for the current version of each of the panels used
    in the original analysis.

    args:
        case [dict]

    returns:
        case [dict]
    """

    for panel in case['initial_analysis']['panels']:

        panel_data = get_panelapp_panel(panel['id'])

        case['reanalysis']['panels'].append({
            'name': panel_data['name'],
            'id': panel_data['id'],
            'version': panel_data['version'],
            'data': panel_data})

    return case


def get_case_solved_info(json_data, case):
    """ Get information on whether the case is solved, and any known or
    potential causal variants.

    args:
        json_data [data]
        case [dict]

    returns:
        case [dict]
    """

    variants = []
    var_types = [
        'variants',
        'structuralVariants',
        'shortTandemRepeats',
        'chromosomalRearrangements']

    for report in json_data['clinical_report']:
        if report['valid']:
            if report['exit_questionnaire']:

                solved = report['exit_questionnaire'][
                    'exit_questionnaire_data']['familyLevelQuestions'][
                    'caseSolvedFamily']

                if solved == 'yes':
                    case['initial_analysis']['solved']['fully'] = True
                    case['initial_analysis']['solved']['partially'] = True

                elif solved == 'unknown':
                    case['initial_analysis']['solved']['partially'] = True

            # identify potentially causal variants

            for var_type in var_types:
                try:
                    for var in report['clinicalReportData'][var_type]:

                        chrom = var['variantCoordinates']['chromosome']
                        pos = var['variantCoordinates']['position']
                        ref = var['variantCoordinates']['reference']
                        alt = var['variantCoordinates']['alternate']

                        coords = f"{chrom}:{pos}:{ref}:{alt}"

                        if coords not in variants:
                            variants.append(coords)

                except Exception:
                    pass

    case['initial_analysis']['variants']['reported'] = variants

    return case


""" create bed file """


@time_execution
def create_bed(case, bed_template, bed_output):
    """ Given a set of PanelApp panel IDs, retrieve data for the current
    versions of those panels and parse out their genes and regions. Use
    these to create a bed file via the Ensembl BioMart (GRCh38) API.
    Save this at the location specified by bed_output.

    args:
        case [dict]
        bed_template [fp]: location of template for api call
        bed_output [fp]: location to write output bed file to
    """

    # list all unique panel genes and regions

    genes = []
    regions = []

    for panel in case['reanalysis']['panels']:

        panel_genes = get_panel_genes(panel['data'])
        panel_regions = get_panel_regions(panel['data'])

        for gene in panel_genes:
            if gene not in genes:
                genes.append(gene)

        for region in panel_regions:
            if region not in regions:
                regions.append(region)

    # read in template biomart query and insert string of HGNC ids

    with open(bed_template, 'r') as reader:
        contents = reader.read()

    hgnc_string = ','.join([ele for ele in genes])
    query = contents.replace('PLACEHOLDER', hgnc_string)

    # query biomart using subprocess

    subprocess.run(['wget', '-O', bed_output, query])

    # add panel regions to the new bed file

    with open(bed_output, 'r') as reader:
        bed_data = pd.read_csv(reader, sep='\t', header=None)

    bed_data.columns = ['chrom', 'start', 'end']
    region_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    updated_bed = pd.concat([bed_data, region_df])

    # sort gene/region dataframe by chromosome and position

    sorted = updated_bed.sort_values(by=['chrom', 'start'])

    with open(bed_output, 'w') as writer:
        sorted.to_csv(writer, sep='\t', header=True, index=False)


def get_panel_genes(panel_data):
    """ Given data for a single panel retrieved from PanelApp, parse out
    and return a list of all unique green gene target HGNC ids.

    args:
        panel_data [dict]: output of query to PanelApp

    returns:
        genes [list]
    """

    genes = []

    for gene in panel_data['genes']:

        if gene['confidence_level'] == '3' and \
            gene['gene_data']['hgnc_id'] and \
            (gene['gene_data']['hgnc_id'] not in genes):

            genes.append(gene['gene_data']['hgnc_id'])

    return genes


def get_panel_regions(panel_data):
    """ Given data for a single panel retrieved from PanelApp, parse out
    and return a list of all unique green regions.

    args:
        panel_data [dict]: output of query to PanelApp

    returns:
        regions [list]: each ele has form [<chrom>, <start>, <end>]
    """

    regions = []

    for region in panel_data['regions']:

        if region['confidence_level'] == '3':

            locus = [
                region['chromosome'],
                region['grch38_coordinates'][0],
                region['grch38_coordinates'][1]]

            regions.append(locus)

    return regions


""" standardise vcfs """


@time_execution
def fix_chroms(input_vcf, output_vcf, chrom_map):
    """ Modify a VCF file so that any UCSC-style CHROM values are
    converted to standard notation.

    args:
        input_vcf [fp]
        output_vcf [fp]
        chrom_map [fp]: map of different versions of chromosome notations
    """

    subprocess.run([
        'bcftools', 'annotate',
        '--rename-chrs', chrom_map,
        '-o', output_vcf,
        input_vcf])


@time_execution
def identify_genome(vcf):
    """ Identify the reference genome used to produce a VCF file.

    args:
        vcf [fp]: the file to examine

    returns:
        ref_genome [str]: GRCh37 or GRCh38
    """

    ref_genome = None

    if vcf.endswith('.gz'):
        with gzip.open(vcf, "rt") as reader:
            contents = reader.read()
    else:
        with open(vcf, 'r') as reader:
            contents = reader.read()

    all_refs = re.findall("GRCh37|GRCh38", contents)

    assert len(set(all_refs)) == 1, f"{vcf} mentions multiple ref genomes."

    ref_genome = all_refs[0]

    return ref_genome


@time_execution
def mark_pcvs(input_vcf, output_vcf, pcvs):
    """ For a VCF file requiring liftover to GRCh38, tag potential
    causal variants identified in the original analysis so that their
    coordinates can be identified following liftover.

    args:
        input_vcf [fp]: vcf to be lifted over
        output_vcf [fp]: path to save output marked vcf to
        pcvs [list]: list of variants in form '<chrom>:<pos>:<ref>:<alt>'

    returns:
        pcvs_present [list]: pcvs which are present in this vcf
    """

    if input_vcf.endswith('.gz'):
        with gzip.open(input_vcf, "rt") as reader:
            lines = reader.readlines()
    else:
        with open(input_vcf, 'r') as reader:
            lines = reader.readlines()

    new_lines = []
    pcvs_present = []
    pcv_header_line = '##INFO=<ID=PCV,Number=1,Type=Integer,' \
        'Description="Potential causal variant">\n'

    for idx, line in enumerate(lines):

        if idx % 1000000 == 0:
            print(f"{idx} lines processed")

        # include all header lines in the output vcf
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                new_lines.append(pcv_header_line)
            new_lines.append(line)

        # modify body lines if for pcvs, otherwise include unmodified
        else:
            split = line.split('\t')
            chrom = split[0].strip()
            pos = split[1].strip()
            ref = split[3].strip()
            alt = split[4].strip()
            var = f'{chrom}:{pos}:{ref}:{alt}'

            if var in pcvs:
                split[7] == 'PCV=1'
                pcvs_present.append(var)
                new_line = '\t'.join(split)
                new_lines.append(new_line)
            else:
                new_lines.append(line)

    pcvs_present.sort()

    with gzip.open(output_vcf, 'wt') as writer:
        for line in new_lines:
            writer.write(line)

    return pcvs_present


@time_execution
def lift_over_vcf(vcf_in, vcf_out, genome_file, chain_file, pcvs):
    """ Liftover a b37 VCF file to b38 using CrossMap. If potential
    causal variants were identified in the original analysis, confirm
    that they are still present after liftover.

    args:
        vcf_in [fp]
        vcf_out [fp]
        genome_file [fp]: reference genome fasta, required for liftover
        chain_file [fp]: required for liftover
        pcvs [list]: variants as chrom:pos:ref:alt (b37 notation)

    returns:
        pcvs_b38 [list]: variants as chrom:pos:ref:alt (b38 notation)
    """

    pcvs_b38 = []

    if vcf_in.endswith('.gz'):
        file_type = 'gvcf'
    else:
        file_type = 'vcf'

    subprocess.run([
        'CrossMap.py',
        file_type,
        chain_file,
        vcf_in,
        genome_file,
        vcf_out])

    if pcvs:
        pcvs_b38 = check_liftover_pcvs(vcf_out)

    assert pcvs == pcvs_b38, f"PCVs lost during liftover of {vcf_in}."

    return pcvs_b38


def check_liftover_pcvs(vcf):
    """ Following VCF liftover, identify all variants with an INFO field
    value of 'PCV=1' (indicating that they were identified as potential
    causal variants in the original analysis).

    args:
        vcf [fp]

    returns:
        pcvs_present [list]
    """

    pcvs_present = []

    if vcf.endswith('.gz'):
        with gzip.open(vcf, "rt") as reader:
            lines = reader.readlines()
    elif vcf.endswith('.vcf'):
        with open(vcf, 'r') as reader:
            lines = reader.readlines()

    for line in lines:
        if line[0] != '#':
            split = line.split('\t')
            if split[7] == 'PCV=1':

                chrom = split[0].strip()
                pos = split[1].strip()
                ref = split[3].strip()
                alt = split[4].strip()

                pcvs_present.append(f'{chrom}:{pos}:{ref}:{alt}')

    pcvs_present.sort()

    return pcvs_present


@time_execution
def filter_on_bed(input_vcf, bed, qual, out_prefix):
    """  Use vcftools to filter variants in a VCF file on QUAL value and
    a bed file.

    args:
        vcf_file [str]: input file to be filtered
        bed_file [str]: listing regions to filter on
        qual [str]: minimum QUAL value to retain a variant
        output_prefix [str]: prefix for output files

    returns:
        filtered_vcf [fp]
    """

    if input_vcf.endswith('.gz'):
        file_type = '--gzvcf'
    else:
        file_type = '--vcf'

    subprocess.run([
        'vcftools',
        file_type,
        input_vcf,
        '--bed', bed,
        '--minQ', qual,
        '--recode',
        '--recode-INFO-all',
        '--out', out_prefix,
        ])

    return f'{out_prefix}.recode.vcf'


@time_execution
def sort_vcf(vcf_fp, output_fp):
    """ Use bcftools to sort a VCF file.

    args:
        vcf_fp [str]: input VCF file
        output_fp [str]: path to sorted output file
        pcvs_b38 [list]: in format chrom:pos:ref:alt
    """

    subprocess.run([
        'bcftools',
        'sort',
        '-o', output_fp,
        vcf_fp,
        ])


""" annotate & filter VCFs """


@time_execution
def cb_client_annotation(vcf_fp):
    """ Use the CellBase client to retrieve variant information.

    args:
        vcf_fp [fp]: variants to annotate

    returns:
        anno [list]: list of dicts, each holds annotation for 1 variant
    """

    var_list = []

    with open(vcf_fp, 'r') as reader:
        lines = reader.readlines()

    for line in lines:
        if line[0] != '#':  # ignore header lines
            line_data = [value.strip() for value in line.split('\t')]
            chrom = line_data[0]
            pos = line_data[1]
            ref = line_data[3]
            alt = line_data[4]
            string = f'{chrom}:{pos}:{ref}:{alt}'
            var_list.append(string)

    cc = ConfigClient("config.json")
    cbc = CellBaseClient(cc)
    var_client = cbc.get_variant_client()

    anno = var_client.get_annotation(var_list, include=[
        # 'conservation',
        # 'cytoband',
        # 'hgvs',
        # 'variation',  # equivalent to 'id'
        # 'mirnaTargets',  # 'geneMirnaTargets'
        # 'drugInteraction',  # 'geneDrugInteraction'
        'functionalScore',
        'geneConstraints',
        'populationFrequencies',  # 'populationFrequencies' + 'id'
        'consequenceType',  # 'consequenceTypes' + 'displayConsequenceType'
        'geneDisease'])  # 'geneTraitAssociation'

    return anno


@time_execution
def condense_annotation(anno_in, pcvs):
    """ Extract specific data values from an annotated list of variants
    if present.

    args:
        anno_in [list]: full annotation of all SNVs from sorted VCF
        pcvs [list of str]: variants identified in original analysis

    returns:
        anno_out [list]: all variants, but only necessary data kept
        anno_pcvs [list]: annotation for originally found variants
    """

    anno_out = []
    anno_pcvs = []

    for variant in anno_in:

        results = variant['results'][0]

        var_dict = {
            'variant': variant['id'].strip(),
            'consequence': None,
            'allele_frequency': None,
            'scaled_cadd': None,
            'affected_transcripts': {}}

        # get the general variant consequence (e.g. missense_variant)
        try:
            var_dict['consequence'] = results['displayConsequenceType']

        except KeyError:
            pass

        # get allele frequency (global gnomAD AF)
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

        anno_out.append(var_dict)

        if var_dict['variant'] in pcvs:
            anno_pcvs.append(var_dict)

    return anno_out, anno_pcvs


@time_execution
def apply_filters(variants, parameters):
    """ Given annotated variant data and a set of filtering parameters,
    return only the variants which meet filtering thresholds.

    args:
        variants [list]: 1 dict of annotation per variant
        parameters [dict]: values to use for each filter

    returns:
        output_variants [dict]: subset which passed filtering
    """

    output_variants = {}
    high_risk = ['frameshift_variant', 'stop_gained', 'stop_lost']

    for variant in variants:

        # variant must affect a transcript
        if variant['affected_transcripts']:

            # variant AF must be below threshold, or not observed
            if (variant['allele_frequency'] <= parameters['max_af']) or \
                not variant['allele_frequency']:

                # rare transcript-affecting frameshifts always pass
                if variant['consequence'] in high_risk:
                    output_variants[variant['id']] = variant

                # CADD score must be above threshold
                elif variant['scaled_cadd'] and \
                    (variant['scaled_cadd'] >= parameters['min_cadd']):

                        output_variants[variant['id']] = variant

    return output_variants


# de novo analysis
# segregation analysis

""" return output """

@time_execution
def create_excel(case, parameters, output_fp):
    """  """


""" coordinate & execute """


def get_cases_list(args):
    """ Parse the arguments from the call to caerus.py to obtain a list
    of cases to process.

    args:
        args [list]: has the value of sys.argv[1:]

    returns:
        cases [list]: ids of cases to process
    """

    assert args, \
        "Call to caerus.py requires at least one argument. Examples:\n" \
        "\tA single case ID\n" \
        "\tA space-separated list of case IDs\n" \
        "\tFilepath of a .txt file containing a case ID on each line\n\n" \
        "Alternatively, use the argument 'all' to process all cases for " \
        "which files are available in 003_230124_caerus.\n\n" \
        "Case IDs should include their version, e.g. SAP-55997-1."

    if len(args) == 1:

        # process all cases for which we have a json and at least one vcf
        if args[0] == 'all':

            subprocess.run(['list_project_cases.sh'])

            with open('usable_project_cases.txt', 'r') as reader:
                cases = [case.strip() for case in reader.readlines()
                    if case.strip()]

        # process a list of cases supplied as lines of a text file
        elif args[0].endswith('.txt'):

            with open(args[0], 'r') as reader:
                cases = [case.strip() for case in reader.readlines()
                    if case.strip()]

        # otherwise assume the arg is a single case id
        else:
            cases = [args[0]]

    # if there are multiple args, assume they're a list of case ids
    else:
        cases = args[0:]

    return cases


def initialise_case(case_id):
    """ Coordinate the functions to initialise the object which holds
    all the relevant case information.

    args:
        case_id [str]

    returns:
        case [dict]
    """

    json_fp = download_files(case_id)

    with open(json_fp, 'r') as reader:
        json_data = json.load(reader)

    check_rd_case(json_data)

    case = setup_case_dict(json_data)
    case = get_family_info(json_data, case)
    case = get_original_panels(json_data, case)
    case = get_current_panels(case)
    case = get_case_solved_info(json_data, case)

    case_vcfs = case['family']['vcfs']['original']
    print(f"Case: {case_id}\nJSON: {json_fp}\nVCFs: {case_vcfs}")

    return case


def process_case(case, chrom_map, genome_file, chain_file, parameters):

    """ Process the case's VCFs. Involves standardisation of chromosome
    notation and reference gneomes; annotation; and filtering and
    prioritising of variants.

    Creates an output excel workbook with details of the analysis and
    the final prioritised variants.

    args:
        case [dict]
        chrom_map [str]: file
        genome_file
        chain_file
        parameters
    """

    case_id = case['id']
    case_vcfs = case['family']['vcfs']['original']
    initial_pcvs = case['initial_analysis']['variants']['reported']

    # process each family member's vcf

    for relation, vcf in case_vcfs.items():

        assert vcf in os.listdir(), \
            f"{vcf} ({relation}) not present in directory"

        fixed_vcf = f"{case_id}_2a_fixed_{relation}.vcf.gz"
        marked_vcf = f"{case_id}_2b_marked_{relation}.vcf.gz"
        lifted_vcf = f"{case_id}_3_b38_{relation}.vcf.gz"
        vcf_prefix = f"{case_id}_4_filtered_{relation}.vcf.gz"
        sorted_vcf = f"{case_id}_5_sorted_{relation}.vcf.gz"
        json_vars = f'{case_id}_6_{relation}_ann.json'

        # fix chromosome notation, identify ref genome used

        fix_chroms(vcf, fixed_vcf, chrom_map)
        ref_genome = identify_genome(fixed_vcf)

        # b37 vcfs require liftover

        if ref_genome == "GRCh37":

            initial_vcf_pcvs = mark_pcvs(
                fixed_vcf, marked_vcf, initial_pcvs)

            lifted_vcf_pcvs = lift_over_vcf(
                marked_vcf, lifted_vcf, genome_file, chain_file, initial_pcvs)

            case['family']['vcfs']['b38'][relation] = lifted_vcf

            if relation == 'proband':

                assert initial_vcf_pcvs == lifted_vcf_pcvs, \
                    f"PCVs lost during liftover.\n" \
                    f"Before liftover: {initial_vcf_pcvs}\n" \
                    f"After liftover: {lifted_vcf_pcvs}"

                case['initial_analysis'][
                    'variants']['b38'] = lifted_vcf_pcvs

        # b38 are fine as they are

        elif ref_genome == "GRCh38":
            case['family']['vcfs']['b38'][relation] = fixed_vcf

            if relation == 'proband':

                assert initial_vcf_pcvs == initial_pcvs, \
                    f"PCVs missing from fixed proband VCF.\n" \
                    f"Reported: {initial_pcvs}\n" \
                    f"After fixing: {initial_vcf_pcvs}"

                case['initial_analysis'][
                    'variants']['b38'] = initial_vcf_pcvs

        b38_vcf = case['family']['vcfs']['b38'][relation]
        b38_pcvs = case['initial_analysis']['variants']['b38']

        # filter b38 vcfs on bed file and qual score, then sort

        bed_filtered_vcf = filter_on_bed(
            b38_vcf, f"{case_id}.bed", parameters['qual'], vcf_prefix)

        sort_vcf(bed_filtered_vcf, sorted_vcf)

        # annotate and filter variants in proband vcf

        if relation == 'proband':

            try:  # saving variant annotations reduces API calls during testing
                with open(json_vars, 'r') as reader:
                    ann_original = json.load(reader)

            except FileNotFoundError:
                ann_original = cb_client_annotation(sorted_vcf)
                with open(json_vars, 'w') as writer:
                    json.dump(ann_original, writer)

            ann_less = condense_annotation(ann_original, b38_pcvs)
            ann_final = apply_filters(ann_less, parameters)
            case['reanalysis']['variants'] = ann_final

    # # create output file
    # output_fp = f"{dt.now()}_reanalysis_{case_id}.xlsx"
    # create_excel(case, parameters, output_fp)

    output_fp = f"{dt.now()}_reanalysis_{case_id}.json"
    with open(output_fp, 'w') as writer:
        json.dump(case, writer)


def main():
    # define static files

    res_dir = "resources/home/dnanexus/"
    bed_template = f"{res_dir}bed_template.txt"
    chrom_map = f"{res_dir}chrom_map.txt"
    genome_file = f"{res_dir}GCF_000001405.40_GRCh38.p14_genomic.fna"
    chain_file = f"{res_dir}hg19ToHg38.over.chain.gz"

    # define filtering parameters

    parameters = {
        'qual': '20',  # for initial QUAL value filtering
        'max_af': 0.05,  # seen in <=5% of the population (gnomAD)
        'min_cadd': 10}  # CADD >= 10 implies p(var not observed) >= 0.9

    # process cases

    cases = get_cases_list(sys.argv[1:])
    print(f"Cases to be processed: {cases}")

    for case_id in cases:

        case = initialise_case(case_id)

        bed_output = f"{case_id}.bed"
        create_bed(case, bed_template, bed_output)

        process_case(case_id, chrom_map, genome_file, chain_file, parameters)

    ## MAKE THE EXCEL OUTPUT JAY
    ## get tier of original variants


if __name__ == "__main__":
    main()
