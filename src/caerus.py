#!usr/bin/env python

"""
Program to process a single EGLH 100k rare disease case, filter and
prioritise variants from the original VCF output, and return potential
causal variants.

INPUTS AND USAGE

caerus.py takes several possible inputs via the command line:

    caerus.py <case id>
    caerus.py <case id 1> <case id 2> <case id 3> ...
    caerus.py <file>.txt
    caerus.py all

The first options two analyse one or multiple specified 100k cases, the
third reads in a file containing one case id per line, and the last
processes every case in DNAnexus for which a JSON file and at least one
SNV VCF are available.

INTERMEDIATE AND OUTPUT FILES

    Case JSON file - downloaded from DNAnexus 003_230124_caerus project
    Case SNV VCFs - as with JSON, number will vary with case
    BED file - generated based on panels listed in JSON via BioMart

Each VCF analysed will generate several intermediates:

    Standardised VCF - chromosome notation fixed, reported variants marked
    Sorted VCF (1)
    (Liftover VCF) - if original VCF was against GRCh37
    Filtered VCF - BED file applied and variants filtered on call quality
    Sorted VCF (2)

Example of duration and output

    caerus.py SAP-32023-1

Processing this case using only the proband VCF, which does not require
liftover, takes ~7.5 minutes and generates files totalling ~2.5GB. THIS
DID NOT INCLUDE THE ANNOTATION STEP OR FILTERING ON ANNOTATION.

"""


import os
import re
import sys
import json
import subprocess
import pandas as pd

from functools import wraps
from datetime import datetime as dt
from panelapp import api, Panelapp
from pycellbase.cbclient import CellBaseClient
from pycellbase.cbconfig import ConfigClient


""" utility functions """


def time_execution(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        t_start = dt.now()
        print(f'\n{t_start} Starting function: {func.__name__}')
        output = func(*args, **kwargs)
        t_end = dt.now()
        t_diff = t_end - t_start
        print(f'{t_end} {func.__name__} completed, duration: {t_diff}\n')
        return output
    return wrapper


def generate_date():
    """ Returns current date as string in the form YYYYMMDD """

    current_date = dt.today()
    date = str(dt.strftime(current_date, "%Y%m%d"))

    return date


def count_variant_lines(vcf_fp):
    """ Given the path to a VCF file, return the number of variants
    (i.e. non-header lines).

    args:
        vcf_fp [str]

    returns:
        var_count [int]
    """

    result = subprocess.run(
        f"bcftools view -H {vcf_fp} | wc -l",
        shell=True,
        capture_output=True,
        text=True).stdout

    count = int(result)
    print(f"{vcf_fp} contains {count} variants")

    return count


def get_variant_lines(vcf_fp):
    """ Given the path to a VCF file, return a list consisting of the
    lines describing variants (i.e. the VCF body).

    args:
        vcf_fp [str]

    returns:
        var_lines [list]
    """

    vcf_body = subprocess.run(
        f"bcftools view -H {vcf_fp}",
        shell=True,
        capture_output=True,
        text=True).stdout

    var_lines = [line.strip() for line in vcf_body.split('\n') if line.strip()]

    return var_lines


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


""" create a bed file """


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

    for panel_id, panel_dict in case['panels'].items():

        panel_genes = get_panel_genes(panel_dict['retrieved_data'])
        panel_regions = get_panel_regions(panel_dict['retrieved_data'])

        for gene in panel_genes:
            if gene not in genes:
                genes.append(gene)

        for region in panel_regions:
            if region not in regions:
                regions.append(region)

    # biomart can't deal with lots of genes at once, so split into groups

    if len(genes) > 500:
        groups = [genes[x:x+500] for x in range(0, len(genes), 500)]

    else:
        groups = [genes]

    print(f"Creating bed file of {len(genes)} genes in {len(groups)} groups")

    # get bed output for each group of genes

    bed_dfs = []

    for gene_group in groups:

        with open(bed_template, 'r') as reader:
            contents = reader.read()

        hgnc_string = ','.join([ele for ele in gene_group])
        query = contents.replace('PLACEHOLDER', hgnc_string)

        subprocess.run(['wget', '-O', 'resources/temp_bed_file.bed', query])

        with open('resources/temp_bed_file.bed', 'r') as reader:
            bed_data = pd.read_csv(reader, sep='\t', header=None)

        bed_data.columns = ['chrom', 'start', 'end']
        bed_dfs.append(bed_data)

    # concatenate gene and region dataframes into a single bed file

    region_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    bed_dfs.append(region_df)

    updated_bed = bed_dfs[0]

    for df in bed_dfs[1:]:
        updated_bed = pd.concat([updated_bed, df])

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


""" process a case json """


@time_execution
def download_files(case_id):
    """ Download the JSON payload and any SNV VCFs associated with the
    specified case from the 003_230124_caerus DNAnexus project.

    args:
        case_id [str]

    returns:
        json_fp [str]
    """

    subprocess.run(f"src/file_downloader.sh {case_id}", shell=True)

    for file in os.listdir():
        if file.endswith('.json') and (case_id in file):
            json_fp = file

    assert json_fp, f"Error: JSON file for {case_id} not downloaded"

    return json_fp


def check_rd_case(json_data):
    """ Confirm that the case is part of the rare disease program (rather
    than cancer). """

    assert json_data['program'].strip().lower() == 'rare_disease', \
        f"{json_data['case_id']} is not a rare disease case."


def check_proband_vcf(json_data):  # CHECK PROBAND VCF EXISTS BEFORE DOING ANYTHING ELSE
    """  """


def setup_case_dict(json_data, params, static):
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
                'original': {},
                'b38': {},
                'filtered': {}}},
        'solved': {
            'fully': False,
            'partially': False},
        'reanalysis': {
            'params': params,
            'static': static,
            'variants': {}},
        'panels':{},
        'variants': {
            'original': {},
            'b38': [],
            'final': {}}}

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

            # get the proband's disorder

            # disorder for disorder in:
            # data['interpretation_request_data']['json_request']['pedigree']['diseasePenetrances'][x]['specificDisease']

            # find variants and their classifications via:
            # data['clinical_report'][all]['exit_questionnaire']['exit_questionnaire_data']['variantGroupLevelQuestions'][all]['variantLevelQuestions'][all]

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

        case['panels'][panel['panelName']] = {
            'original_name': panel['specificDisease'],
            'original_version': panel['panelVersion'],
            'retrieved_id': '',
            'retrieved_name': '',
            'retrieved_version': ''}

    assert len(case['panels'].keys()) >= 1, "Case has no panels"

    return case


def get_current_panels(case):
    """ Retrieve data for the current version of each of the panels used
    in the original analysis.

    args:
        case [dict]

    returns:
        case [dict]
    """

    for panel_id, panel_dict in case['panels'].items():

        panel_data = get_panelapp_panel(panel_id)

        if panel_data:
            panel_dict['retrieved_id'] = panel_data['id']
            panel_dict['retrieved_name'] = panel_data['name']
            panel_dict['retrieved_version'] = panel_data['version']
            panel_dict['retrieved_data'] = panel_data

        else:
            panel_dict['retrieved_id'] = None
            panel_dict['retrieved_name'] = None
            panel_dict['retrieved_version'] = None
            panel_dict['retrieved_data'] = None

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

    original_vars = {}

    for report in json_data.get('clinical_report'):
        if report and report['valid']:

            exit_q = report.get('exit_questionnaire').get(
                'exit_questionnaire_data')

            # identify whether case is fully or partially solved

            solved = exit_q.get('familyLevelQuestions').get('caseSolvedFamily')

            if solved == 'yes':
                case['solved']['fully'] = True
                case['solved']['partially'] = True

            elif solved == 'unknown':
                case['solved']['partially'] = True

            # identify variants reported in the original analysis

            if exit_q.get('variantGroupLevelQuestions'):
                for group in exit_q.get('variantGroupLevelQuestions'):

                    for lst in ['variantLevelQuestions',
                        'structuralVariantLevelQuestions',
                        'shortTandemRepeatLevelQuestions']:

                        for var in group.get(lst):

                            intpn = var.get('acmgClassification')

                            chrom = var['variantCoordinates']['chromosome']
                            pos = var['variantCoordinates']['position']
                            ref = var['variantCoordinates']['reference']
                            alt = var['variantCoordinates']['alternate']

                            coords = f"{chrom}:{pos}:{ref}:{alt}"

                            if coords not in original_vars.keys():
                                original_vars[coords] = {'acmg': intpn}

    case['variants']['original'] = original_vars

    return case


def check_proband_vcf(case):
    """ Confirm that the proband's VCF has been downloaded. """

    proband_vcf = case['family']['vcfs']['original']['proband']

    assert proband_vcf in os.listdir(), \
        f"Proband VCF {proband_vcf} not present for case {case['id']}"


""" standardise a vcf """


@time_execution
def sort_vcf(input_vcf, output_vcf):
    """ Use bcftools to sort a VCF file.

    args:
        vcf_fp [str]: input VCF file
        output_fp [str]: path to sorted output file
        pcvs_b38 [list]: in format chrom:pos:ref:alt
    """

    if output_vcf not in os.listdir():
        subprocess.run(
            f"bcftools sort -o {output_vcf} {input_vcf}", shell=True)

    else:
        print(f'Skipping sorting as {output_vcf} already exists')


@time_execution
def bgzip_vcf(vcf):
    """  """

    subprocess.run(f"gunzip {vcf}", shell=True)
    subprocess.run(f"bgzip {vcf[:-3]}", shell=True)


@time_execution
def normalise_vcf(input_vcf, output_vcf):
    """  """

    if output_vcf not in os.listdir():
        subprocess.run(
            f"bcftools norm -m - -o {output_vcf} {input_vcf}", shell=True)

    else:
        print(f'Skipping normalisation as {output_vcf} already exists')


@time_execution
def rename_vcf_chroms(input_vcf, output_vcf, chrom_map):
    """  """

    if output_vcf not in os.listdir():
        subprocess.run(
            f"bcftools annotate --rename-chrs {chrom_map} " \
            f"-o {output_vcf} {input_vcf}",
            shell=True)

    else:
        print(f'Skipping renaming as {output_vcf} already exists')


@time_execution
def tag_original_vars(case, input_vcf, output_vcf):
    """ Given a list of variants identified in the original analysis,
    create a tab-delimited text file of the original variants, and then
    annotate those variants within the VCF.

    args:
        case [dict]
        input_vcf
        output_vcf
    """

    if output_vcf not in os.listdir():

        # create bgzipped annotation file

        original_vars = case['variants']['original'].keys()
        anno_lines = []
        anno_fp = f"{case['id']}_mark_pcvs.tab"

        for var in original_vars:
            split = [e.strip() for e in var.split(':')]
            line = f"{split[0]}\t{split[1]}"
            anno_lines.append(line)

        with open(anno_fp, 'w') as writer:
            writer.write('#CHROM\tPOS\n')
            writer.write('\n'.join(anno_lines))

        subprocess.run(f'bgzip -f {anno_fp}', shell=True)

        # create tabix index for annotation file

        subprocess.run(f'tabix -s1 -b2 -e2 {anno_fp}.gz', shell=True)

        # annotate originally reported variants

        subprocess.run(
            f"bcftools annotate " \
            f"-a {anno_fp}.gz " \
            f"-c CHROM,POS " \
            f"-h resources/header.txt " \
            f"-m PCV=1 " \
            f"-o {output_vcf} {input_vcf}",
            shell=True)

    else:
        print(f'Skipping tagging as {output_vcf} already exists')


@time_execution
def get_tagged_vars(vcf):
    """  """

    # get lines where the tag is present

    tagged = subprocess.run(
        f"bcftools view -H -i'PCV=1' {vcf}",
        shell=True, capture_output=True, text=True).stdout

    tagged_lines = [val.strip() for val in tagged.split('\n') if val.strip()]

    # get the variant id from the vcf line

    tagged_vars = []

    for line in tagged_lines:
        split = [e.strip() for e in line.split('\t')]
        chrom = split[0]
        pos = split[1]
        ref = split[3]
        alt = split[4]

        tagged_vars.append(f"{chrom}:{pos}:{ref}:{alt}")

    tagged_vars.sort()

    return tagged_vars


@time_execution
def identify_genome(vcf):
    """ Identify the reference genome used to produce a VCF file.

    args:
        vcf [fp]: the file to examine

    returns:
        ref_genome [str]: GRCh37 or GRCh38
    """

    ref_genome = None

    # look for mentions of either ref genome in the vcf

    b37_lines = subprocess.run(
        'zgrep -cio "GRCh37" LP3000402-DNA_C06_LP3000402-DNA_C06.vcf.gz',
        shell=True, capture_output=True, text=True).stdout

    b38_lines = subprocess.run(
        'zgrep -cio "GRCh38" LP3000402-DNA_C06_LP3000402-DNA_C06.vcf.gz',
        shell=True, capture_output=True, text=True).stdout

    # there should be at least one mention of the build

    assert (int(b37_lines) > 0) or (int(b38_lines) > 0), \
        f"Error identifying genome for {vcf}: no ref genome specified"

    # identify which ref build is used

    if (int(b37_lines) > 0) and (int(b38_lines) == 0):
        ref_genome = "GRCh37"

    elif (int(b38_lines) > 0) and (int(b37_lines) == 0):
        ref_genome = "GRCh38"

    # both ref builds shouldn't be used for the same vcf

    assert ref_genome, \
        f"Error identifying genome for {vcf}: both ref genomes mentioned"

    print(f"{vcf} uses {ref_genome}")

    return ref_genome


@time_execution
def lift_over_vcf(input_vcf, output_vcf, chain, genome):
    """ Liftover a b37 VCF file to b38 using CrossMap. If potential
    causal variants were identified in the original analysis, confirm
    that they are still present after liftover.

    args:
        case [dict]
        input_vcf [fp]
        output_vcf [fp]

    returns:
        b38_vars [list]: variants as chrom:pos:ref:alt (b38 notation)
    """

    if output_vcf not in os.listdir():

        if input_vcf.endswith('.gz'):
            file_type = 'gvcf'
        else:
            file_type = 'vcf'

        subprocess.run(
            f"CrossMap.py {file_type} {chain} {input_vcf} {genome} {output_vcf}",
            shell=True)

    else:
        print(f'Skipping liftover as {output_vcf} already exists')


@time_execution
def filter_on_bed(input_vcf, bed, params, out_prefix):
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

    if f'{out_prefix}.recode.vcf' not in os.listdir():

        if input_vcf.endswith('.gz'):
            file_type = '--gzvcf'
        else:
            file_type = '--vcf'

        subprocess.run(
            f"vcftools {file_type} {input_vcf} --bed {bed} " \
            f"--minQ {params['qual']} --minGQ {params['gq']} " \
            f"--minDP {params['depth']} --recode --recode-INFO-all " \
            f"--out {out_prefix}",
            shell=True)

    else:
        print(f'Skipping filtering as {out_prefix}.recode.vcf already exists')

    return f'{out_prefix}.recode.vcf'


@time_execution
def filter_on_affection(case):
    """ Uses bedtools intersect to identify variants present in the VCFs
    of all affected family members.

    args:
        case [dict]

    returns:
        case [dict]: paths of VCFs filtered on segregation added
    """

    initial_vcf = case['family']['vcfs']['filtered']['proband']
    output_filename = f"{case['id']}_intersect.vcf.gz"

    # don't use family filtering if proband is adopted

    if not case['family']['adopted']['proband']:
        for member, member_vcf in case['family']['vcfs']['filtered'].items():

            # don't bother intersecting the proband's vcf against itself

            if member == 'proband':
                pass

            # perform intersect on each available vcf in turn

            elif (case['family']['affected'][member] == 'True') and \
                not case['family']['adopted'][member]:

                output_filename = (output_filename[:-7] + f"_{member}.vcf.gz")

                subprocess.run(
                    f"bedtools intersect -a {initial_vcf} -b {member_vcf} " \
                    f"| gzip > {output_filename}",
                    shell=True)

                initial_vcf = output_filename

    case['family']['vcfs']['intersect'] = output_filename

    return case


def process_single_individual(case, relation):
    """ Standardise the VCF for a single individual including fixing
    chromosome notation and lifting over to GRCh38, before filtering on
    variant call quality metrics and a bed file.

    args:
        case [dict]
        relation [str]

    returns:
        case [dict]
    """

    print(f"\nProcessing {relation} VCF for case {case['id']}")

    case_id = case['id']
    params = case['reanalysis']['params']
    static = case['reanalysis']['static']
    vcf = case['family']['vcfs']['original'][relation]

    original_vars = [key for key in case['variants']['original']]
    original_vars.sort()

    sorted_vcf_1 = f"{case_id}_1_sorted_{relation}.vcf.gz"
    normalised_vcf = f"{case_id}_2_normalised_{relation}.vcf.gz"
    renamed_vcf = f"{case_id}_3_renamed_{relation}.vcf.gz"
    tagged_vcf = f"{case_id}_4_tagged_{relation}.vcf.gz"
    lifted_vcf = f"{case_id}_5_lifted_{relation}.vcf.gz"
    prefix = f"{case_id}_6_filtered_{relation}"
    sorted_vcf_2 = f"{case_id}_7_sorted_{relation}.vcf.gz"

    # standardise vcf

    sort_vcf(vcf, sorted_vcf_1)
    bgzip_vcf(sorted_vcf_1)
    normalise_vcf(sorted_vcf_1, normalised_vcf)
    rename_vcf_chroms(normalised_vcf, renamed_vcf, static['chrom_map'])

    # if the were originally reported variants, tag them in the vcf

    if original_vars and (relation == 'proband'):
        tag_original_vars(case, renamed_vcf, tagged_vcf)
        tagged_vars = get_tagged_vars(tagged_vcf)
        fixed_vcf = tagged_vcf

        assert tagged_vars == original_vars, \
            f"Original: {original_vars}\nTagged: {tagged_vars}"

    else:
        fixed_vcf = renamed_vcf

    # identify the reference genome used, lift over if it was b37

    ref_genome = identify_genome(renamed_vcf)

    if ref_genome == "GRCh37":
        lift_over_vcf(fixed_vcf, lifted_vcf, static['chain'], static['genome'])
        case['family']['vcfs']['b38'][relation] = lifted_vcf

        # check there are as many b38 vars as were originally reported
        # (can't compare them directly because they have different coords)

        if original_vars and (relation == 'proband'):
            b38_vars = get_tagged_vars(lifted_vcf)
            case['variants']['b38'] = b38_vars

            assert len(b38_vars) == len(original_vars), \
                f"{len(original_vars)} original tagged vars " \
                f"but {len(b38_vars)} tagged b38 vars"

        elif not original_vars:
            case['variants']['b38'] = original_vars

    # b38 are fine as they are

    elif ref_genome == "GRCh38":
        case['family']['vcfs']['b38'][relation] = fixed_vcf
        case['variants']['b38'] = original_vars

    b38_vcf = case['family']['vcfs']['b38'][relation]

    # filter on bed file and variant call quality, then sort again

    bed_filtered_vcf = filter_on_bed(b38_vcf, f"{case_id}.bed", params, prefix)
    sort_vcf(bed_filtered_vcf, sorted_vcf_2)

    case['family']['vcfs']['filtered'][relation] = sorted_vcf_2

    return case


""" annotate & filter VCFs """


@time_execution
def cb_client_annotation(vcf_fp, cb_config):
    """ Use the CellBase client to retrieve variant information.

    args:
        vcf_fp [fp]: variants to annotate

    returns:
        anno [list]: list of dicts, each holds annotation for 1 variant
    """

    var_list = []
    lines = get_variant_lines(vcf_fp)

    for line in lines:

        c_p = re.findall(r'^(.*?)\t(.*?)\t', line)  # chrom & pos
        r_a = re.findall(r'\t([ACGNT]*)\t([ACGNT]*)\t', line)  # ref & alt

        variant = f"{c_p[0][0]}:{c_p[0][1]}:{r_a[0][0]}:{r_a[0][1]}"
        var_list.append(variant)

    cc = ConfigClient(cb_config)
    cbc = CellBaseClient(cc)
    var_client = cbc.get_variant_client()

    anno = var_client.get_annotation(var_list, include=[
        # 'conservation',
        # 'cytoband',
        # 'hgvs',
        # 'variation',  # equivalent to 'id'
        # 'mirnaTargets',  # 'geneMirnaTargets'
        # 'drugInteraction',  # 'geneDrugInteraction'
        # 'geneConstraints',
        # 'geneDisease'  # 'geneTraitAssociation'
        'functionalScore',
        'populationFrequencies',  # 'populationFrequencies' + 'id'
        'consequenceType'])  # 'consequenceTypes' + 'displayConsequenceType'

    return anno


@time_execution
def condense_annotation(input_vars):
    """ Extract specific data values from an annotated list of variants
    if present.

    args:
        input_vars [list]: full annotation of all SNVs from sorted VCF

    returns:
        output_vars [list]: dicts for same variants, but reduced info
    """

    output_vars = []

    for variant in input_vars:

        results = variant['results'][0]

        var_dict = {
            'id': variant['id'].strip(),
            'consequence': None,
            'gnomad_af': None,
            'scaled_cadd': None,
            'transcripts': {},
            'original': False,
            'reanalysis': False}

        # get the general variant consequence (e.g. missense_variant)
        try:
            var_dict['consequence'] = results['displayConsequenceType']

        except KeyError:
            pass

        # get allele frequency (global gnomAD AF)
        try:
            for dct in results['populationFrequencies']:
                if dct['population'] == 'ALL':
                    var_dict['gnomad_af'] = dct['altAlleleFreq']

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

            var_dict['transcripts'] = genes

        except KeyError:
            pass

        output_vars.append(var_dict)

    return output_vars


@time_execution
def apply_filters(input_vars, params, original_vars):
    """ Given annotated variant data and a set of filtering parameters,
    return only the variants which meet filtering thresholds.

    args:
        input_vars [list]: 1 dict of annotation per variant
        params [dict]: values to use for each filter
        original_vars [list]

    returns:
        output_vars [list]: subset which passed filtering
    """

    output_vars = []
    high_risk = ['frameshift_variant', 'stop_gained']

    for var in input_vars:

        if var['id'] in original_vars:
            var['original'] = True

        # affects a transcript
        if var['transcripts']:

            # low population af or not observed
            if (var['gnomad_af'] <= params['max_af']) or \
                not var['gnomad_af']:

                # frameshift, or high cadd score
                if var['consequence'] in high_risk or \
                    (var['scaled_cadd'] and \
                    (var['scaled_cadd'] >= params['min_cadd'])):

                    var['reanalysis'] = True

        if var['original'] or var['reanalysis']:
            output_vars.append(var)

    return output_vars


""" return output """

@time_execution
def create_excel(case, output_fp):
    """ Generate an Excel workbook to display the results of a reanalysis.

    args:
        case [dict]: containing reanalysis information
        output_fp [str]: to write xlsx file to
    """

    df_dict = create_dfs(case)

    writer = pd.ExcelWriter(output_fp, engine='xlsxwriter')

    # write the info, parameters and family dfs to the first sheet

    row = 0
    col = 0

    df_dict['info'].to_excel(
        writer,
        sheet_name='case_info',
        startrow=row,
        startcol=col,
        header=False)

    for df_row in df_dict['info'].iterrows():
            row += 1
    row += 1

    df_dict['param_df'].to_excel(
        writer,
        sheet_name='case_info',
        startrow=row,
        startcol=col,
        header=False)

    for df_row in df_dict['info'].iterrows():
            row += 1
    row += 1

    df_dict['family_df'].to_excel(
        writer,
        sheet_name='case_info',
        startrow=row,
        startcol=col,
        index=False)

    # write the variant df to the second sheet

    df_dict['var_df'].to_excel(
        writer,
        sheet_name='variants',
        index=False)

    # write the panels df to the third sheet

    df_dict['panels_df'].to_excel(
        writer,
        sheet_name='panels',
        index=False)


def create_dfs(case):
    """ Create dataframes to hold the data which will go in the output.

    args:
        case [dict]: holding all the data

    returns:
        output_dicts [dict]: dfs holding the info to go in the workbook
    """

    info = {
        'reanalysis_date': dt.today(),
        'case_id': case['id'],
        'solved_fully': case['solved']['fully'],
        'solved_partially': case['solved']['partially']}

    family = {'relation': [], 'affected': [], 'adopted': []}

    for relation, affected in case['family']['affected'].items():
        family['relation'].append(relation)
        family['affected'].append(affected)
        family['adopted'].append(case['family']['adopted'][relation])

    vars = {'id': [], 'consequence': [], 'gnomad_af': [], 'scaled_cadd': [],
        'transcripts': [], 'original': [], 'reanalysis': []}

    for var in case['variants']['final']:
        vars['id'].append(var['id'])
        vars['consequence'].append(var['consequence'])
        vars['gnomad_af'].append(var['gnomad_af'])
        vars['scaled_cadd'].append(var['scaled_cadd'])
        vars['transcripts'].append(var['transcripts'])
        vars['original'].append(var['original'])
        vars['reanalysis'].append(var['reanalysis'])

    panels = {'original_id': [], 'original_name': [], 'original_version': [],
        'retrieved_id': [], 'retrieved_name': [], 'retrieved_version': []}

    for panel_id, panel_dict in case['panels'].items:
        panels['original_id'].append(panel_id)
        panels['original_name'].append(panel_dict['original_name'])
        panels['original_version'].append(panel_dict['original_version'])
        panels['retrieved_id'].append(panel_dict['retrieved_id'])
        panels['retrieved_name'].append(panel_dict['retrieved_name'])
        panels['retrieved_version'].append(panel_dict['retrieved_version'])

    output_dicts = {
        'info_df': pd.Series(info),
        'param_df': pd.Series(case['reanalysis']['params']),
        'family_df': pd.DataFrame(family),
        'var_df': pd.DataFrame(vars),
        'panels_df': pd.DataFrame(panels)}

    return output_dicts


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
        "\nCall to caerus.py requires at least one argument. Examples:\n" \
        "\tA single case ID\n" \
        "\tA space-separated list of case IDs\n" \
        "\tFilepath of a .txt file containing a case ID on each line\n\n" \
        "Alternatively, use the argument 'all' to process all cases for " \
        "which files are available in 003_230124_caerus.\n\n" \
        "Case IDs should include their version, e.g. SAP-55997-1."

    if len(args) == 1:

        # process all cases for which we have a json and at least one vcf
        if args[0] == 'all':

            subprocess.run('list_project_cases.sh', shell=True)

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


def initialise_case(case_id, params, static):
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

    case = setup_case_dict(json_data, params, static)
    case = get_family_info(json_data, case)
    case = get_original_panels(json_data, case)
    case = get_current_panels(case)
    case = get_case_solved_info(json_data, case)

    check_proband_vcf(case)

    case_vcfs = case['family']['vcfs']['original']
    print(f"Case: {case_id}\nJSON: {json_fp}\nVCFs: {case_vcfs}")

    return case


def process_case(case, use_family):

    """ Process the case's VCFs. Involves standardisation of chromosome
    notation and reference genomes; annotation; and filtering and
    prioritising of variants.

    Creates an output excel workbook with details of the analysis and
    the final prioritised variants.

    args:
        case [dict]
        use_family [bool]
    """

    # process the proband's vcf

    case = process_single_individual(case, 'proband')
    vcf_to_annotate = case['family']['vcfs']['filtered']['proband']

    # optionally process family vcfs and intersect vcfs

    if use_family:
        for relation in case['family']['vcfs']['original'].keys():

            # only process vcfs for affected individuals who aren't adopted

            if (relation != 'proband'):
                if (case['family']['affected'][relation] == 'True'):
                    if not case['family']['adopted'][relation]:

                        case = process_single_individual(case, relation)

                    else:
                        print(f"{relation} VCF excluded (adopted)")
                else:
                    print(f"{relation} VCF excluded (not affected)")

        if len(case['family']['vcfs']['filtered'].keys()) > 1:

            case = filter_on_affection(case)
            vcf_to_annotate = case['family']['vcfs']['intersect']

    # annotate and filter variants in specified vcf
    # (saving variant annotations as jsons reduces API calls during testing)

    cb_config = case['reanalysis']['static']['cb_config']
    json_vars = f"{case['id']}_6_annotation.json"

    # try:
    #     with open(json_vars, 'r') as reader:
    #         ann_original = json.load(reader)

    # except FileNotFoundError:
    #     ann_original = cb_client_annotation(vcf_to_annotate, cb_config)
    #     with open(json_vars, 'w') as writer:
    #         json.dump(ann_original, writer)

    # ann_less = condense_annotation(ann_original)
    # ann_final = apply_filters(ann_less, case['reanalysis']['params'], b38_pcvs)
    # case['variants']['final'] = ann_final

    # # create output file
    # output_fp = f"{dt.now()}_reanalysis_{case_id}.xlsx"
    # create_excel(case, params, output_fp)

    output_fp = f"{case['id']}_temp_output.json"
    with open(output_fp, 'w') as writer:
        json.dump(case, writer)


def main():
    # define filtering parameters and static files

    static = {
        'bed': "resources/bed_template.txt",
        'chrom_map': "resources/chrom_map.txt",
        'genome': "resources/GCF_000001405.40_GRCh38.p14_genomic.fna",
        'chain': "resources/hg19ToHg38.over.chain.gz",
        'cb_config': 'resources/cb_config.json'}

    params = {
        'qual': '20',  # variant call quality
        'gq': '20',  # genotype quality
        'depth': '10',  # read depth at variant position
        'max_af': 0.01,  # seen in <=5% of the population (gnomAD)
        'min_cadd': 10}  # CADD >= 10 implies p(var not observed) >= 0.9

    # process cases

    cases = get_cases_list(sys.argv[1:])
    print(f"Cases to be processed: {cases}")

    for case_id in cases:

        case = initialise_case(case_id, params, static)
        bed = f"{case_id}.bed"

        if bed not in os.listdir():
            create_bed(case, static['bed'], bed)
        else:
            print(f'\nSkipping bed file creation as {bed} already exists')

        process_case(case, use_family=True)


if __name__ == "__main__":
    main()
