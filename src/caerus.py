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
from openpyxl import load_workbook
from openpyxl.styles import NamedStyle, Font, Alignment
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
    """ Returns current date as string in the form YYYY-MM-DD """

    current_date = dt.today()
    date = str(dt.strftime(current_date, "%Y-%m-%d"))

    return date


def get_vcf_section_lines(vcf, section):
    """ Given the path to a VCF file, return the header or variant lines
    as a list.

    args:
        vcf [str]: path to input vcf file
        section [str]: 'header' or 'variant'

    returns:
        lines [list]: each element is one line of the vcf header/body
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
        vcf [str]: path to input vcf file
        section [str]: 'header' or 'variant'

    returns:
        count [int]: number of lines in specific vcf section
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
def create_bed(case, bed_template, bed_output, output_dir):
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

    all_genes = []
    regions = []

    for panel_id, info in case['panels'].items():
        if info['retrieved_data']:

            panel_genes = get_panel_genes(info['retrieved_data'])
            panel_regions = get_panel_regions(info['retrieved_data'])

            info['current_genes'] = len(panel_genes)
            info['current_regions'] = len(panel_regions)

            all_genes += panel_genes

            if panel_regions:
                for region in panel_regions:
                    if region not in regions:
                        regions.append(region)

    if all_genes:
        genes = list(set(all_genes))

    assert genes or regions, \
        f"No current versions of original panels available for {case['id']}"

    print(f"Creating bed: {len(genes)} genes and {len(regions)} regions\n")

    bed_dfs = []

    # retrieve gene positions from ensembl biomart (GRCh38), make dataframes

    if genes:

        if len(genes) > 500:
            groups = [genes[x:x+500] for x in range(0, len(genes), 500)]
        else:
            groups = [genes]

        print(f"Genes processed in {len(groups)} group(s)")
        temp_bed_file = f'{output_dir}temp_bed_file.bed'

        for gene_group in groups:

            with open(bed_template, 'r') as reader:
                contents = reader.read()

            hgnc_string = ','.join([ele for ele in gene_group])
            query = contents.replace('PLACEHOLDER', hgnc_string)

            subprocess.run(['wget', '-O', temp_bed_file, query])

            with open(temp_bed_file, 'r') as reader:
                bed_data = pd.read_csv(reader, sep='\t', header=None)

            bed_data.columns = ['chrom', 'start', 'end']
            bed_dfs.append(bed_data)

    # turn regions into a dataframe

    if regions:
        region_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
        bed_dfs.append(region_df)

    # concatenate all dataframes into a single bed file

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

        if file.startswith('InterpretationDetail') and \
            file.endswith('.json') and \
            (case_id in file):

            json_fp = file

    assert json_fp, f"Error: JSON file for {case_id} not downloaded"

    return json_fp


def check_rd_case(json_data):
    """ Confirm that the case is part of the rare disease program (rather
    than cancer).

    args:
        json_data [dict]
    """

    assert json_data['program'].strip().lower() == 'rare_disease', \
        f"{json_data['case_id']} is not a rare disease case."


def setup_case_dict(json_data, params, static, date):
    """ Initialise the dictionary which will hold all of the case
    information.

    args:
        json_data [dict]
        params [dict]: filtering parameters to use
        static [dict]: paths to static files
        date [str]: date of reanalysis

    returns:
        case [dict]
    """

    case = {
        'id': json_data['case_id'],
        'assembly': json_data['assembly'].strip(),
        'disorders': [],
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
            'date': date,
            'params': params,
            'static': static},
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

    # identify whether the proband is affected and/or adopted

    adopted = True

    for person in family:
        if person['isProband']:

            assert person['affectionStatus'] == 'AFFECTED'

            if person['adoptedStatus'].lower() == 'notadopted':
                adopted = False

    # get the proband's disorder

    disorders = json_data['interpretation_request_data']['json_request']['pedigree']['diseasePenetrances']

    for dct in disorders:
        case['disorders'].append({
            'disorder': dct['specificDisease'],
            'penetrance': dct['penetrance']})

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
                        f"{case['id']} {relation} has multiple SNV VCFs"

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

    for panel_id, info in case['panels'].items():

        info['retrieved_id'] = None
        info['retrieved_name'] = None
        info['retrieved_version'] = None
        info['retrieved_data'] = None
        info['current_genes'] = None
        info['current_regions'] = None

        print(f"Retrieving panel {panel_id}: {info['original_name']}")

        panel_data = get_panelapp_panel(panel_id)

        if panel_data:
            info['retrieved_id'] = panel_data['id']
            info['retrieved_name'] = panel_data['name']
            info['retrieved_version'] = panel_data['version']
            info['retrieved_data'] = panel_data

    return case


def get_case_solved_info(json_data, case):
    """ Get information on whether the case is solved, and any
    originally reported variants along with their classifications.

    args:
        json_data [data]
        case [dict]

    returns:
        case [dict]
    """

    original_vars = {}

    var_types = ['variants', 'structuralVariants', 'shortTandemRepeats',
        'chromosomalRearrangements']

    for report in json_data.get('clinical_report'):
        if report and report['valid']:

            # identify whether case is fully or partially solved

            solved = None
            questions = None
            exit_q = report.get('exit_questionnaire')

            if exit_q:

                solved = exit_q.get('exit_questionnaire_data').get(
                    'familyLevelQuestions').get('caseSolvedFamily')

                questions = exit_q.get('variantGroupLevelQuestions')

            if solved == 'yes':
                case['solved']['fully'] = True
                case['solved']['partially'] = True

            elif solved == 'unknown':
                case['solved']['partially'] = True

            # get classifications of originally reported variants

            if questions:
                for group in exit_q.get('variantGroupLevelQuestions'):

                    for q_lst in ['variantLevelQuestions',
                        'structuralVariantLevelQuestions',
                        'shortTandemRepeatLevelQuestions']:

                        if group.get(q_lst):
                            for var in group.get(q_lst):

                                chrom = var['variantCoordinates']['chromosome']
                                pos = var['variantCoordinates']['position']
                                ref = var['variantCoordinates']['reference']
                                alt = var['variantCoordinates']['alternate']

                                var_id = f"{chrom}:{pos}:{ref}:{alt}"
                                intpn = var.get('acmgClassification')

                                if var_id not in original_vars.keys():
                                    original_vars[var_id] = {'acmg': intpn}

            # get info for originally reported variants

            data = report.get('clinical_report_data')

            if data:
                for var_type in var_types:
                    var_data = data.get(var_type)

                    if var_data:
                        for var in var_data:

                            # get variant id

                            chrom = var['variantCoordinates']['chromosome']
                            pos = var['variantCoordinates']['position']
                            ref = var['variantCoordinates']['reference']
                            alt = var['variantCoordinates']['alternate']

                            var_id = f"{chrom}:{pos}:{ref}:{alt}"

                            if var_id not in original_vars.keys():
                                original_vars[var_id] = {
                                    'tier': None,
                                    'explains_phen': None,
                                    'clin_sig': None,
                                    'clin_rev': None}

                            events = var.get('reportEvents')
                            attrs = var.get('variantAttributes')

                            # get variant tier / if it explains a phenotype

                            if events:
                                tier = list(set([e['tier'] \
                                    for e in events if e['tier']]))

                                exp = list(set([e['fullyExplainsPhenotype'] \
                                    for e in events \
                                    if e['fullyExplainsPhenotype']]))

                                original_vars[var_id]['tier'] = tier
                                original_vars[var_id]['explains_phen'] = exp

                            # get previous clinvar interpretations

                            if attrs:
                                anns = attrs.get(
                                    'additionalTextualVariantAnnotations')

                                if anns:
                                    rev = anns.get('clinvar_reviewStatus')
                                    sig = anns.get(
                                        'clinvar_clinicalSignificances')

                                    original_vars[var_id]['clin_sig'] = sig
                                    original_vars[var_id]['clin_rev'] = rev

                            if 'acmg' not in original_vars[var_id].keys():
                                original_vars[var_id]['acmg'] = None

    case['variants']['original'] = original_vars

    # solved cases should have at least 1 reported variant

    if solved == 'yes':
        assert case['variants']['original'], \
            f"{case['id']} is solved but has no originally reported variants"

    return case


def check_proband_vcf(case):
    """ Confirm that the proband's VCF has been downloaded.

    args:
        case [dict]
    """

    proband_vcf = case['family']['vcfs']['original']['proband']

    assert proband_vcf in os.listdir(), \
        f"Proband VCF {proband_vcf} not present for case {case['id']}"


""" standardise a vcf """


@time_execution
def sort_vcf(input_vcf, output_vcf):
    """ Use bcftools to sort a VCF file.

    args:
        input_vcf [str]
        output_vcf [str]
    """

    if not os.path.exists(output_vcf):
        subprocess.run(
            f"bcftools sort -o {output_vcf} {input_vcf}", shell=True)

    else:
        print(f'Skipping sorting ({output_vcf} already exists)')


@time_execution
def bgzip_vcf(vcf):
    """ Un-gzip a vcf, then bgzip it.

    args:
        vcf [str]
    """

    subprocess.run(f"gunzip {vcf}", shell=True)
    subprocess.run(f"bgzip {vcf[:-3]}", shell=True)


@time_execution
def normalise_vcf(input_vcf, output_vcf):
    """ Use bcftools to normalise a VCF and decompose multiallelic
    variants.

    args:
        input_vcf [str]
        output_vcf [str]
    """

    if not os.path.exists(output_vcf):
        subprocess.run(
            f"bcftools norm -m - -o {output_vcf} {input_vcf}", shell=True)

    else:
        print(f'Skipping normalisation ({output_vcf} already exists)')


@time_execution
def rename_vcf_chroms(input_vcf, output_vcf, chrom_map):
    """ Use bcftools to rename CHROM field values to standard notation.

    args:
        input_vcf [str]
        output_vcf [str]
        chrom_map [str]: path to file defining mapping between notations
    """

    if not os.path.exists(output_vcf):
        subprocess.run(
            f"bcftools annotate --rename-chrs {chrom_map} " \
            f"-o {output_vcf} {input_vcf}",
            shell=True)

    else:
        print(f'Skipping chromosome renaming ({output_vcf} already exists)')


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

    if not os.path.exists(output_vcf):

        # create bgzipped annotation file

        original_vars = case['variants']['original'].keys()
        anno_lines = []
        anno_fp = f"intermediate_files/{case['id']}_mark_pcvs.tab"

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
        print(f'Skipping variant tagging ({output_vcf} already exists)')


@time_execution
def get_tagged_vars(vcf):
    """ Identify variant lines in a VCF where the INFO field contains
    the tag PCV=1 (for putative causal variant).

    args:
        vcf [str]

    returns:
        tagged_vars [list]: ids of variants which have the tag
    """

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
        vcf [str]

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
        input_vcf [fp]
        output_vcf [fp]
        chain [str]: path to chain file for liftover
        genome [str]: path to fasta.gz for target genome (i.e. b38)
    """

    if not os.path.exists(output_vcf):

        if input_vcf.endswith('.gz'):
            file_type = 'gvcf'
        else:
            file_type = 'vcf'

        subprocess.run(
            f"CrossMap.py " \
            f"{file_type} {chain} {input_vcf} {genome} {output_vcf}",
            shell=True)

    else:
        print(f'Skipping liftover as {output_vcf} already exists')


@time_execution
def filter_on_bed(input_vcf, output_prefix, bed, params):
    """  Use vcftools to filter variants on variant call quality, and
    against the regions specified in a supplied bed file. Quality
    parameters filtered on are QUAL score (minQ), genotype quality (GQ),
    depth (DP), and whether FILTER = PASS (remove-filtered-all).

    args:
        input_vcf [str]
        output_prefix [str]: prefix for output files
        bed [str]: path to bed file of genomic regions
        params [dict]: contains filtering parameters

    returns:
        output_file [str]
    """

    output_file = f'{output_prefix}.recode.vcf'

    if not os.path.exists(output_file):

        if input_vcf.endswith('.gz'):
            file_type = '--gzvcf'
        else:
            file_type = '--vcf'

        subprocess.run(
            f"vcftools {file_type} {input_vcf} " \
            f"--bed {bed} " \
            f"--remove-filtered-all " \
            f"--minQ {params['min_qual']} " \
            f"--minGQ {params['min_gq']} " \
            f"--minDP {params['min_depth']} " \
            f"--recode --recode-INFO-all " \
            f"--out {output_prefix}",
            shell=True)

    else:
        print(
            f'Skipping filtering on bed file ({output_file} already exists)')

    return output_file


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


def process_single_individual(case, relation, output_dir):
    """ Standardise the VCF for a single individual including fixing
    chromosome notation and lifting over to GRCh38, before filtering on
    variant call quality metrics and a bed file.

    args:
        case [dict]
        relation [str]: relationship of individual to proband

    returns:
        case [dict]
    """

    print(f"\nProcessing {relation} VCF for case {case['id']}")

    case_id = case['id']
    ref_genome = case['assembly']
    params = case['reanalysis']['params']
    static = case['reanalysis']['static']
    vcf = case['family']['vcfs']['original'][relation]

    original_vars = [key for key in case['variants']['original']]
    original_vars.sort()

    bed = f"{output_dir}{case_id}.bed"
    sorted_vcf_1 = f"{output_dir}{case_id}_1_sorted_{relation}.vcf.gz"
    normalised_vcf = f"{output_dir}{case_id}_2_normalised_{relation}.vcf.gz"
    renamed_vcf = f"{output_dir}{case_id}_3_renamed_{relation}.vcf.gz"
    tagged_vcf = f"{output_dir}{case_id}_4_tagged_{relation}.vcf.gz"
    lifted_vcf = f"{output_dir}{case_id}_5_lifted_{relation}.vcf.gz"
    prefix = f"{output_dir}{case_id}_6_filtered_{relation}"
    sorted_vcf_2 = f"{case_id}_7_sorted_{relation}.vcf.gz"

    # standardise vcf

    sort_vcf(vcf, sorted_vcf_1)
    bgzip_vcf(sorted_vcf_1)
    normalise_vcf(sorted_vcf_1, normalised_vcf)
    rename_vcf_chroms(normalised_vcf, renamed_vcf, static['chrom_map'])

    # # identify the reference genome

    # ref_genome = identify_genome(renamed_vcf)

    # if ref genome is b37...

    if ref_genome.lower() == "grch37":

        # if the were originally reported variants, tag them in the vcf

        if original_vars and (relation == 'proband'):
            tag_original_vars(case, renamed_vcf, tagged_vcf)
            tagged_vars = get_tagged_vars(tagged_vcf)
            fixed_vcf = tagged_vcf

            assert tagged_vars == original_vars, \
                f"Original: {original_vars}\nTagged: {tagged_vars}"

        else:
            fixed_vcf = renamed_vcf

        # lift over to b38

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

    elif ref_genome.lower() == "grch38":
        case['family']['vcfs']['b38'][relation] = renamed_vcf
        case['variants']['b38'] = original_vars

    b38_vcf = case['family']['vcfs']['b38'][relation]

    # filter on bed file and variant call quality, then sort again

    bed_filtered_vcf = filter_on_bed(b38_vcf, prefix, bed, params)
    sort_vcf(bed_filtered_vcf, sorted_vcf_2)

    case['family']['vcfs']['filtered'][relation] = sorted_vcf_2

    return case


""" annotate & filter VCFs """


def split_large_vcfs(case_id, vcf, anno_dir):
    """ split large numbers of variants into smaller groups (can have
    memory issues otherwise)

    gets all variant lines from vcf
    splits into groups with set maximum size
    for each group, gets list of variant ids
    saves lists of variant ids to text files, 1 per group
    """

    # split variants into subgroups

    print('Splitting variants into subgroups...')

    lines = get_vcf_section_lines(vcf, 'variant')

    if len(lines) > 20000:
        groups = [lines[x:x+20000] for x in range(0, len(lines), 20000)]

    else:
        groups = [lines]

    group_count = len(groups)

    print(f"{len(lines)} variants to annotate in {group_count} groups")
    print("Splitting variants into groups...")

    # for each group, make a file listing the variant ids

    for idx, group in enumerate(groups):

        fp = f'{anno_dir}{case_id}_group_{idx}.txt'

        if not os.path.exists(fp):

            vars = []

            for line in group:

                split = [x.strip() for x in line.split('\t')]
                var_id = f"{split[0]}:{split[1]}:{split[3]}:{split[4]}"

                if var_id not in vars:
                    vars.append(var_id)

            with open(fp, 'w') as writer:
                writer.write('\n'.join(vars))

    return group_count


def get_annotation(idx, vars_fp, anno_fp, cb_config):
    """ get annotation for the variants listed in the supplied text file
    """

    assert os.path.exists(vars_fp), \
        f"Can't annotate {vars_fp}: file doesn't exist."

    # if the annotation file doesn't exist, create it

    if not os.path.exists(anno_fp):

        print(f"Retrieving annotation for group {idx}")

        with open(vars_fp, 'r') as reader:
            vars = [x.strip() for x in reader.readlines()]

        cc = ConfigClient(cb_config)
        cbc = CellBaseClient(cc)
        var_client = cbc.get_variant_client()

        # options for get_annotation 'include':

        # 'hgvs', 'cytoband', 'conservation', 'geneConstraints',
        # 'functionalScore','variation' (='id'),
        # 'mirnaTargets' (='geneMirnaTargets'),
        # 'geneDisease' (='geneTraitAssociation')
        # 'drugInteraction' (='geneDrugInteraction'),
        # 'populationFrequencies' (= 'populationFrequencies' + 'id')
        # 'consequenceType' (='consequenceTypes', 'displayConsequenceType')

        annotations = var_client.get_annotation(vars, include=[
            'functionalScore', 'consequenceType', 'populationFrequencies'])

        condensed = condense_annotation(annotations)

        with open(anno_fp, 'w') as writer:
            json.dump(condensed, writer)

    else:
        print(f"Skipping retrieval for group {idx}")


def condense_annotation(variants):
    """ Extract specific data values from an annotated list of variants
    if present.

    args:
        case_id [str]

    returns:
        output_vars [list]: dicts of selected variant info
    """

    output_vars = []

    for variant in variants:

        results = variant['results'][0]

        var_dict = {
            'id': variant['id'].strip(),
            'consequences': None,
            'gnomad_af': None,
            'scaled_cadd': None,
            'transcripts': None,
            'original': False,
            'reanalysis': False}

        # get allele frequency (global gnomAD AF)

        af = results.get('populationFrequencies')

        if af:
            for dct in af:
                if dct.get('population') == 'ALL':
                    var_dict['gnomad_af'] = dct.get('altAlleleFreq')

        # get scaled CADD scores

        functional = results.get('functionalScore')

        if functional :
            for dct in functional:
                if dct.get('source') == 'cadd_scaled':
                    var_dict['scaled_cadd'] = dct.get('score')

        # get affected genes and transcripts, and consequence types

        genes = {}
        consequences = []
        types = results.get('consequenceTypes')

        if types:
            for dct in types:

                gene = dct.get('geneName')
                transcript = dct.get('transcriptId')
                terms = dct.get('sequenceOntologyTerms')

                if gene:
                    if gene not in genes.keys():
                        genes[gene] = []

                    if transcript:
                        if transcript not in genes[gene]:
                            genes[gene].append(transcript)

                if terms:
                    for ele in terms:
                        term = (ele.get('name'), ele.get('accession'))

                        if term and (term not in consequences):
                            consequences.append(term)

        var_dict['transcripts'] = genes
        var_dict['consequences'] = consequences

        output_vars.append(var_dict)

    return output_vars


def combine_files(case_id, anno_dir, output_fp):
    """  """

    if not os.path.exists(output_fp):

        print('Combining annotations into single file...')

        variants = []

        for file in os.listdir(anno_dir):
            if file.startswith(case_id) and file.endswith('_annotation.json'):
                with open(f"{anno_dir}{file}", 'r') as reader:
                    data = json.load(reader)

                variants += data

        with open(output_fp, 'w') as writer:
            json.dump(variants, writer)

        return variants

    else:
        print(f"Skipping file combination ({output_fp} already exists)")

        with open(output_fp, 'r') as reader:
            variants = json.load(reader)

        return variants


def retrieve_cadd(variant):
    """ retrieve and return the cadd score for an SNV variant """

    cadd = None

    chrom, pos, ref, alt = variant['id'].split(':')
    var = f"{chrom}:{pos}_{ref}_{alt}"

    result = subprocess.run(
        'curl -i ' \
        f'https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.6/{var}',
        shell=True, capture_output=True, text=True).stdout

    if result:

        pattern = r'(?<="PHRED":")(.*?)(?=")'
        regex = re.findall(pattern, result)

        if regex and (len(regex) == 1):
            cadd = float(regex[0])

    return cadd


@time_execution
def annotation_processing(case_id, vcf, cb_config, output_dir):
    """  """

    # specify the location to store annotation output

    anno_dir = f"{output_dir}temp_anno_files/"
    output_fp = f"{anno_dir}{case_id}_combined_annotation.json"

    if not os.path.exists(anno_dir):
        subprocess.run(f"mkdir {anno_dir}", shell=True)

    # process variants in smaller subgroups for memory purposes

    if not os.path.exists(output_fp):

        var_groups = split_large_vcfs(case_id, vcf, anno_dir)

        for i in range(var_groups):

            vars_fp = f"{anno_dir}{case_id}_group_{i}.txt"
            anno_fp = f"{anno_dir}{case_id}_group_{i}_annotation.json"

            get_annotation(i, vars_fp, anno_fp, cb_config)

        ann_vars = combine_files(case_id, anno_dir, output_fp)

    else:
        print(f"Skipping annotation ({output_fp} already exists)")

    return ann_vars


@time_execution
def apply_filters(params, all_vars, original_vars):
    """ Given annotated variant data and a set of filtering parameters,
    return only the variants which meet filtering thresholds.

    args:
        input_vars [list]: annotated variants after bed/quality filtering
        params [dict]: values to use for each filter
        original_vars [list]: those reported in original analysis

    returns:
        output_vars [list]: variants after filtering on annotation
    """

    output_vars = []

    consequences = ['SO:0001893', 'SO:0001574', 'SO:0001575',
        'SO:0001589', 'SO:0001578', 'SO:0001582', 'SO:0001889', 'SO:0001821',
        'SO:0001822', 'SO:0001583', 'SO:0001630', 'SO:0001626', 'SO:0001587']

    for var in all_vars:

        # originally reported variants should always be returned

        if var['id'] in original_vars:
            var['original'] = True

        # identify whether variant is an snv

        snv = False
        split = var['id'].split(':')

        if (len(split[2]) == 1) and (len(split[3]) == 1):
            snv = True

        # filter 1: low or unrecorded allele frequency

        if not var['gnomad_af'] or (var['gnomad_af'] <= params['max_af']):

            # filter 2: high/medium consequence type

            for cons in var['consequences']:
                if (cons[1] in consequences):

                    if not snv and (var['gnomad_af'] or var['scaled_cadd']):
                        var['reanalysis'] = True

                    # filter 3 (SNVs only): high scaled cadd score

                    elif snv:
                        if not var['scaled_cadd']:
                            var['scaled_cadd'] = retrieve_cadd(var)

                        if var['scaled_cadd'] and \
                            (var['scaled_cadd'] >= params['min_cadd']):

                            var['reanalysis'] = True

        # add to output

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

    df_dict['info_df'].to_excel(
        writer,
        sheet_name='case_info',
        startrow=row,
        startcol=col,
        header=False)

    row += 7

    df_dict['param_df'].to_excel(
        writer,
        sheet_name='case_info',
        startrow=row,
        startcol=col,
        header=False)

    row += 6

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

    writer.save()
    format_workbook(output_fp)


def create_dfs(case):
    """ Create dataframes to hold the data which will go in the output.

    args:
        case [dict]: holding all the data

    returns:
        output_dicts [dict]: dfs holding the info to go in the workbook
    """

    original = case['variants']['original']

    # define case summary df

    summary = {
        'reanalysis_date': case['reanalysis']['date'],
        'case_id': case['id'],
        'disorders': case['disorders'],
        'required_liftover': False,
        'solved_fully': case['solved']['fully'],
        'solved_partially': case['solved']['partially']}

    if case['assembly'].lower() == 'grch37':
        summary['required_liftover'] = True

    # define family info df

    family = {'relation': [], 'affected': [], 'adopted': []}

    for relation in case['family']['vcfs']['original'].keys():
        family['relation'].append(relation)
        family['affected'].append(case['family']['affected'][relation])

        if case['family']['adopted'][relation]:
            family['adopted'].append('True')
        else:
            family['adopted'].append('False')

    # define variant df

    vars = {'id': [], 'consequences': [], 'gnomad_af': [], 'scaled_cadd': [],
        'transcripts': [], 'reanalysis': [], 'original': [], 'acmg': [],
        'tier': [], 'explains_phenotype': [], 'clinvar_rank': [],
        'clinvar_status': []}

    for var in case['variants']['final']:
        vars['id'].append(var['id'])
        vars['consequences'].append(var['consequences'])
        vars['gnomad_af'].append(var['gnomad_af'])
        vars['scaled_cadd'].append(var['scaled_cadd'])
        vars['transcripts'].append(var['transcripts'])

        for key in 'reanalysis', 'original':
            if var[key]:
                vars[key].append('True')
            else:
                vars[key].append('False')

        if var['id'] in original.keys():
            var_dct = original[var['id']]

            vars['acmg'].append(var_dct['acmg'])
            vars['tier'].append(var_dct['tier'])
            vars['explains_phenotype'].append(var_dct['explains_phen'])
            vars['clinvar_rank'].append(var_dct['clin_sig'])
            vars['clinvar_status'].append(var_dct['clin_rev'])

        else:
            for key in ['acmg', 'tier', 'explains_phenotype', \
                'clinvar_rank', 'clinvar_status']:

                vars[key].append('')

    # define panels df

    panels = {'original_id': [], 'original_name': [], 'original_version': [],
        'retrieved_id': [], 'retrieved_name': [], 'retrieved_version': [],
        'current_genes': [], 'current_regions': []}

    for panel_id, info in case['panels'].items():
        panels['original_id'].append(panel_id)
        panels['original_name'].append(info['original_name'])
        panels['original_version'].append(info['original_version'])
        panels['retrieved_id'].append(info['retrieved_id'])
        panels['retrieved_name'].append(info['retrieved_name'])
        panels['retrieved_version'].append(info['retrieved_version'])
        panels['current_genes'].append(info['current_genes'])
        panels['current_regions'].append(info['current_regions'])

    # collect into output

    output_dicts = {
        'info_df': pd.Series(summary),
        'param_df': pd.Series(case['reanalysis']['params']),
        'family_df': pd.DataFrame(family),
        'var_df': pd.DataFrame(vars),
        'panels_df': pd.DataFrame(panels)}

    return output_dicts


def format_workbook(output_fp):
    """ Apply formatting to the excel workbook output.

    args:
        output_fp [str]: path to workbook
    """

    wb = load_workbook(filename=output_fp)

    ws1 = wb['case_info']
    ws2 = wb['variants']
    ws3 = wb['panels']

    # define normal and header styles

    normal_font = NamedStyle(name="normal_font")
    normal_font.font = Font(name='Arial', size=10)

    header_font = NamedStyle(name="header_font")
    header_font.font = Font(name='Arial', size=10, bold=True)

    for font in normal_font, header_font:

        font.alignment = Alignment(horizontal='left', vertical='center')
        wb.add_named_style(font)

    # apply normal style across workbook

    for ws in wb.worksheets:
        for row in ws.rows:
            for cell in row:
                cell.style = normal_font

    # set column widths and apply header style to headings

    for row in range(1, 21):
        ws1[f'A{row}'].style = header_font

    for col in 'ABC':
        ws1[f'{col}14'].style = header_font
        ws1.column_dimensions[col].width = 20

    for col in 'ABCDEFGHIJKL':
        ws2[f'{col}1'].style = header_font
        ws2.column_dimensions[col].width = 20

    for col in 'ABCDEFGHIJK':
        ws3[f'{col}1'].style = header_font
        ws3.column_dimensions[col].width = 20

    wb.save(filename=output_fp)


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


def initialise_case(case_id, params, static, date):
    """ Coordinate the functions to initialise the object which holds
    all the relevant case information.

    args:
        case_id [str]
        params [dict]: parameter values to filter variants against
        static [dict]: paths to static files
        date [str]: date or reanalysis

    returns:
        case [dict]
    """

    json_fp = download_files(case_id)

    with open(json_fp, 'r') as reader:
        json_data = json.load(reader)

    check_rd_case(json_data)

    case = setup_case_dict(json_data, params, static, date)
    case = get_family_info(json_data, case)
    case = get_original_panels(json_data, case)
    case = get_current_panels(case)
    case = get_case_solved_info(json_data, case)

    check_proband_vcf(case)

    case_vcfs = case['family']['vcfs']['original']
    print(f"\nCase: {case_id}\nJSON: {json_fp}\nVCFs: {case_vcfs}")

    return case


def process_case(case, output_dir, use_family):

    """ Process the case's VCFs. Involves standardisation of chromosome
    notation and reference genomes; annotation; and filtering and
    prioritising of variants.

    Creates an output excel workbook with details of the analysis and
    the final prioritised variants.

    args:
        case [dict]
        use_family [bool]: whether to use family VCFs in filtering
    """

    # process the proband's vcf

    case = process_single_individual(case, 'proband', output_dir)

    vcf_to_annotate = case['family']['vcfs']['filtered']['proband']
    cb_config = case['reanalysis']['static']['cb_config']
    b38_vars = case['variants']['b38']

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

    # get variant annotations and use them to filter variants

    ann_vars = annotation_processing(
        case['id'], vcf_to_annotate, cb_config, output_dir)

    ann_final = apply_filters(case['reanalysis']['params'], ann_vars, b38_vars)

    case['variants']['final'] = ann_final

    # create output file

    output_fp = f"{case['id']}_reanalysis_{case['reanalysis']['date']}.xlsx"
    create_excel(case, output_fp)


def main():
    # define filtering parameters and static files

    static = {
        'bed': "resources/bed_template.txt",
        'chrom_map': "resources/chrom_map.txt",
        'genome': "resources/GCF_000001405.40_GRCh38.p14_genomic.fna",
        'chain': "resources/hg19ToHg38.over.chain.gz",
        'cb_config': 'resources/cb_config.json'}

    params = {
        'min_qual': '20',  # variant call quality
        'min_gq': '20',  # genotype quality
        'min_depth': '20',  # read depth at variant position
        'max_af': 0.001,  # seen in <=5% of the population (gnomAD)
        'min_cadd': 20}  # CADD >= 10 implies p(var not observed) >= 0.9

    # identify cases to process, set up output dict

    cases = get_cases_list(sys.argv[1:])
    print(f"Cases to be processed: {cases}")

    date = generate_date()

    cwd = os.getcwd()
    output_dir = f"{cwd}/intermediate_files/"

    if not os.path.exists(output_dir):
        subprocess.run('mkdir intermediate_files', shell=True)

    # process cases

    for case_id in cases:

        case = initialise_case(case_id, params, static, date)
        bed = f"{output_dir}{case_id}.bed"

        if not os.path.exists(bed):
            create_bed(case, static['bed'], bed, output_dir)
        else:
            print(f'\nSkipping bed file creation ({bed} already exists)')

        process_case(case, output_dir, use_family=True)


if __name__ == "__main__":
    main()
