#!usr/bin/env python

"""
Program to process a single EGLH 100k rare disease case.
Called within caerus.sh.

Inputs
    case_id [str]
    json [fp]: 100k JSON containing case data
    vcfs [fp array]: SNV vcf.gz files available for this case

Intermediate files


Outputs
    variants [fp]: excel workbook of variant data
"""


import gzip
import json
import os
import pandas as pd
import re
import shutil
import subprocess
import sys

from datetime import datetime as dt
from functools import wraps
from panelapp import api, Panelapp, queries
from pycellbase.cbconfig import ConfigClient
from pycellbase.cbclient import CellBaseClient


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


""" preliminary checks """


def check_same_case(case_id, json_data, vcfs):
    """ confirm case ids match between all inputs """

    json_id = json_data['case_id']

    assert case_id == json_id, \
        f"JSON case ID {json_id} doesn't match {case_id}."

    for vcf in vcfs:

        vcf_id = do something here

        assert vcf_id == case_id, \
            f"{vcf} case ID {vcf_id} doesn't match {case_id}."


def check_rd_case(case_id, json_data):
    """ confirm rare disease case """

    assert json_data['program'].strip().lower() == 'rare_disease', \
        f"{case_id} is not a rare disease case."


""" process json """


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
            'samples': {},
            'vcfs': {}},
        'initial_analysis': {
            'panels': [],
            'variants': {
                'reported': [],
                'b38_lifted': [],
                'annotated': {}},
            'solved': {
                'fully': False,
                'partially': False}},
        'reanalysis': {
            'panels': [],
            'parameters': {},
            'variants': {
                'filtered': [],
                'annotated': {}}}}

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

    # identify whether the proband is adopted

    adopted = True

    for person in json_data['interpretation_request_data']['json_request'][
        'pedigree']['members']:

        if person['isProband']:
            if person['adoptedStatus'].lower() == 'notadopted':
                adopted = False
            elif person['adoptedStatus'].lower() == 'adopted':
                adopted = True

    # get info for (a) proband (b) all family if proband is not adopted

    for person in json_data['interpretation_request_data']['json_request'][
        'pedigree']['members']:

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
                case['family']['samples'][relation] = samples
                case['family']['vcfs'][relation] = {}
                vcf_paths = []

                for sample in case['family']['samples'][relation]:

                    vcfs = [vcf['uriFile'] for vcf in json_data[
                        'interpretation_request_data']['json_request']['vcfs']
                        if (sample in vcf['sampleId'])
                        and (vcf['fileType'] == 'VCF_small')]  # only SNVs

                    for vcf in vcfs:
                        vcf_paths.append(vcf)

                if vcf_paths:
                    assert len(vcf_paths) == 1, \
                        f"{case['id']} {relation} has multiple VCF paths"

                    case['family']['vcfs'][relation]['path'] = vcf_paths[0]

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

    for report in json_data['clinicalReports']:
        if report['valid']:
            if report['exitQuestionnaire']:

                solved = report['exitQuestionnaire'][
                    'familyLevelQuestions']['caseSolvedFamily']

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
def create_bed(case, bed_fp):
    """ Given a set of PanelApp panel IDs, retrieve data for the current
    versions of those panels and parse out their genes and regions. Use
    these to create a bed file via the Ensembl BioMart (GRCh38) API.
    Save this at the location specified by bed_fp.

    args:
        case [dict]
        bed_fp [str]: fp to write bed file to
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

    with open('bed_template.txt', 'r') as reader:
        contents = reader.read()

    hgnc_string = ','.join([ele for ele in genes])
    query = contents.replace('PLACEHOLDER', hgnc_string)

    # query biomart using subprocess

    subprocess.run(['wget', '-O', bed_fp, query])

    # add panel regions to the new bed file

    with open(bed_fp, 'r') as reader:
        bed_data = pd.read_csv(reader, sep='\t', header=None)

    bed_data.columns = ['chrom', 'start', 'end']
    region_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])
    updated_bed = pd.concat([bed_data, region_df])

    # sort gene/region dataframe by chromosome and position

    sorted = updated_bed.sort_values(by=['chrom', 'start'])

    with open(bed_fp, 'w') as writer:
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
    converted to standard notation. In addition, change the INFO field
    for any potential causal variants as a means of tracking them.

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

    if ref_genome == 'GRCh37':
        print(f'{vcf} requires liftover to GRCh38')
    elif ref_genome == 'GRCh38':
        print(f'{vcf} does not require liftover to GRCh38')

    return ref_genome


@time_execution
def mark_pcvs(input_vcf, output_vcf, pcvs):
    """ For a VCF file requiring liftover to GRCh38, tag potential
    causal variants identified in the original analysis so that their
    coordinates can be identified following liftover.

    args:
        input_vcf [fp]: vcf to be lifted over
        output_vcf [fp]: path to save output to
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
    pcv_line = '##INFO=<ID=PCV,Number=1,Type=Integer,' \
        'Description="Potential causal variant">\n'

    for idx, line in enumerate(lines):

        if idx % 1000000 == 0:
            print(f"{idx} lines processed")

        if line.startswith('#'):
            if line.startswith('#CHROM'):
                new_lines.append(pcv_line)
            new_lines.append(line)

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

    pcvs_present.sort()

    with gzip.open(output_vcf, 'wt') as writer:
        for line in new_lines:
            writer.write(line)

    return pcvs_present


@time_execution
def lift_over_vcf(vcf_in, vcf_out, genome_file, chain_file, pcvs):
    """ Liftover a b37 VCF file to b38 using CrossMap.

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
    """ Use the CB client to retrieve variant information.

    args:
        vcf_fp [fp]: variants to annotate

    returns:
        anno [list]: each element is a dict of annotation for 1 variant
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
def condense_annotation(anno, pcvs):
    """ Extract specific data values from an annotated list of variants
    if present.

    args:
        anno [list]: full annotation of all SNVs from sorted VCF
        pcvs [list of str]: variants identified in original analysis

    returns:
        anno_less [list]: all variants, but only necessary data kept
        anno_pcvs [list]: annotation for originally found variants
    """

    anno_less = []
    anno_pcvs = []

    for variant in anno:

        results = variant['results'][0]

        var_dict = {
            'variant': variant['id'].strip(),
            'consequence': None,
            'allele_frequency': None,
            'scaled_cadd': None,
            'affected_transcripts': {},
            'constraint': {},
            'vda': []}

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
            for assn in variant['geneTraitAssociation']:

                var_dict['vda'].append({
                    'id': assn['id'],
                    'name': assn['name'],
                    'score': assn['score'],
                    'sources': assn['sources']})

        except KeyError:
            pass

        anno_less.append(var_dict)

        if var_dict['variant'] in pcvs:
            anno_pcvs.append(var_dict)

    return anno_less, anno_pcvs


# filter
# de novo analysis
# segregation analysis

""" return output """

@time_execution
def create_excel(case, parameters, output_fp):
    """  """

    # create an excel file with 3 sheets

    # sheet 1: case data
    #   solved fully/partially?
    #   family members with vcfs: relation to proband, affected, adopted

    # sheet 2: variant data
    #   variant info: id, type, AF, CADD score, transcripts, constraint
    #   which family members it occurred in
    #   highlighted in original analysis or reanalysis

    # sheet 3: analysis data
    #   panel(s) used, list original and current versions
    #   all genes & regions analysed, and whether in original/reanalysis
    #   filtering parameters

    return something


""" coordinate & execute """


def main():

    # Show Liv a run-through of the pipeline
    # Get the VCFs and upload them too
    # Set up an app which runs the program
    # Case information
    #   How many solved (fully or partially) and unsolved cases
    #   Average number of highlighted variants
    #       all cases, solved cases, unsolved cases (presumably zero)
    #   Average number of snv vcfs per case
    #   proportion with vcfs from both parents
    #       among all cases, solved cases, unsolved cases
    # Fix de novo filtering
    # Fix segregation filtering
    # Check all cases can be processed through the pipeline
    # Optimise filtering parameters - check previously highlighted variants are always identified, without returning lots of others…

    # define static files

    chrom_map = 'chrom_map.txt'
    gen_dir = 'outside_github/genomes/'
    genome_file = f'{gen_dir}GCF_000001405.40_GRCh38.p14_genomic.fna'
    chain_file = f'{gen_dir}hg19ToHg38.over.chain.gz'

    # define filtering parameters

    parameters = {
        'qual': '20',  # for initial QUAL value filtering
        'max_af': 0.05,  # seen in <=5% of the population (gnomAD)
        'min_cadd': 10,  # CADD >= 10 implies p(var not observed) >= 0.9
        'min_pli': 0.9,  # higher implies less tolerant to truncation (gnomAD)
        'max_oe_lof': 0.35,  # lower implies less tolerant to variants (gnomAD)
        'max_oe_mis': 0.35,
        'max_oe_syn': 0.35,
        'min_vda': 0.75}  # score=0.75 ~= 2+ curated sources (DisGeNET)

    # read in arguments to process a single case

    case_id = sys.argv[1]
    json_fp = sys.argv[2]
    vcf_str = sys.argv[3]
    vcfs = vcf_str.split(' ')

    with open(json_fp, 'r') as reader:
        json_data = json.load(reader)

    # define output filenames





if __name__ == "__main__":
    main()
