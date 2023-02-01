#!usr/bin/env python

"""
Notes on current pipeline:
- Only looks at variants within coding regions

Run once per case:

    process json		<1s
    create bed		    <1s

Run once per VCF:

    identify genome		~15s
    fix chrom notation	~1m 30s
    (liftover)
    filter on bed		~10s
    sort vcf		    <1s
    annotate variants	~25s
    condense annotation <1s
    apply filters 		<1s
    create output       <1s

Run once per case if family data available:
    filter de novo
    filter segregation

Running the pipeline for a case with two b38 VCFs should take ~5m
"""


import gzip
import json
import os
import pandas as pd
import re
import subprocess

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


""" process a json file """


@time_execution
def process_case_json(json_fp):
    """ Given the path to a 100k case JSON file, extract the information
    needed for the reanalysis pipeline.

    args:
        json_fp [str]

    returns:
        data [dict]: parsed case information
    """

    case = None

    with open(json_fp, 'r') as reader:
        case_json = json.load(reader)

    prog = case_json['program']

    if prog == 'rare_disease':
        case = initialise_case_dict(case_json)
        case = get_case_panels(case_json, case)
        case = collect_family_info(case_json, case)
        case = case_solved_data(case_json, case)

    elif prog != 'rare_disease':
        print(f"{case_json['case_id']} is a {prog} case, not rare disease.")

    return case


def initialise_case_dict(json):
    """ Set up a dict to hold relevant information for a single case.

    args:
        json [str]: filepath to json dump for a single case

    returns:
        case [dict]: blank dict to hold required case data
    """

    case = {
        'case_id': json['case_id'],
        'program': json['program'],
        'genome': json['assembly'],  # PROBABLY UNNECESSARY
        'no_samples': json['number_of_samples'],  # PROBABLY UNNECESSARY
        'panels': [],
        'affected': {},
        'adopted': {},
        'samples': {},  # PROBABLY UNNECESSARY
        'vcfs': {
            'paths': {},  # -HOPEFULLY- UNNECESSARY
            'initial': {},
            'b38': {},
            'b38_sorted': {}},
        'solved': {
            'fully': False,
            'partially': False},
        'pcvs': {  # potential causal variants identified in original analysis
            'initial': [],
            'b38': []},  # initial pcv notation after liftover
        'reanalysis': {}}

    return case


def get_case_panels(json, case):
    """ Update a case dict with information on which panels were used.

    args:
        json [str]: filepath to json dump for a single case
        case [dict]: dict to hold required case data

    returns:
        case [dict]: updated with panels info
    """

    for panel in json['interpretation_request_data']['json_request'][
        'pedigree']['analysisPanels']:

        case['panels'].append({
            'name': panel['specificDisease'],
            'id': panel['panelName'],
            'version': panel['panelVersion']})

    return case


def collect_family_info(json, case):
    """ Get information on whether each family member is affected or
    adopted, and the IDs of relevant DNA samples and VCFs.

    args:
        json [str]: filepath to json dump for a single case
        case [dict]: dict to hold required case data

    returns:
        case [dict]: updated with family info
    """

    # identify whether the proband is adopted

    adopted = True

    for person in json['interpretation_request_data']['json_request'][
        'pedigree']['members']:

        if person['isProband']:
            if person['adoptedStatus'].lower() == 'notadopted':
                adopted = False

    # get info for (a) proband (b) all family if proband is not adopted

    for person in json['interpretation_request_data']['json_request'][
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

            case['adopted'][relation] = True

            if person['adoptedStatus'] == 'notadopted':
                case['adopted'][relation] = False

            if ('sibling' in relation) and case['adopted'][relation]:
                continue  # omit data for any adopted siblings

            # get affection status

            affect_map = {
                'AFFECTED': 'True',
                'UNCERTAIN': 'Unknown',
                'UNAFFECTED': 'False'}

            case['affected'][relation] = affect_map[person['affectionStatus']]

            # get IDs of all DNA samples

            samples = [s['sampleId'] for s in person['samples']
                if s['product'] == 'DNA']

            if samples:
                case['samples'][relation] = samples

                # get paths to SNV VCFs created from those samples

                vcf_paths = []

                for sample in case['samples'][relation]:

                    vcfs = [vcf['uriFile'] for vcf in json[
                        'interpretation_request_data']['json_request']['vcfs']
                        if (sample in vcf['sampleId'])
                        and (vcf['fileType'] == 'VCF_small')]  # only SNVs

                    for vcf in vcfs:
                        vcf_paths.append(vcf)

                if vcf_paths:
                    assert len(vcf_paths) == 1, \
                        f"{relation} has multiple VCF paths"

                    case['vcfs']['paths'][relation] = vcf_paths

    return case


def case_solved_data(json, case):
    """ Get information on whether the case is solved, and any known or
    potential causal variants.

    args:
        json [str]: filepath to json dump for a single case
        case [dict]: dict to hold required case data

    returns:
        case [dict]: updated with outcome and variant info
    """

    case['pcvs']['initial'] = []

    for report in json['clinicalReports']:  # a case can have multiple reports
        if report['valid']:

            # identify whether case is fully or partially solved
            if report['exitQuestionnaire']:

                solved = report['exitQuestionnaire'][
                    'familyLevelQuestions']['caseSolvedFamily']

                if solved == 'yes':
                    case['solved']['fully'] = True
                    case['solved']['partially'] = True

                elif solved == 'unknown':
                    case['solved']['partially'] = True

            # identify potentially causal variants

            var_types = [
                'variants',
                'structuralVariants',
                'shortTandemRepeats',
                'chromosomalRearrangements']

            for var_type in var_types:

                try:
                    for var in report['clinicalReportData'][var_type]:

                        chrom = var['variantCoordinates']['chromosome']
                        pos = var['variantCoordinates']['position']
                        ref = var['variantCoordinates']['reference']
                        alt = var['variantCoordinates']['alternate']

                        coords = f"{chrom}:{pos}:{ref}:{alt}"

                        if coords not in case['pcvs']['initial']:
                            case['pcvs']['initial'].append(coords)

                except Exception:
                    pass

    return case


""" functions to create bed file """


def get_panelapp_panel(panel_id, panel_version=None):
    """ Retrieve panel object representing specified version of a
    PanelApp panel ('Panelapp.Panel' doesn't always work, because some
    older panel versions don't contain the hgnc_symbol element)

    args:
        panel_id [str/int]: the panel's PanelApp ID

        panel_version [str/float] (optional): version of the panel to use - if
            not supplied, retrieves current version

    returns:
        result: dict of data for specified version of specified panel
    """

    # easiest way is to use the panelapp package

    try:
        result = Panelapp.Panel(str(panel_id), panel_version).get_data()

    # except it doesn't work for some older panel versions

    except Exception:  # can fail for different reasons
        path = ["panels", str(panel_id)]
        param = {"version": panel_version}

        url = api.build_url(path, param)
        result = api.get_panelapp_response(url)

    return result


@time_execution
def create_bed(panels, output_fp):
    """ Using the output of get_case_panels(), create bed files for the
    entities (genes and STRs) covered by the case's panels.

    Bed files are created by querying Ensembl BioMart (GRCh38) with
    lists of HGNC IDs.

    args:
        case_id [str]
        panels [dict]: data for panels used in original analysis
    """

    # get panel genes and regions

    genes = []
    regions = []

    for panel in panels:

        data = get_panelapp_panel(panel['id'])

        print(f"Query: {panel['id']} {panel['name']} {panel['version']}")
        print(f"Retrieved: {data['id']} {data['name']} {data['version']}\n")

        panel_genes = get_panel_genes(data)
        panel_regions = get_panel_regions(data)

        for gene in panel_genes:
            if gene not in genes:
                genes.append(gene)

        for region in panel_regions:
            if region not in regions:
                regions.append(region)

    # read in template biomart query and insert string of HGNC ids

    hgnc_string = ','.join([ele for ele in genes])

    with open('bed_template.txt', 'r') as reader:
        contents = reader.read()

    query = contents.replace('PLACEHOLDER', hgnc_string)

    # query biomart using subprocess

    subprocess.run(['wget', '-O', output_fp, query])

    # add regions into bed file

    region_df = pd.DataFrame(regions, columns=['chrom', 'start', 'end'])

    with open(output_fp, 'r') as reader:
        bed_data = pd.read_csv(reader, sep='\t', header=None)

    bed_data.columns = ['chrom', 'start', 'end']

    with_regions = pd.concat([bed_data, region_df])

    # sort gene/region dataframe by chromosome and position

    sorted = with_regions.sort_values(by=['chrom', 'start'])

    with open(output_fp, 'w') as writer:
        sorted.to_csv(writer, sep='\t', header=True, index=False)


def get_panel_genes(panel_data):
    """ Given data from a PA panel, identify all gene HGNC IDs and
    return as a list.

    args:
        panel_data [dict]: data from a PA panel

    returns:
        genes [list]: list of gene HGNC IDs
    """

    genes = []

    for gene in panel_data['genes']:

        if gene['confidence_level'] == '3' and \
            gene['gene_data']['hgnc_id'] and \
            (gene['gene_data']['hgnc_id'] not in genes):

            genes.append(gene['gene_data']['hgnc_id'])

    return genes


def get_panel_regions(panel_data):
    """ Given data from a PA panel, identify all standalone genomic
    regions covered by that panel and return as a list.

    args:
        panel_data [dict]: data from a PA panel

    returns:
        regions [list]: of lists, each with form ['chrom', 'start', 'end']
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


""" functions to standardise vcf files """


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


def check_vcf_pcvs(vcf):
    """ Check which variants of interest are present in a VCF file
    using their INFO field value of 'PCV=1'.

    args:
        vcf [fp]

    returns:
        pcvs_present [list]: in format chrom:pos:ref:alt
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
def lift_over_vcf(vcf_in, vcf_out, genome_file, chain_file, pcvs_initial):
    """ Liftover a b37 VCF file to b38 using CrossMap.

    args:
        vcf_in [fp]
        vcf_out [fp]
        genome_file [fp]: reference genome fasta, required for liftover
        chain_file [fp]: required for liftover
        pcvs_initial [list]: variants as chrom:pos:ref:alt (b37 notation)

    returns:
        pcvs_b38 [list]: variants as chrom:pos:ref:alt (b38 notation)
    """

    # lift b37 vcf file over to b38

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

    # return a list of the PCVs present in the output file

    pcvs_b38 = []

    if pcvs_initial:
        pcvs_b38 = check_vcf_pcvs(vcf_out)

    return pcvs_b38


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
def sort_vcf(vcf_file, output_fp):
    """ Use bcftools to sort a VCF file.

    args:
        vcf_file [str]: input VCF file
        output_fp [str]: path to sorted output file
        pcvs_b38 [list]: in format chrom:pos:ref:alt
    """

    subprocess.run([
        'bcftools',
        'sort',
        '-o', output_fp,
        vcf_file,
        ])


""" functions to annotate variants """


@time_execution
def cbtools_annotation(input_fp, output_fp):
    """ Use the CellBase API to annotate a VCF file.

    Args:
        vcf_file [str]: input VCF file
        output_fp [str]: path to sorted output file
    """

    # according to the pycellbase pypi page
    # subprocess.run([ 'cbtools.py', 'annotation', input_fp, '>', output_fp])

    # according to cbtools.py annotate -h

    subprocess.run([
        'cbtools.py', 'annotate',
        '--config', 'config.json',
        '-o', output_fp,
        input_fp])


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


@time_execution
def apply_filters(variants, parameters):
    """ Given annotated variant data and a set of filtering parameters,
    identify which variants meet filtering thresholds.

    args:
        variants [list]: 1 dict of annotation per variant
        parameters [dict]: values to use for each filter

    returns:
        output_variants [dict]: subset of 'variants' which passed filtering
    """

    output_variants = {}
    lof_vars = ['stop_gained', 'start_lost']
    exclude = ['intron_variant', '5_prime_UTR_variant',
        'non_coding_transcript_exon_variant']

    for variant in variants:

        variant['filters'] = {
            'affected_transcripts': False,
            'af': False,
            'cadd': False,
            'has_vda': False,
            'constraint': False,
            'passed_filters': 0}

        # does the variant affect any transcripts

        if variant['affected_transcripts']:
            variant['filters']['affected_transcripts'] = True
            variant['filters']['passed_filters'] += 1

        # compare allele frequency to threshold

        if not variant['allele_frequency'] or \
            (variant['allele_frequency'] and
            (variant['allele_frequency'] <= parameters['max_af'])):

            variant['filters']['af'] = True
            variant['filters']['passed_filters'] += 1

        # compare CADD score to threshold

        if variant['scaled_cadd'] and \
            (variant['scaled_cadd'] >= parameters['min_cadd']):

            variant['filters']['cadd'] = True
            variant['filters']['passed_filters'] += 1

        # does the variant have any DisGeNET variant-disease associations
        # for info only - not currently used for filtering

        if variant['vda']:
            variant['filters']['has_vda'] = True

        # compare constraint to relevant thresholds
        # for info only - not currently used for filtering

        consq = variant['consequence']
        constraint = variant['constraint']

        if consq and constraint:

            if (consq == 'synonymous_variant') and \
                (constraint['oe_syn'] <= parameters['max_oe_syn']):

                variant['filters']['constraint'] = True
                variant['filters']['passed_filters'] += 1

            elif (consq == 'missense_variant') and \
                (constraint['oe_mis'] <= parameters['max_oe_mis']):

                variant['filters']['constraint'] = True
                variant['filters']['passed_filters'] += 1

            elif (consq in lof_vars) and \
                ((constraint['oe_lof'] <= parameters['max_oe_lof']) or
                (constraint['exac_pLI'] >= parameters['min_pli'])):

                variant['filters']['constraint'] = True

        # add variant to output depending on which filters passed

        if (consq not in exclude) and \
            (variant['filters']['af']) and \
            (variant['filters']['cadd']) and \
            (variant['filters']['affected_transcripts']):

            output_variants[variant['variant']] = variant

    return output_variants


@time_execution
def create_output(case, person, parameters, output_fp):
    """ From a list of annotated and filtered variants identified in a
    specific family member, create an output text file describing the results.

    args:
        case [dict]
        person [str]
        parameters [dict]: parameters used for filtering
        output_fp [fp]
    """

    pcvs_initial = case['pcvs']['initial']
    pcvs_b38 = case['pcvs']['b38']
    pcvs_anno = case['pcvs']['annotated']

    output = case['reanalysis'][person]
    output_strs = [output[v]['variant'] for v in case['reanalysis'][person].keys()]

    original = len([v for v in output_strs if v in pcvs_initial])
    novel = len([v for v in output_strs if v not in pcvs_initial])

    with open(output_fp, 'w') as writer:

        writer.write(f"{dt.now()}\n\nResults for {person}\n"
            f"\tAffected: {case['affected'][person]}\n"
            f"\tAdopted: {case['adopted'][person]}\n"
            f"\tCase solved:\n"
            f"\t\tFully: {case['solved']['fully']}\n"
            f"\t\tPartially: {case['solved']['partially']}\n\n"
            "Parameters applied:\n")

        for key, value in parameters.items():
            writer.write(f"\t{key}: {value}\n")

        writer.write(f"\n{len(output_strs)} variants passed filtering\n"
            f"\t{original} identified in the original analysis\n"
            f"\t{novel} not identified in the original analysis\n")

        if [v for v in pcvs_initial if v in output_strs] == pcvs_initial:
            writer.write(
                f"\n\tAll originally identified SNVs passed filtering\n\n")
        else:
            writer.write(
                f"\n\tNot all originally identified SNVs passed filtering\n\n")

        writer.write('\nVARIANTS IDENTIFIED IN ORIGINAL ANALYSIS\n'
            'Kept after filtering:\n\n')

        for var_dict in pcvs_anno:
            if var_dict['variant'] in output_strs:
                for key, value in var_dict.items():
                    writer.write(f"{key}: {value}\n")
                writer.write('\n')

        writer.write('\nExcluded after filtering:\n\n')

        for var_dict in pcvs_anno:
            if var_dict['variant'] not in output_strs:
                for key, value in var_dict.items():
                    writer.write(f"{key}: {value}\n")
                writer.write('\n')

        writer.write(f"\n\nVARIANTS NOT IDENTIFIED IN ORIGINAL ANALYSIS\n\n")

        for var in output_strs:
            if var not in pcvs_b38:
                for key, value in output[var].items():
                    writer.write(f"{key}: {value}\n")
                writer.write('\n')


@time_execution
def filter_on_segregation(case):
    """ Compare the VCFs of other family members against that of the
    proband to identify variants which segregate in accordance with the
    phenotype.

    args:
        case [dict]

    returns:
        compare [fp]: list of variants which segregate with phenotype
    """

    # ALSO NEEDS TO TAKE PENETRANCE INTO ACCOUNT

    compare = case['vcfs']['b38_sorted']['proband']
    i = 0

    for person, vcf in case['vcfs']['b38_sorted'].items():
        if person != 'proband':

            # if person is affected and not adopted, get intersect
            if (case['affected'][person] == 'True') and \
                (case['adopted'][person] == 'False'):

                output = f"outside_github/vcfs/{case['case_id']}" \
                    f"_4_segregation_{i}.vcf"

                subprocess.run(
                    ['bedtools', 'intersect',
                    '-a', compare, '-b', vcf, '>', output])

                compare = output

            # if person is NOT affected, get proband vars not seen in person
            elif case['affected'][person] == 'False':

                output = f"outside_github/vcfs/{case['case_id']}" \
                    f"_4_segregation_{i}.vcf"

                subprocess.run(
                    ['bedtools', 'intersect',
                    '-a', compare, '-b', vcf, '-v', '>', output])

                compare = output

            i += 1

    return compare


@time_execution
def filter_de_novo(case):
    """ Compare the VCFs of both parents against that of the proband to
    identify de novo variants.

    args:
        case [dict]

    returns:
        compare [fp]: list of variants present in proband but not parents
    """

    # ALSO NEEDS TO TAKE PENETRANCE INTO ACCOUNT

    compare = case['vcfs']['b38_sorted']['proband']
    parents = ['mother', 'father']
    i = 0

    for person, vcf in case['vcfs']['b38_sorted'].items():
        if person in parents:

            output = f"outside_github/vcfs/{case['case_id']}" \
                f"_4_denovo_{i}.vcf"

            subprocess.run(
                ['bedtools', 'intersect',
                '-a', compare, '-b', vcf, '-v', '>', output])

            compare = output
            i += 1

    return compare


def run_whole_case(case_id, chrom_map, genome_file, chain_file, parameters):
    """  """

    json_fp = f'outside_github/inputs/{case_id}.json'
    vcf_initial = f'outside_github/inputs/{case_id}'
    vcf_prefix = f'outside_github/intermediates/{case_id}'
    output_dir = 'outside_github/outputs/'

    # parse case json for required data

    case = process_case_json(json_fp)
    assert case, f"Error parsing data for case {case}"

    case['pcvs']['initial'].sort()
    pcvs_initial = case['pcvs']['initial']

    assert 'proband' in case['vcfs']['paths'].keys(), \
        'Error: Proband VCF path required'

    bed_file = f'outside_github/bed_files/{case_id}.bed'
    create_bed(case['panels'], bed_file)

    # obtain a copy of each VCF from its path

    for person, vcf in case['vcfs']['paths'].items():

        assert len(vcf) == 1, f'{person} has multiple VCFs'

        if vcf[0].endswith('.vcf.gz'):
            local_path = f'{vcf_initial}_1_{person}.vcf.gz'
        elif vcf[0].endswith('.vcf'):
            local_path = f'{vcf_initial}_1_{person}.vcf'

        case['vcfs']['initial'][person] = local_path

    assert 'proband' in case['vcfs']['initial'].keys(), \
        'Error: Proband initial VCF required'

    # standardise CHROM notation and reference genome across VCFs

    for person, vcf in case['vcfs']['initial'].items():

        fixed_vcf = f'{vcf_prefix}_2_fixed_{person}.vcf.gz'
        marked_vcf = f'{vcf_prefix}_3a_marked_{person}.vcf.gz'
        b38_vcf = f'{vcf_prefix}_3b_marked_b38_{person}.vcf.gz'

        fix_chroms(vcf, fixed_vcf, chrom_map)
        fixed_pcvs = mark_pcvs(fixed_vcf, marked_vcf, pcvs_initial)

        if person == 'proband':
            assert fixed_pcvs == pcvs_initial, \
                f'Error: PCVs lost during fixing\n' \
                f'Initial PCVS: {pcvs_initial}\n' \
                f'Retained after fixing: {fixed_pcvs}'

        genome = identify_genome(marked_vcf)

        if genome == 'GRCh37':
            lifted_pcvs = lift_over_vcf(
                marked_vcf, b38_vcf, genome_file, chain_file, pcvs_initial)

            if person == 'proband':
                assert len(lifted_pcvs) == len(pcvs_initial), \
                    f'Error: PCVs lost during liftover\n' \
                    f'Initial PCVS: {pcvs_initial}\n' \
                    f'Retained after liftover: {lifted_pcvs}'

                case['pcvs']['b38'] = lifted_pcvs

        elif genome == 'GRCh38':
            b38_vcf = marked_vcf

            if person == 'proband':
                case['pcvs']['b38'] = pcvs_initial
                pcvs_b38 = case['pcvs']['b38']

        case['vcfs']['b38'][person] = b38_vcf

    assert len(pcvs_b38) == len(pcvs_initial), \
        f'Error: Different number of b37 and b38 PCVs\n' \
        f'Initial PCVs: {pcvs_initial}\nb38 PCVS: {pcvs_b38}'

    assert 'proband' in case['vcfs']['b38'].keys(), \
        'Error: Proband b38 VCF required'

    # filter, sort, annotate and prioritise each b38 vcf

    for person, vcf in case['vcfs']['b38'].items():

        filter_prefix = f'{vcf_prefix}_4_filtered_{person}'
        sorted = f'{vcf_prefix}_5_sorted_{person}.vcf'
        case['vcfs']['b38_sorted'][person] = sorted
        output_fp = f'{output_dir}{case_id}_{person}.txt'

        filter_1 = filter_on_bed(
            vcf, bed_file, parameters['qual'], filter_prefix)

        sort_vcf(filter_1, sorted)

        # saving a VCF's annotated variants reduces API calls in testing
        json_vars = f'outside_github/intermediates/{case_id}_6_{person}_ann.json'
        try:
            with open(json_vars, 'r') as reader:
                anno = json.load(reader)
        except FileNotFoundError:
            anno = cb_client_annotation(sorted)
            with open(json_vars, 'w') as writer:
                json.dump(anno, writer)

        anno_less, anno_pcvs = condense_annotation(anno, pcvs_b38)
        case['pcvs']['annotated'] = anno_pcvs
        anno_filt = apply_filters(anno_less, parameters)

        if person == 'proband':
            pcvs_final = [anno_filt[v]['variant'] for v in anno_filt.keys()]
            pcvs_final.sort()
            print(f'All initial PCVs: {pcvs_initial}\n'
                f'b38 PCVs: {pcvs_b38}\n'
                f'Final PCVs in proband: {pcvs_final}\n')

        case['reanalysis'][person] = anno_filt
        create_output(case, person, parameters, output_fp)

    # # filter based on available family vcfs

    # family_vcfs = {
    #     'proband': case['vcfs']['b38_sorted']['proband'],
    #     'segregated': None,
    #     'de_novo': None}

    # if case['reanalysis'].keys() != ['proband']:

    #     family_vcfs['segregated'] = filter_on_segregation(case)

    #     if ('mother' in case['reanalysis'].keys()) and \
    #         ('father' in case['reanalysis'].keys()):

    #         family_vcfs['de_novo'] = filter_de_novo(case)

    return case


def run_single_vcf(
    case_id, person, chrom_map, genome_file, chain_file, parameters):
    """  """

    input_vcf = f'outside_github/inputs/{case_id}_1_{person}.vcf.gz'
    fixed_vcf = f'outside_github/intermediates/{case_id}_2_fixed_{person}.vcf.gz'
    marked_vcf = f'outside_github/intermediates/{case_id}_2_marked_{person}.vcf.gz'
    b38_vcf = f'outside_github/intermediates/{case_id}_2_marked_b38_{person}.vcf.gz'
    prefix = f'outside_github/intermediates/{case_id}_3_filtered_{person}'
    sorted_vcf = f'outside_github/intermediates/{case_id}_4_sorted_{person}.vcf'
    var_json = f'outside_github/outputs/{case_id}_output_vars_{person}.json'
    output_fp = f'outside_github/outputs/{case_id}_{person}.txt'

    case = process_case_json(f'outside_github/inputs/{case_id}.json')
    pcvs = case['pcvs']['initial']

    fix_chroms(input_vcf, fixed_vcf, chrom_map)
    fixed_pcvs = mark_pcvs(fixed_vcf, marked_vcf, pcvs)
    genome = identify_genome(fixed_vcf)

    if genome == 'GRCh37':
        pcvs_38 = lift_over_vcf(
            fixed_vcf, b38_vcf, genome_file, chain_file, pcvs)

    else:
        pcvs_38 = pcvs
        b38_vcf = fixed_vcf

    filter_vcf = filter_on_bed(b38_vcf,
        f'outside_github/bed_files/{case_id}.bed', '20', prefix)

    sort_vcf(filter_vcf, sorted_vcf)

    variants = cb_client_annotation(sorted_vcf)
    with open(var_json, 'w') as writer:
        json.dump(variants, writer)

    # with open(var_json, 'r') as reader:
    #     variants = json.load(reader)

    condensed = condense_annotation(variants)
    output = apply_filters(condensed, parameters)

    if person == 'proband':
        case['pcvs']['final'] = [v['variant'] for v in output]

    create_output(case, person, output, parameters, output_fp)


def main():
    """ file paths """

    all_case_ids_file = 'all_case_ids.txt'

    cases_with_vcfs = [
        'outside_github/inputs/SAP-55997-1.json',  # b38 solved (2), p & twin (A)
        'outside_github/inputs/SAP-56251-1.json']  # b38 solved (2), p & parents (NA)

    cases_no_vcfs = [
        'outside_github/inputs/SAP-48034-1.json',  # unsolved, no vars
        'outside_github/inputs/SAP-56130-1.json',  # unsolved, no vars
        'outside_github/inputs/SAP-56069-1.json',  # unsolved, 1 var?
        'outside_github/inputs/SAP-56172-1.json']  # partially solved, v1, 9 vars

    gen_dir = 'outside_github/genomes/'
    genome_file = f'{gen_dir}GCF_000001405.40_GRCh38.p14_genomic.fna'
    chain_file = f'{gen_dir}hg19ToHg38.over.chain.gz'

    chrom_map = 'chrom_map.txt'

    """ define filtering parameters """

    parameters = {
        'qual': '20',  # for initial QUAL value filtering
        'max_af': 0.05,  # seen in <=5% of the population (gnomAD)
        'min_cadd': 10,  # CADD >= 10 implies p(var not observed) >= 0.9
        'min_pli': 0.9,  # higher implies less tolerant to truncation (gnomAD)
        'max_oe_lof': 0.35,  # lower implies less tolerant to variants (gnomAD)
        'max_oe_mis': 0.35,
        'max_oe_syn': 0.35,
        'min_vda': 0.75}  # score=0.75 ~= 2+ curated sources (DisGeNET)

    """ run the pipeline """

    case_id = 'SAP-55997-1'
    run_whole_case(case_id, chrom_map, genome_file, chain_file, parameters)
    # run_single_vcf(
    #     case_id, person, chrom_map, genome_file, chain_file, parameters)

    """ end """


if __name__ == '__main__':
    main()
