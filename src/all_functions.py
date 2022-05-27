#!usr/bin/env python

"""
oc.<client>.search()

-gets data on ALL instances in that client
-e.g. oc.projects.search() gets data on all projects

oc.<client>.info(<client> = <client instance id>)

-gets data on a SUBSET of instances in that client
-e.g. oc.projects.info(projects = <project id>) gets data from one project
"""


import datetime as dt
import pandas as pd
import subprocess

from bs4 import BeautifulSoup as bs

from panelapp import api
from panelapp import Panelapp
from panelapp import queries

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpencgaClient


""" Generic functions """


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

    config_dict = {'rest': {'host': host}}
    config = ClientConfiguration(config_dict)

    oc = OpencgaClient(config)
    oc.login(user)

    token = oc.token
    oc = OpencgaClient(configuration=config, token=token)

    return oc


def get_hgnc_id(hgnc_df, gene_symbol):
    """ Get the HGNC ID for a supplied gene symbol, if one exists.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        gene_symbol [str]: any gene symbol

    returns:
        hgnc_id [str], or None: HGNC ID of gene associated with entity
    """

    hgnc_id = None

    # If a row exists where this gene is official gene symbol,
    try:

        # Get that row's index and associated HGNC ID
        target_index = hgnc_df.index[hgnc_df['symbol'] == gene_symbol]
        hgnc_id = hgnc_df.loc[target_index[0], 'hgnc_id']

    # If this gene isn't an official symbol,
    except IndexError:

        i = 0

        # Look through 'previous official gene symbols' column
        for value in hgnc_df['prev_symbol']:

            # If gene appears in this column, get that row's HGNC ID
            if gene_symbol in str(value):

                hgnc_id = hgnc_df.iloc[i].loc['hgnc_id']
                break

            i += 1

    return hgnc_id


def print_dict_keys(some_dict):
    """ For a supplied dict, print each key and the data type of its
    value """

    for key, value in some_dict.items():
        print('{} [{}]'.format(key, type(value)))


def write_dict_contents(some_dict, filename):
    """ Write out entire contents of a dict to specified filename """

    with open(filename, 'w') as writer:

        for key, value in some_dict.items():
            writer.write('{}: {}\n'.format(key, value))


""" Bed file functions """


def get_hgnc_string(entity_list):
    """

    args:
        entity_list [list]:

    returns:
        hgnc_string [str]:
    """

    hgnc_string = ''

    for entity in entity_list:
        if entity['associated_gene']:

            hgnc = entity['associated_gene'][6:]

            if entity != entity_list[-1]:
                hgnc_string += '{},'.format(hgnc)

            else:
                hgnc_string += '{}'.format(hgnc)

    return hgnc_string


def get_biomart_output(case_id, panel_list):
    """ Using the output of get_case_panels(), create bed files for the
    entities (genes and STRs) covered in the original and current
    versions of the panels in the case.

    Bed files are created by querying Ensembl BioMart (GRCh37) with
    lists of HGNC IDs.

    Can't currently include regions in bed files, because they appear to
    have neither an associated gene nor a GRCh37 location...

    args:
        case_id: ID of the OpenCGA case (for file naming purposes)

        panel_list: output of get_case_panels() - each element is a dict
            containing details on one panel (id, name, original &
            current versions, original & current entities)

    returns:
        list of bed file filename strings
    """

    # given an HGNC ID, how does BioMart turn that into a bed file region?
    # does the region include all exons and introns? what about padding?

    # get lists of all entities in original and current versions of panels

    original = []
    current = []

    for panel in panel_list:

        for entity in panel['original_entities']:
            if entity not in original:
                original.append(entity)

        for entity in panel['current_entities']:
            if entity not in current:
                current.append(entity)

    # define destinations for output bed files

    arg_dicts = [
        {'input' : original,
        'output' : 'data/bed_files/case_{}_original.bed'.format(case_id)},

        {'input' : current,
        'output' : 'data/bed_files/case_{}_current.bed'.format(case_id)}
        ]

    # query BioMart API using subprocess for each list of entities

    for dict in arg_dicts:

        hgnc_string = get_hgnc_string(dict['input'])
        filename = dict['output']

        query = (
            'https://grch37.ensembl.org/biomart/martservice?query=<?xml '
            'version="1.0" encoding="utf-8"?><!DOCTYPE Query>'
            '<Query count="" datasetConfigVersion="0.6" formatter="TSV" '
            'header="0" uniqueRows="0" virtualSchemaName="default">'
            '<Dataset interface="default" name="hsapiens_gene_ensembl">'
            '<Filter name="hgnc_id" value="{}"/>'
            '<Attribute name="chromosome_name"/>'
            '<Attribute name="start_position"/>'
            '<Attribute name="end_position"/>'
            '</Dataset></Query>'.format(hgnc_string))

        biomart_call = [
            'wget',
            '-O',
            '{}'.format(filename),
            '{}'.format(query)
            ]

        subprocess.run(biomart_call)

    return arg_dicts[0]['output'], arg_dicts[1]['output']


def write_bed_query_template(template_file, case_id, gene_string):
    """ Construct an XML query file to supply to BioMart for a specific case

    args:
        template_file: name of XML file containing template BioMart query
        case_id: ID of OpenCGA case
        gene_string: comma-separated list of HGNC IDs (lacking 'HGNC:' prefix)

    returns:
        output_file: path to XML BioMart query file for specified case
    """

    with open(template_file, 'r') as reader:
        template_contents = reader.read()

    soup = bs(template_contents, 'xml')

    # Replace empty tag attribute with list of HGNC IDs
    soup.Filter['value'] = gene_string

    output_file = 'data/case_{}_bed_query.xml'.format(case_id)

    with open(output_file, 'w') as writer:
        writer.write(str(soup))

    return output_file


""" PanelApp or PanelApp panel object functions """


def write_current_pa_list():
    """ Get a list of IDs of currently existing PanelApp panels """

    panels = queries.get_all_panels()  # gets current panel ids & objects
    current_panels = sorted([key for key in panels.keys()])  # gets panel ids

    with open('data/refs/current_pa_panel_ids.txt', 'w') as writer:

        for panel_id in current_panels:
            writer.write('{}\n'.format(panel_id))

    return current_panels


def get_panelapp_panel(id, version = None):
    """ Retrieve panel object representing specified version of a
    PanelApp panel ('Panelapp.Panel' doesn't always work, because some
    older panel versions don't contain the hgnc_symbol element)

    args:
        id [str/int]: the panel's PanelApp ID

        version [str/float] (optional): version of the panel to use - if
            not supplied, retrieves current version

    returns:
        panel: dict of PanelApp data on specified version of panel
        None: if id does not correspond to a current PA panel
    """

    if version:
        path = ["panels", str(id)]
        param = {"version": str(version)}

        url = api.build_url(path, param)
        panel = api.get_panelapp_response(url)

        return panel

    else:
        # If version not supplied, check panel currently exists
        with open('data/refs/current_pa_panel_ids.txt') as reader:
            contents = reader.readlines()

        current_panels = [line.strip() for line in contents]

        if id in current_panels:
            panel = Panelapp.Panel(str(id)).get_data()

            return panel

        else:
            print('{} is not a current PanelApp panel ID')

            return None


def write_pa_panel(id, version = None):
    """ Dump out the contents of a PanelApp panel to a text file

    args:
        id [str/int]: PanelApp ID for a panel
        version [str/float] (optional): specific panel version to get
    """

    panel = get_panelapp_panel(id, version)

    filename = 'pa_dump_panel_{}_version_{}.txt'.format(id, panel['version'])

    with open(filename, 'w') as writer:
        writer.write(str(panel))


def get_all_green_entities(hgnc_df, id, version = None):
    """ Get details of all 'green' genes, regions and STRs (i.e.
    confidence level 3) from a PanelApp panel.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        id [str/int]: PanelApp ID for a panel
        version [str/float] (optional): panel version, defaults to current

    returns:
        panel_version [str]: version of panel in PA
        panel_entities [list of dicts]: green genes, regions and STRs
    """

    try:
        panel = get_panelapp_panel(id, version)

        panel_version = panel['version']  # in case version not supplied as arg
        panel_entities = []

        for type in ['genes', 'regions', 'strs']:
            for entity in panel['{}'.format(type)]:
                if entity['confidence_level'] == '3':

                    info = get_entity_dict(hgnc_df, entity)
                    panel_entities.append(info)

        return panel_version, panel_entities

    except AttributeError as error:

        print('Error for panel {}: {}'.format(id, error))

        return None, None


def get_green_genes(id, version = None):
    """ Get a list of the green genes in a PanelApp panel

    args:
        hgnc_df: pandas df containing dump of HGNC site
        panel_object: PanelApp panel object

    returns:
        panel_genes: list of dicts, one for each green panel gene
            {entity_type,
            entity_name,
            associated_gene,
            chromosome,
            grch37_coordinates}
    """

    panel = get_panelapp_panel(id, version)

    panel_genes = []

    for gene in panel['genes']:
        if gene['confidence_level'] == '3':

            gene_symbol = gene['gene_symbol']

            if gene_symbol not in panel_genes:
                panel_genes.append(gene_symbol)

    return panel_genes


def compare_green_genes(id_1, id_2, version_1 = None, version_2 = None):
    """ Compare the lists of green genes in two PanelApp panel objects
    (can be two versions of the same panel)

    args:
        pa_panel_1: dict representing a PanelApp panel
        pa_panel_2: dict representing a PanelApp panel

    returns:
        bool for 'the two gene lists are identical'
    """

    panel_1 = get_panelapp_panel(id_1, version_1)
    panel_2 = get_panelapp_panel(id_2, version_2)

    gene_lists = []

    for panel in [panel_1, panel_2]:

        genes = get_green_genes(panel)
        genes.sort()
        gene_lists.append(genes)

    if gene_lists[0] == gene_lists[1]:
        return True

    else:
        return False


def get_entity_dict(hgnc_df, entity):
    """ Get a dict of details about a gene/region/STR from a PA panel

    args:
        hgnc_df: pandas df containing dump of HGNC site
        entity [dict]: single gene/region/STR from a PA panel object

    returns:
        info [dict]: details about that entity
            {type,
            name,
            chromosome,
            grch37_coordinates,
            associated_gene}
    """

    # Initialise output dict

    info = {'type' : entity['entity_type']}

    # Get human-readable entity name (different for regions)

    if entity['entity_type'] == 'region':
        info['name'] = entity['verbose_name']

    else:
        info['name'] = entity['entity_name']

    # Get entity's genomic coordinates (different for genes)

    if entity['entity_type'] == 'gene':

        try:
            location = entity['gene_data']['ensembl_genes']['GRch37']['82'][\
                'location']

            split_location = location.split(':')
            split_coordinates = split_location[1].split('-')

            info['chromosome'] = split_location[0]

            info['grch37_coordinates'] = [
                int(split_coordinates[0]),
                int(split_coordinates[1])]

        except KeyError:

            info['chromosome'] = None
            info['grch37_coordinates'] = None

    else:
        info['chromosome'] = entity['chromosome']
        info['grch37_coordinates'] = entity['grch37_coordinates']

    # Get HGNC ID for entity's associated gene (if there is one)

    if entity['gene_data'] != None:

        gene_symbol = entity['gene_data']['gene_symbol']

        info['associated_gene'] = get_hgnc_id(hgnc_df, gene_symbol)

    else:
        info['associated_gene'] = None

    return info


def get_entity_hgnc(hgnc_df, gene_symbol):
    """ Get the HGNC IF for a supplied gene symbol, if one exists.

    args:
        hgnc_df: pandas df containing dump of HGNC site
        gene_symbol [str]: any gene symbol

    returns:
        hgnc_id [str], or None: HGNC ID of gene associated with entity
    """

    hgnc_id = None

    # If a row exists where this gene is official gene symbol,
    try:

        # Get that row's index and associated HGNC ID
        target_index = hgnc_df.index[hgnc_df['symbol'] == gene_symbol]
        hgnc_id = hgnc_df.loc[target_index[0], 'hgnc_id']

    # If this gene isn't an official symbol,
    except IndexError:

        i = 0

        # Look through 'previous official gene symbols' column
        for value in hgnc_df['prev_symbol']:

            # If gene appears in this column, get that row's HGNC ID
            if gene_symbol in str(value):

                hgnc_id = hgnc_df.iloc[i].loc['hgnc_id']
                break

            i += 1

    return hgnc_id


""" Project-level functions """


def print_all_projects(oc):
    """ List id of all OpenCGA projects to which user has access """

    all_projects = oc.projects.search().get_results()

    for project in all_projects:
        print(project['id'])


def print_project_studies(oc, project_id):
    """ List ids of all studies in a single project """

    project_studies = oc.studies.search(project = project_id).get_results()

    for study in project_studies:
        print(study['id'])


""" Study-level functions """


def get_single_case_data(oc, study_id, case_id):
    """ Retrieve all data from specified case in specified study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
        case_id [str]: id specifying case within that study

    returns:
        case: dict containing all case data
    """

    case = oc.clinical.info(
        study = study_id,
        clinical_analysis = case_id
        ).get_results()[0]

    return case


def get_multi_case_data(oc, study_id, case_limit):
    """ Retrieve the first <case_limit> cases from a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
        case_id [str]: id specifying case within that study

    returns:
        cases: list of case dicts
    """

    cases = oc.clinical.search(
        study = study_id,
        limit = case_limit
        ).get_results()

    return cases


def get_study_case_count(oc, study_id):
    """ Print the number of cases in a study
    (cases are huge, so only retrieves their IDs to count them)

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        len(cases): number of cases in the study
    """

    cases = oc.clinical.search(
        study = study_id,
        include = 'id',
        ).get_results()

    print('{} contains {} cases'.format(study_id, len(cases)))

    return len(cases)


def get_case_list(oc, study_id):
    """ Get a list of all case IDs in a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        case_ids [list]: list of case ID strings
    """

    cases = oc.clinical.search(study = study_id, include = 'id').get_results()

    case_ids = [case['id'] for case in cases]
    case_ids.sort()

    return case_ids


def write_study_case_ids(oc, study_id):
    """ Write out the IDs of all cases within a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
    """

    case_ids = get_case_list(oc, study_id)

    with open('all_case_ids.txt', 'w') as writer:
        for id in case_ids:
            writer.write(id + '\n')


def get_study_file_count(oc, study_id):
    """ Print the number of files in a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        len(files): number of files in the study
    """

    files = oc.files.search(
        study = study_id,
        include = 'id',
        ).get_results()

    print('{} contains {} files'.format(study_id, len(files)))

    return len(files)


def get_study_individual_count(oc, study_id):
    """ Print the number of individuals in a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        len(individuals): number of individuals in the study
    """

    individuals = oc.individuals.search(
        study = study_id,
        include = 'id'
        ).get_results()

    print('{} contains {} individuals'.format(study_id, len(individuals)))

    return len(individuals)


def get_study_sample_count(oc, study_id):
    """ Print the number of samples in a study

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        len(samples): number of samples in the study
    """

    samples = oc.samples.search(
        study = study_id,
        include = 'id',
        ).get_results()

    print('{} contains {} samples'.format(study_id, len(samples)))

    return len(samples)


def check_panels_all_cases(oc, study_id):
    """ Iterate over every case in a study and check that the panels
    covered in 'case > panels' are the same as those covered by
    'interpretation > panels' across all interpretations

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study

    returns:
        all_panels_match [bool]: True/False for whether all cases
            satisfy "'case > panels' and 'interpretations > panels' are
            equivalent"
    """

    case_ids = get_case_list(oc, study_id)

    case_count = len(case_ids)
    matching_cases = 0

    for case_id in case_ids:

        case = get_single_case_data(oc, study_id, case_id)
        check = check_single_case_panels(case)

        if check == True:
            matching_cases += 1

    if matching_cases == case_count:  # if panels match for all cases
        print('Case panels match interpretation panels for all cases')

    else:
        mismatches = case_count - matching_cases

        print('Panel mismatches in {} of {} cases'.format(
            mismatches,
            case_count))


def get_grs_cases(oc, study_id, grs_panels):
    """ Given a list of panels which cover at least one gene, region and
    STR each, return a list of cases in the specified study which
    include at least one of those panels.

    args:
        oc: OpenCGa client
        study_id [str]: defines OpenCGA project and study to use
        grs_panels [filename]: newline-separated list of panel IDs

    returns:
        grs_cases [list]:
    """

    cases = oc.clinical.search(
        study = study_id,
        include = 'id,panels',
        ).get_results()

    with open(grs_panels) as reader:
        id_list = reader.readlines()

    grs_panel_ids = sorted([id.strip() for id in id_list])
    grs_cases = []

    for case in cases:
        for panel in case['panels']:

            panel_id = ((panel['id']).split('-'))[1].strip()

            if (panel_id in grs_panel_ids) and (case['id'] not in grs_cases):

                grs_cases.append(case['id'])

    grs_cases.sort()

    return grs_cases


def print_smallest_multipanel_grs_case(oc, study_id, grs_panels):
    """ Given a list of panels which cover at least one gene, region and
    STR each, identify cases in the specified study which include at
    least one of those panels, and print the ID of the smallest.

    args:
        oc: OpenCGa client
        study_id [str]: defines OpenCGA project and study to use
        grs_panels [filename]: newline-separated list of panel IDs
    """

    grs_cases = get_grs_cases(oc, study_id, grs_panels)

    smallest_count = 10000
    smallest_case = ''

    for case in grs_cases:

        data = get_single_case_data(oc, study_id, case)
        panel_count = len(data['panels'])

        if panel_count > 1:

            count = get_case_entity_count(data)

            if count < smallest_count:

                smallest_count = count
                smallest_case = case

    print('{} has {} entities'.format(smallest_case, smallest_count))


""" Case-level functions """


def get_all_case_intpns(case):
    """ Return a list of all interpretation dicts in a case

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        all_intpns: list of dicts, each dict is one case interpretation
    """

    all_intpns = []

    all_intpns.append(case['interpretation'])

    for interpretation in case['secondaryInterpretations']:
        all_intpns.append(interpretation)

    return all_intpns


def get_case_panels(case, hgnc_dump):
    """ Retrieve details about each panel used in the specified case.
    Includes panel ID and name; version of original panel; list of all
    genes, regions and STRs with confidence level 3 in the original
    version; and the same for the current version.

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        panels_list [list]: list of panel dicts-
            {'panel_id',
            'panel_name',
            'original_version',
            'original_entities',
            'current_version',
            'current_entities'}
    """

    panels_list = []

    for panel in case['panels']:

        id = panel['source']['id']
        name = panel['name']

        original_version, original_entities = get_all_green_entities(
            hgnc_dump,
            id,
            panel['source']['version'])

        current_version, current_entities = get_all_green_entities(
            hgnc_dump,
            id)

        panel_info = {
            'panel_id' : id,
            'panel_name' : name,
            'original_version' : original_version,
            'original_entities' : original_entities,
            'current_version' : current_version,
            'current_entities' : current_entities
            }

        panels_list.append(panel_info)

    return panels_list


def get_case_phenotypes(case):
    """ Return a list of the phenotypes associated with a case.

    A case may have multiple interpretations. Each interpretation has
    its own set of variants. Each variant can have multiple 'evidences',
    each of which may be associated with multiple phenotypes. Finding
    all the phenotypes associated with a case therefore requires
    iterating over all of its intepretations, variants, and evidences.

    args:
        case [dict]: all data from a single case in an OpenCGA study

    returns:
        phenotypes [list]: each element is a str of one phenotype
    """

    phenotypes = []
    interpretations = get_all_case_intpns(case)

    for intpn in interpretations:
        for variant in intpn['primaryFindings']:
            for evidence in variant['evidences']:
                for phenotype in evidence['phenotypes']:

                    phen_string = phenotype['id']
                    phen_list = phen_string.split(',')

                    for phen in phen_list:
                        if phen.strip() not in phenotypes:

                            phenotypes.append(phen.strip())

    return phenotypes



def check_single_case_panels(case):
    """ Check whether 'case > panels' contains the same panels as are
    covered by all the case's interpretations

    args:
        case: dict representing one case in an OpenCGA study

    returns:
        panels_match [bool]: True/False for whether "all panels in
            case>panels" and "all panels contained within the case's
            interpretations" are equivalent
    """

    # Get list of panels in 'case > panels'

    case_panels = []

    for panel in case['panels']:

        case_panels.append({
            'id' : panel['source']['id'],
            'name' : panel['name'],
            'version' : panel['source']['version']})

    # Get list of all panels wihtin all case interpretations

    interpretations = get_all_case_intpns(case)

    all_intpn_panels = []

    for intpn in interpretations:
        for panel in intpn['panels']:

            all_intpn_panels.append({
                'id' : panel['source']['id'],
                'name' : panel['name'],
                'version' : panel['source']['version']})

    # Check whether each intpn panel is in the case panels list

    intpn_panel_count = len(all_intpn_panels)
    matches = 0

    for panel in all_intpn_panels:
        if panel in case_panels:
            matches += 1

    if matches == intpn_panel_count:
        return True

    else:
        return False


def print_case_panels(case):
    """ Print the id, name and version of each panel in a case

    args:
        case: dict representing one case in an OpenCGA study
    """

    for panel in case['panels']:
        print('PanelApp panel {}: {} (version: {})'.format(
            panel['source']['id'],
            panel['name'],
            panel['source']['version']))


def write_all_case_variants(case):
    """ Write out a text file with details of all variants in a case
    (all primaryFindings across all interpretations)

    args:
        case: dict representing one case in an OpenCGA study
    """

    # initialise file with a header

    filename = 'case_{}_all_variant_info.txt'.format(case['id'])

    with open(filename, 'w') as writer:

        writer.write('Interpretation and variant info for case {}\n'.format(
            case['id']))

    # add info from each interpretation

    interpretations = get_all_case_intpns(case)

    for intpn in interpretations:

        if intpn == interpretations[0]:  # primary intpn is 1st element
            with open(filename, 'a') as writer:

                writer.write('\nInterpretation {} (primary)\n'.format(
                    intpn['id']))

        else:  # any subsequent elements are secondary intpns
            with open(filename, 'a') as writer:

                writer.write('\nInterpretation {} (secondary)\n'.format(
                    intpn['id']))

        panels = get_single_intpn_panels(intpn)
        panel_genes = get_single_intpn_panel_genes(intpn)
        var_genes = get_single_intpn_variant_genes(intpn)

        primary_vars = intpn['primaryFindings']

        with open(filename, 'a') as writer:

            # Number of panels in interpretation
            writer.write('\n{} applied {} panels:\n'.format(
                intpn['id'],
                len(panels)))

            # Details of each panel
            for panel in panels:
                writer.write('\tPanel {}: {}, version {}\n'.format(
                    panel['id'],
                    panel['name'],
                    panel['version']))

            # Total uniques genes across all panels
            writer.write('\nPanels covered {} unique genes: {}\n'.format(
                str(len(panel_genes)),
                str(panel_genes)))

            # Total unique genes in which primary findings occurred
            writer.write('\nVariants in {} unique genes: {}\n'.format(
                str(len(var_genes)),
                str(var_genes)))

            # Overlap between panel genes and primary findings genes
            strings = compare_panel_and_variant_genes(intpn)

            for string in strings:
                writer.write('\n{}'.format(string))

            # Selected details of each variant
            writer.write('\n\n{} primary findings:\n'.format(
                len(primary_vars)))

            for variant in primary_vars:

                genes, tiers = get_variant_dict(variant)

                writer.write('{} ({}): {}\n'.format(
                    variant['id'],
                    str(genes),
                    str(tiers)))


""" Interpretation-level functions """


def get_single_intpn_panels(intpn):
    """ Get a list of {id, name, version} for each panel in a case
    interpretation

    args:
        intpn: dict representing one intepretation in an OpenCGA case

    returns:
        output_list: list of dicts, each dict is {id, name, version} for
            one panel of the interpretation
    """

    output_list = []

    for panel in intpn['panels']:

        output_list.append({
            'id' : panel['source']['id'],
            'name' : panel['name'],
            'version' : panel['source']['version']
            })

    return output_list


def get_single_intpn_panel_genes(intpn):
    """ Get a list of the unique genes covered by all the panels in an
    interpretation. Genes are green genes obtained by querying PanelApp
    with each panel's id and version.

    args:
        intpn: dict representing one intepretation in an OpenCGA case

    returns:
        all_genes: list of green genes from the intepretation's panels
    """

    all_genes = []

    for cga_panel in intpn['panels']:

        pa_genes = get_green_genes(
            cga_panel['source']['id'],
            cga_panel['source']['version'])

        for gene in pa_genes:
            if gene not in all_genes:
                all_genes.append(gene)

    all_genes.sort()

    return all_genes


def get_single_intpn_variant_genes(intpn):
    """ Get a list of all genes in an interpretation where at least one
    variant was identified

    args:
        intpn: dict representing one intepretation in an OpenCGA case

    returns:
        all_genes: list of gene symbols
    """

    all_genes = []

    primary_genes = intpn['stats']['primaryFindings']['geneCount']
    secondary_genes = intpn['stats']['primaryFindings']['geneCount']

    for lst in [primary_genes, secondary_genes]:
        for gene in lst.keys():
            if gene not in all_genes:
                all_genes.append(gene)

    all_genes.sort()

    return all_genes


def compare_panel_and_variant_genes(intpn):
    """ Compare the list of unique genes covered by an interpretation's
    panels with the list of genes in which variants were identified

    args:
        intpn: dict representing one intepretation in an OpenCGA case

    returns:
        strings: list of strings describing overlap between genes in
            interpretation panels, and genes in which variants were
            identified
    """

    variant_genes = get_single_intpn_variant_genes(intpn)
    panel_genes = get_single_intpn_panel_genes(intpn)

    strings = []

    if panel_genes in variant_genes:
        strings.append('All panel genes had at least 1 variant identified')

    else:
        strings.append('Some panel genes had no variants identified')

    if (variant_genes == panel_genes) or (variant_genes in panel_genes):
        strings.append('Variants were only identified within panel genes')

    else:
        strings.append('Variants were identified within non-panel genes')

    return strings


def print_intpn_variant_info(intpn):
    """ Get number primary findings in an interpretation, and number of
    variants assigned to each tier

    args:
        intpn: dict representing one intepretation in an OpenCGA case
    """

    primary_tiers = intpn['stats']['primaryFindings']['tierCount']

    primary_tier_list = [
        '{} with tier {}'.format(
            value,
            key
            ) for key, value in primary_tiers.items()]

    print('Interpretation {} ({}) has {} primary findings: {}'.format(
            intpn['id'],
            intpn['description'],
            len(intpn['primaryFindings']),
            primary_tier_list,
        ))


""" Variant-level functions """


def get_variant_dict(variant):
    """ For a single variant, in a single interpretation of an OpenCGA
    case, get a dict of the variant's id, affected genes, and tiers.

    Each variant has as many 'evidences' as its parent interpretation
    has panels. Each evidence is associated with an affected gene and
    the variant's assigned tier. Therefore, a single variant can have
    multiple affected genes and assigned tiers. A variants' tiers -may-
    all be the same, but not necessarily.

    args:
        variant: dict representing a single variant from a single
            interpretation in an OpenCGA case

    returns:
        var_genes: list of unique genes affected by that variant
        var_tiers: list of unique tiers across all variant evidences
    """

    var_genes = []
    var_tiers = []

    for evidence in variant['evidences']:
        if evidence['genomicFeature']['geneName'] not in var_genes:
            var_genes.append(evidence['genomicFeature']['geneName'])

        if evidence['classification']['tier'] not in var_tiers:
            var_tiers.append(evidence['classification']['tier'])

    return var_genes, var_tiers


""" Oddly specific functions """


def print_awkward_variant(oc):
    """ Print out the data for the awkward variant which has two
    different tiers """

    case_info = oc.clinical.info(
        study = 'emee-glh@reanalysis:rd37',
        clinical_analysis = 'SAP-3904-1'
        ).get_results()[0]

    variants = case_info['interpretation']['primaryFindings']

    awkward_variant = '2:47394827:CGTCTCA:C'

    for variant in variants:
        if variant['id'] == awkward_variant:
            for key, value in variant.items():
                print('{}: {}'.format(key, value))


""" Functions on VCFs """


def chromosome_notation(chr_file, vcf_file):
    """
    Convert RefSeq contig notation into UCSC format (needed for filtering)

    Args:
        chr_file [string]: contains chromosome notation info
        sam_file [string]: name of input .sam file

    Output: variants.vcf (approx. 53.6 MB), ~3s
    """

    start_time = dt.now()
    print('{} Fixing chromosome notation'.format(start_time.strftime('%H:%M')))

    # Generate a list of contig name conversions

    chr_names = []

    with open(chr_file, 'r') as reader:
        lines = reader.readlines()

    for line in lines:
        elements = line.split('\t')

        values = [
            elements[0].strip(),  # RefSeq contig name
            elements[1].strip(),  # UCSC contig name
            elements[2].strip(),  # Alternative sequence name
            ]

        chr_names.append(values)

    # Use sed to find and replace strings in the target file

    for names in chr_names:

        # Some contigs don't have a UCSC name, so use the alternative one
        if names[1] == 'na':

            subprocess.run([
                'sed',
                '-i',  # change file in-place
                's/{}/{}/g'.format(names[0], names[2]),
                vcf_file,
                ])

        # Otherwise use the UCSC name
        else:

            subprocess.run([
                'sed',
                '-i',
                's/{}/{}/g'.format(names[0], names[1]),
                vcf_file,
                ])


    end_time = dt.now()
    duration = end_time - start_time
    print('Chromosome notation took {}\n'.format(str(duration)))

    return vcf_file


def filter_vcf(vcf_file, bed_file, minQ, output_prefix):
    """  Filter variants in a .vcf file using QUAL values

    Args:
        vcf_file [string]: .vcf file of variants to be filtered
        bed_file [string]: name of bed file to filter regions on
        minQ [string]: minimum value of QUAL for variants to pass
        output_prefix [string]: prefix for output files

    Command line:   vcftools \
                        --vcf vcf_file/variants.vcf \
                        --bed cardiac_bed.bed \
                        --minQ 20 \
                        --recode \
                        --out vcf_file/filtered

    Outputs:      filtered.vcf (approx. 41.2kB), ~0.2s
    """

    start_time = dt.now()
    print('{} Filtering variants'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'vcftools',
        '--vcf',
        vcf_file,
        '--bed',
        bed_file,
        '--minQ',
        minQ,
        '--recode',
        '--out',
        output_prefix,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Variant filtering took {}\n'.format(str(duration)))

    return '{}.recode.vcf'.format(output_prefix)


def sort_variants(vcf_file, output_file):
    """
    Use 'bcftools sort' to sort a .vcf file

    Args:
        vcf_file [string]: name of input .bam file
        output_file [string]: name for output .bam file

    Command line:   bcftools sort \
                        -o vcf_file/sorted.vcf \
                        vcf_file/filtered.recode.vcf

    Outputs:        sorted.vcf (approx. 41.2 kB), <1s
    """

    start_time = dt.now()
    print('{} Sorting .vcf file'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bcftools',
        'sort',
        '-o',
        output_file,
        vcf_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Sorting took {}\n'.format(str(duration)))

    return output_file


""" Main """


def main():

    """ Session setup """

    user = 'jmiles'  # terminal will prompt for password
    host = 'https://uat.eglh.app.zettagenomics.com/opencga'


    """ Case details """

    project = 'reanalysis'
    study = 'rd37'
    study_id = 'emee-glh@{}:{}'.format(project, study)
    # case_id = 'SAP-48034-1'  # case ID for the VCF
    case_id = 'CSA-1001-1'  # small case


    """ Required/useful files """

    vcf_file = 'data/vcfs/SAP-48034-1.vcf'  # VCF file from selected case
    chr_file = 'data/refs/convert_scaffolds'  # refseq/UCSC/alt contig aliases
    hgnc_dump = 'data/refs/hgnc_dump_20210727.txt'  # HGNC website data dump
    current_panels = 'data/refs/current_pa_panel_ids.txt'  # current PA panels


    """ Function implementations """

    # ## LAST UPDATED: Tues 24 May 22
    # current_panels = write_current_pa_list()

    ## Set up an OpenCGA session
    oc = session_setup(user, host)

    # ## Retrieve all data from a specific case
    case = get_single_case_data(oc, study_id, case_id)

    # ## Get info for all panels in a case including original/current genes
    # hgnc_df = pd.read_csv(hgnc_dump, sep='\t')
    # panel_list = get_case_panels(case, hgnc_df)

    # ## Construct bed files for original and current gene lists
    # original_bed, current_bed = get_biomart_output(case_id, panel_list)

    # ## Filter VCF on bed file and QUAL filter
    # fixed_vcf = chromosome_notation(chr_file, vcf_file)
    # filtered_vcf = filter_vcf(fixed_vcf, current_bed, '20', 'data/vcfs/filter')
    # sorted_vcf = sort_variants(filtered_vcf, 'data/vcfs/sorted.vcf')

    ## Get a list of phenotype IDs for the case
    phenotypes = get_case_phenotypes(case)
    print(phenotypes)


if __name__ == '__main__':
    main()
