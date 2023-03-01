#!usr/bin/env python

""" OpenCGA

Project: reanalysis (2427 cases total)
Study rd37 (625 cases)
Study rd38 (1802 cases)


"""

import os
import json
import requests

from datetime import datetime as dt

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpencgaClient


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


def get_project_ids(oc):
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

    for id in project_ids:
        print(f'project: {id}')

    return project_ids


def get_projects_study_ids(oc, project_id):
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

    for id in study_ids:
        print(f'study: {id}')

    return study_ids


def get_studys_case_ids(oc, study_id):
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

    print(f'Study {study_id} contains {len(case_ids)} cases')

    cases_file = f'do_not_upload/data_structure/opencga_{study_id}_cases.txt'

    with open(cases_file, 'w') as writer:
        for id in case_ids:

            writer.write(f'{id}\n')

    return case_ids


def write_all_opencga_cases_list(oc, project_id):

    all_cases = []
    output_fp = 'do_not_upload/data_structure/opencga_all_cases.txt'
    project_studies = oc.studies.search(project=project_id).get_results()
    study_ids = [study['id'] for study in project_studies]

    for study in study_ids:
        cases = oc.clinical.search(study=study, include='id').get_results()
        for case in cases:
            all_cases.append(case['id'])

    all_cases.sort()

    with open(output_fp, 'w') as writer:
        for id in all_cases:
            writer.write(f'{id}\n')


def get_single_case(oc, study_id, case_id):
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


def write_single_case(oc, study_id, case_id):
    """ Retrieve and write out data from a single case.

    args:
        oc: OpenCGA client
        study_id [str]: id specifying project and study
        case_id [str]
    """

    dir = 'do_not_upload/data_structure'
    case_file = f'{dir}/opencga_{study_id[-4:]}_{case_id}.json'
    case_data = get_single_case(oc, study_id, case_id)

    with open(case_file, 'w') as writer:
        json.dump(case_data, writer)


def compare_case_lists(opencga_file):
    """ Compare the list 'cases I have JSONs for' with the list 'cases
    in the OpenCGA reanalysis project'.

    args:
        opencga_file [fp]: lists all cases in all studies in reanalysis project
    """

    with open(opencga_file, 'r') as reader:
        lines = reader.readlines()  # all reanalysis cases in opencga

    opencga_cases = [line.strip() for line in lines if line.strip()]

    cases_100k = []  # all 100k cases that I have jsons for
    cases_cancer = []  # cases which are cancer program
    cases_rd = []  # cases which are rd program
    rd_covered = []  # rd cases which are in opencga reanalysis project
    rd_missing = []  # rd cases which are not

    input = 'do_not_upload/jsons/20230104'
    output = 'do_not_upload/data_structure/ians_jsons.tsv'

    with open(output, 'w') as writer:
        writer.write('100k cases for which JSONs are available:\n')

    for json_file in os.listdir(input):

        json_path = f"{input}/{json_file}"
        with open(json_path, 'r') as reader:
            case_json = json.load(reader)

        case_id = case_json['case_id'].strip()
        program = case_json['program'].strip()

        cases_100k.append(case_id)

        if program == 'cancer':
            cases_cancer.append(case_id)

        elif program == 'rare_disease':
            cases_rd.append(case_id)

            if case_id in opencga_cases:
                rd_covered.append(case_id)

            else:
                rd_missing.append(case_id)

        with open(output, 'a') as writer:
            writer.write(f'{case_id}\t{program}\n')

    print(
        f"{len(cases_100k)} total 100k cases: "
        f"{len(cases_rd)} RD and {len(cases_cancer)} cancer")

    print(
        f"RD cases: "
        f"{len(rd_covered)} in OpenCGA, {len(rd_missing)} not in OpenCGA")

    assert len(cases_rd) + len(cases_cancer) == len(cases_100k)
    assert len(rd_covered) + len(rd_missing) == len(cases_rd)


def main():

    user = 'jmiles'
    host = 'https://uat.eglh.app.zettagenomics.com/opencga'

    project = 'reanalysis'
    studies = ['rd37', 'rd38']
    study = 'rd37'
    study_id = f'emee-glh@{project}:{study}'
    opencga_cases = 'do_not_upload/data_structure/opencga_all_cases.txt'
    case_id = 'CSA-1000-1'

    oc = session_setup(user, host)

    # # get all case ids for opencga rd projects

    # projects = get_project_ids(oc)
    # for project in projects:
    #     print(project)
    #     studies = get_projects_study_ids(oc, project)
    #     if project == 'reanalysis':
    #         for study in studies:
    #             cases = get_studys_case_ids(oc, study)

    # write_all_opencga_cases_list(oc, project)
    # compare_case_lists(opencga_cases)

    write_single_case(oc, study_id, case_id)


if __name__ == '__main__':
    main()
