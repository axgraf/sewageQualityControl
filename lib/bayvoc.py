# Created by alex at 06.10.23

from datetime import datetime
from requests.auth import HTTPBasicAuth
import requests
import json
from base64 import b64encode
import jwt
import pandas as pd
from .config import *


class BayVOC:

    def __init__(self, config: Config):
        self.config = config
        self.__authenticate()
        self.column_map = {
    'collectionDate': 'collectionDate',
    'bem_lab': 'bem_lab',
    'additionalInfo': 'bem_pn',
    'labCopiesMlGeneN1': 'biomarker_N1',
    'labCopiesMlGeneN2': 'biomarker_N2',
    'labCopiesMlGeneN3': 'biomarker_N3',
    'labCopiesMlGeneE': 'biomarker_E',
    'labCopiesMlGeneORF': 'biomarker_ORF',
    'labCopiesMlGeneRDRP': 'biomarker_RDRP',
    'nh4n': 'nh4n',
    'lf': "lf",
    'mean_sewage_flow': 'mean_sewage_flow',
    'crassphage': 'crassphage',
    'pmmov': 'pmmov',
    'trockentag': 'trockentag'
}



    def __authenticate(self):
        try:
            response = requests.post(self.config.authentication_endpoint,
                                     json={
                                         "username": self.config.bay_voc_username,
                                         "password": self.config.bay_voc_password
                                     }, auth=HTTPBasicAuth(self.config.basic_auth_username,
                                                           self.config.basic_auth_password))
            if response.status_code == 200:
                print("new token requested")
                auth_data = json.loads(response.text)
                self.token = auth_data['jwttoken']
                decoded_token = jwt.decode(self.token, options={"verify_signature": False})
                self.expiration_time = datetime.fromtimestamp(decoded_token['exp'])
            else:
                logging.error("Authentication error")
                logging.error(response.text)
                exit(1)
        except requests.exceptions.ConnectionError as e:
            print(e)
            logging.error("Can not authenticate to Bay-VOC Rest API - Is the service down?")

    def __create_authentication_header(self):
        if self.token:
            basic_auth = b64encode(
                bytes(f"{self.config.basic_auth_username}:{self.config.basic_auth_password}", "utf-8")).decode(
                "ascii")
            basic_auth_str = "Basic {}".format(basic_auth)
            bearer_token_str = "Bearer {}".format(self.token)
            return {'Authorization': basic_auth_str,
                    'AuthorizationBearer': bearer_token_str
                    }
        else:
            logging.error("No JWT token obtained from Bay-VOC hub")

    def read_all_sewagesamples_from_db(self):
        response = requests.get(self.config.sewage_get_all_samples_endpoint,
                                headers=self.__create_authentication_header())
        json_sewage_samples = response.json()
        df = pd.DataFrame(json_sewage_samples)
        df.rename(columns=self.column_map, inplace=True)
        df['bem_lab'] = ''
        sample_locs = df.groupby('name')
        sewage_samples_dict = dict()
        for name, group_list in sample_locs:
            group_list = group_list[list(self.column_map.values())]
            if 'A-STADT' in name :
                sewage_samples_dict[name] = group_list

        return sewage_samples_dict
