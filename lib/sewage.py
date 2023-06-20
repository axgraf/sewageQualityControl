# Created by alex at 15.06.23
import sys
import logging
import itertools
from datetime import datetime
from .flags import Flag

class SewageSample:

    def __init__(self, feature):
        #self.latitude = feature['geometry']['x']
        #self.longitude = feature['geometry']['y']
        #self.spatialReference = feature['geometry']['spatialReference']['latestWkid']
        self.arcgisId = feature['attributes']['ObjectId']
        self.location_name = feature['attributes']['NAME']
        #self.description = feature['attributes']['BESCHREIBUNG']
        #self.locationType = feature['attributes']['TYP']
        timestamp = feature['attributes']['ENDE']
        self.collectionDate = None
        if timestamp:
            self.collectionDate = datetime.fromtimestamp(timestamp / 1000).strftime('%Y-%m-%d')
        else:
            logging.warning("Warning: No collection date. Sample location name: '{}' with ARCGIS Id: '{}'".format(self.location_name, self.arcgisId))
        self.mean_sewage_flow = feature['attributes']['VOLUMENSTROM']
        #self.weather = feature['attributes']['WETTER']
        #self.sampleProcedureType = feature['attributes']['PN_ART']
        #self.sampleType = feature['attributes']['P_ART']
        self.crassphage = feature['attributes']['CRASSPHAGE']
        self.pmmov = feature['attributes']['PMMOV']
        self.biomarker_N1 = feature['attributes']['N1_LAB']
        self.biomarker_N2 = feature['attributes']['N2_LAB']
        self.biomarker_N3 = feature['attributes']['N3_LAB']
        self.biomarker_E = feature['attributes']['E_LAB']
        self.biomarker_ORF = feature['attributes']['ORF_LAB']
        self.biomarker_RDRP = feature['attributes']['RDRP_LAB']
        self.labCopiesMWMeanGenes = feature['attributes']['GEN_MW_LAB']  #  Der Mittelwert aus den im Labor bestimmten Genen, also den Einträgen der Attribute X_LAB z.B. N1_LAB.
        self.labCopiesGMWeanGenes = feature['attributes']['GEN_GM_LAB']  #  Der gleitende Mittelwert aus den Befunden der letzten 3 Probenahmetage
        self.trockentag = feature['attributes']['TRO_TAG']
        self.air_temperature = feature['attributes']['L_TEMP']
        self.water_temperature = feature['attributes']['W_TEMP']
        self.nh4n = feature['attributes']['NH4N']  # Ammonium
        self.pH = feature['attributes']['PH']
        self.lf = feature['attributes']['LF']  # Leitfähigkeit
        self.mean_norm_genes = feature['attributes']['GEN_MW_NORM']  # Der Mittelwert aus den normierten Zielgenen von X-normiert z.B. N1_NORM.
        self.rolling_average_3_measurements = feature['attributes']['GEN_GMW_NORM']  # Der gleitende Mittelwert aus dem normierten Mittelwert der Befunde (Zielgene-Mw-Normiert GEN_MW_NORM) der letzten 3 Probenahmetage
        self.rolling_average_5_measurements = feature['attributes']['Q_MW']
        self.bem_qs = feature['attributes']['BEM_QS']
        self.bem_pn = feature['attributes']['BEM_PN']
        self.bem_lab = feature['attributes']['BEM_LAB']
        self.qualified = feature['attributes']['QUALI']
        self.location_id = feature['attributes']['GlobalID']
        self.rain_region_id = feature['attributes']['ID_REGEN_EZG']
        self.loess_regression_value_last_release = None
        self.loess_percentage_difference_seven_days_before_last_release = None
        self.loess_regression_value = None
        self.loess_percentage_difference_seven_days_before = None

        self.flags = []
        self.discarded = False

    def has_collection_date(self):
        return bool(self.collectionDate)

    def has_neccessary_values(self):
        return self.has_collection_date and self.mean_norm_genes and self.qualified != "Nein"

    def is_biomarker_below_threshold(self, biomarker_name, threshold):
        if biomarker_name in self.get_biomarkers_dict():
            biomarker_value = self.get_biomarkers_dict()[biomarker_name]
            if biomarker_value is None or biomarker_value < threshold:
                return True
        return False

    def get_biomarkers_dict(self):
        biomarkers_dict = {
            'N1': self.biomarker_N1,
            'N2': self.biomarker_N2,
            'N3': self.biomarker_N3,
            'E': self.biomarker_E,
            'ORF': self.biomarker_ORF,
            'RDRP': self.biomarker_RDRP,
        }
        return biomarkers_dict

    def get_biomarker_ratios(self, threshold):
        biomarker_dict = self.get_biomarkers_dict()
        biomarker_ratio = dict()
        res = list(itertools.combinations(biomarker_dict, 2))
        for biomarker_tuple in res:
            biomarker1 = biomarker_tuple[0]
            biomarker2 = biomarker_tuple[1]
            test = biomarker_dict[biomarker1]
            if not self.is_biomarker_below_threshold(biomarker1, threshold) and not self.is_biomarker_below_threshold(biomarker2, threshold):
                ratio = biomarker_dict[biomarker1] / biomarker_dict[biomarker2]
            else:
                ratio = None
            biomarker_ratio[biomarker1+"_"+biomarker2] = ratio
        return biomarker_ratio

    def has_flags(self):
        return len(self.flags) > 0

    def set_flag(self, flag: Flag):
        self.flags.append(flag)

