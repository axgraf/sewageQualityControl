# Created by alex at 15.06.23
#import sys
#import logging
#from arcgis.gis import GIS
#from arcgis.features import FeatureLayer
#from .config import Config
#from .sewage import SewageSample



class Arcgis:
    pass

#     def __init__(self, config: Config):
#         self.sewage_plants = dict()
#         self.regions2plants = dict()
#         self.sewage_samples = dict()
#         self.config = config
#         self.__connnect_gis(config.gis_url, config.gis_user, config.gis_password)
#
#     def obtain_sewage_samples(self):
#         self.__connnect_gis(self.config.gis_url, self.config.gis_user, self.config.gis_password)
#         self.__get_sewage_plants()
#         self.__get_messwerte()
#         import pickle
#         import os
#         if not os.path.exists('sewageData.dat'):
#             with open('sewageData.dat', 'wb') as f:
#                 pickle.dump(self.sewage_samples, f)
#         if not os.path.exists('sewagePlantData.dat'):
#             with open('sewagePlantData.dat', 'wb') as f:
#                 pickle.dump(self.sewage_plants, f)
#         return self.sewage_samples, self.sewage_plants
#
#     #        self.__loess_regression()
#
#     def __connnect_gis(self, gis_url, user, password):
#         logging.info("Connect to ARCGIS server...")
#         self.gis = GIS(gis_url, user, password)
# #        groups = self.gis.groups.search('title: LB Bayern (LBBY)')  # obtain group
# #        group_content = groups[0].content()
#
#         self.monitoring_daten = self.gis.content.get('d3b1c622cceb40e48353da110fde73b8')
#         self.messstellen_bayern = self.gis.content.get('2de71e70d30945768f24b71158a38d66')
#
#
#     def __get_sewage_plants(self):
#         logging.info("Obtain sewage plants infos...")
#         self.sewage_plants = dict()
#         sewage_plant_layer = FeatureLayer(
#             'https://services-eu1.arcgis.com/e0dlK9aWS0lF3hlT/arcgis/rest/services/PTKA_DB_Mod_v3_VIEW_ReadONLY_d0bb0/FeatureServer/0')  # get content of group, used to get Ids for the data below
#         sewage_plants = sewage_plant_layer.query().features
#         for feature in sewage_plants:
#             self.sewage_plants[feature.attributes['NAME']] = feature.attributes['TW_ABFLUSS']
#
#
#     def __get_messwerte(self):
#         logging.info("Obtain measurements...")
#         messstellen_url = self.messstellen_bayern.layers[0].url
#         messtellen_feature = FeatureLayer(messstellen_url)
#         names = set()
#         sample_features = messtellen_feature.query().features
#         for feature in sample_features:
#             test = feature.as_dict
#             names.add(test['attributes']['NAME'])
#             sewage_sample = SewageSample(feature.as_dict)
#             #            self.__map_sewage_location(sewage_sample)
#             # replace short location name with long name from sewage plant
#             if 'M_01' in sewage_sample.location_name:
#                 print("here")
#             if not sewage_sample.location_name == "A-LK_02_STADTBERGEN":  # Mix aus Königsbrunn + 20% Stadtbergen und dient als Referenzsstandort (nicht mit aufnehmen)
#                 if sewage_sample.has_collection_date():
#                     if sewage_sample.location_name in self.sewage_samples:
#                         if any(loc in sewage_sample.location_name for loc in
#                                ("BGL_01", "BGL_02", "BGL_03", "BGL_04", "BGL_05")):  # BGL only after March 2023
#                             if sewage_sample.collectionDate >= '2023-03-01':
#                                 self.sewage_samples[sewage_sample.location_name].append(sewage_sample)
#                         else:
#                             self.sewage_samples[sewage_sample.location_name].append(sewage_sample)
#                     else:
#                         if any(loc in sewage_sample.location_name for loc in
#                                ("BGL_01", "BGL_02", "BGL_04", "BGL_05")):  # BGL only after March 2023
#                             if sewage_sample.collectionDate >= '2023-03-01':
#                                 self.sewage_samples[sewage_sample.location_name] = [sewage_sample]
#                         else:
#                             self.sewage_samples[sewage_sample.location_name] = [sewage_sample]
#                 else:
#                     print("no collection date:\t{}".format(sewage_sample.location_name))
#         logging.info("\t\t{} measurements obtained".format(len(self.sewage_samples)))
