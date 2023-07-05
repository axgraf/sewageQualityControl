# Created by alex at 28.06.23
import os
import itertools
import dateutil
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from .constant import *

sns.set_style('darkgrid')


def get_label_colors():
    colors = {'not tested': 'dimgrey',
              'inlier': 'royalblue',
              'outlier': 'orangered'
              }
    return colors


def create_output_folder(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def plot_biomarker_outlier_summary(measurements_df, sample_location, output_folder, outlier_detection_methods):
    create_output_folder(output_folder)

    plot_frame = pd.DataFrame()
    for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
        biomarker_ratio = biomarker1 + "/" + biomarker2
        length = len(measurements_df[biomarker_ratio].dropna())
        if length > 0:
            dat = pd.DataFrame()
            dat['date'] = measurements_df[Columns.DATE.value]
            dat['ratio'] = measurements_df[biomarker_ratio]
            dat['biomarker_ratio'] = biomarker_ratio
            dat['outlier'] = measurements_df[Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)]
            plot_frame = plot_frame.append(dat)
    if plot_frame.shape[0] > 0:
        plot_frame['outlier'] = np.where(
            (SewageFlag.is_flag_set_for_series(plot_frame['outlier'], SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES)),
            'not tested',
            np.where((SewageFlag.is_flag_set_for_series(plot_frame['outlier'], SewageFlag.BIOMARKER_RATIO_OUTLIER)), 'outlier', 'inlier'))

            #(plot_frame['outlier'] & SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES.value > 0), 'not tested',
            #np.where((plot_frame['outlier'] & SewageFlag.BIOMARKER_RATIO_OUTLIER.value > 0), 'outlier', 'inlier'))
        #plot_frame['outlier'] = np.where((plot_frame['outlier'] & SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES.value > 0), 'not tested',
        #                                 np.where((plot_frame['outlier'] & SewageFlag.BIOMARKER_RATIO_OUTLIER.value > 0), 'outlier', 'inlier' ))

        g = sns.FacetGrid(plot_frame, col="biomarker_ratio", col_wrap=1, margin_titles=True, height=2.5, aspect=6, sharey=False)
        min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
        max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
        g.set(xlim=(min_date, max_date))
        g.map_dataframe(sns.scatterplot, x="date", y="ratio", hue="outlier", palette=get_label_colors())
        #g.add_legend(legend_data=get_label_colors())
        g.add_legend(legend_data={
            key: value for key, value in zip(get_label_colors().keys(), g._legend_data.values())
        })
        g.set_titles(row_template='{row_name}', col_template='{col_name}')
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle("Biomarker ratios for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods), fontsize=16)
        plt.show()
        g.savefig(os.path.join(output_folder, "Biomarker_qc_{}.png".format(sample_location)), dpi=300)


def plot_water_quality(measurements_df, sample_location, output_folder, outlier_detection_methods):
    create_output_folder(output_folder)

    plot_frame = pd.DataFrame()
    for qual_type, outlier_flag, not_enough_flag in zip([Columns.AMMONIUM.value, Columns.CONDUCTIVITY.value],
                                                        [SewageFlag.AMMONIUM_OUTLIER.value, SewageFlag.CONDUCTIVITY_OUTLIER.value],
                                                        [SewageFlag.NOT_ENOUGH_AMMONIUM_VALUES.value, SewageFlag.NOT_ENOUGH_CONDUCTIVITY_VALUES.value]):
        dat = pd.DataFrame()
        dat['date'] = measurements_df[Columns.DATE.value]
        dat['value'] = measurements_df[qual_type]
        dat['type'] = qual_type
        dat['outlier'] = measurements_df[Columns.FLAG.value]
        dat['outlier'] = np.where((measurements_df[Columns.FLAG.value] & not_enough_flag > 0), 'not tested',
                                     np.where((measurements_df[Columns.FLAG.value] & outlier_flag > 0), 'outlier', 'inlier'))
        plot_frame = plot_frame.append(dat)

    g = sns.FacetGrid(plot_frame, col="type", col_wrap=1, margin_titles=True, height=3.5, aspect=6, sharey=False)
    min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
    max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
    g.set(xlim=(min_date, max_date))
    g.map_dataframe(sns.scatterplot, x="date", y="value", hue="outlier", palette=get_label_colors())
    g.add_legend()
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Water quality control for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods),  fontsize=16)
    plt.show()
    g.savefig(os.path.join(output_folder, "Water_quality_control_{}.png".format(sample_location)), dpi=300)


