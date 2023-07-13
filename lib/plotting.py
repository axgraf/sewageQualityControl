# Created by alex at 28.06.23
import os
import itertools
import dateutil
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .constant import *


sns.set_style('darkgrid')
sns.set_context("talk")


def get_label_colors():
    colors = {'not tested': 'dimgrey',
              'inlier': 'royalblue',
              'outlier': 'orangered'
              }
    return colors


def create_output_folder(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def plot_biomarker_outlier_summary(measurements_df, sample_location, output_folder, outlier_detection_methods, interactive=False):
    create_output_folder(output_folder)
    plot_frame = pd.DataFrame()
    for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
        biomarker_ratio = biomarker1 + "/" + biomarker2
        length = len(measurements_df[biomarker_ratio].dropna())
        if length > 0:
            dat = pd.DataFrame()
            dat['date'] = measurements_df[Columns.DATE.value]
            dat['log2[ratio]'] = measurements_df[biomarker_ratio]
            dat['biomarker_ratio'] = biomarker_ratio
            dat['outlier'] = measurements_df[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)]
            plot_frame = pd.concat([plot_frame, dat])
    if plot_frame.shape[0] > 0:
        plot_frame['outlier'] = np.where(
            (SewageFlag.is_flag_set_for_series(plot_frame['outlier'], SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES)),
            'not tested',
            np.where((SewageFlag.is_flag_set_for_series(plot_frame['outlier'], SewageFlag.BIOMARKER_RATIO_OUTLIER)), 'outlier', 'inlier'))

        g = sns.FacetGrid(plot_frame, col="biomarker_ratio", col_wrap=1, margin_titles=True, height=2.5, aspect=6, sharey=False, legend_out=True)
        min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
        max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
        g.set(xlim=(min_date, max_date))
        g.map_dataframe(sns.scatterplot, x="date", y="log2[ratio]", hue="outlier", palette=get_label_colors())
        g.add_legend(frameon=True)
        g.set_titles(row_template='{row_name}', col_template='{col_name}')
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle("Biomarker ratios for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods))
        if interactive:
            plt.show()
        g.savefig(os.path.join(output_folder, "Biomarker_qc_{}.png".format(sample_location)), dpi=300)
        plt.cla()
        plt.close()
        return g.fig


def plot_surrogatvirus (measurements_df, sample_location, output_folder, outlier_detection_methods, interactive=False):
    create_output_folder(output_folder)

    plot_frame = pd.DataFrame()

    for sVirus in Columns.get_surrogatevirus_columns():
        dat = pd.DataFrame()
        dat['date'] = measurements_df[Columns.DATE.value]
        dat['value'] = measurements_df[sVirus]
        dat['type'] = sVirus
        dat['outlier'] = np.where(
            (SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.get_surrogate_flag(sVirus)],
                                               SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE)),
            'not tested', np.where((SewageFlag.is_flag_set_for_series(
                measurements_df[CalculatedColumns.get_surrogate_outlier_flag(sVirus)], SewageFlag.SURROGATEVIRUS_OUTLIER)),
                'outlier', 'inlier'))
        plot_frame = pd.concat([plot_frame, dat])

    g = sns.FacetGrid(plot_frame, col="type", col_wrap=1, margin_titles=True, height=3.5, aspect=6, sharey=False, legend_out=True)
    min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
    max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
    g.set(xlim=(min_date, max_date))
    g.map_dataframe(sns.scatterplot, x="date", y="value", hue="outlier", palette=get_label_colors())
    g.add_legend(frameon=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Surrogatvirus quality control for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods))
    if interactive:
        plt.show()
    g.savefig(os.path.join(output_folder, "Surrogatvirus_quality_control_{}.png".format(sample_location)), dpi=300)
    plt.cla()
    plt.close()
    return g.fig


def plot_water_quality(measurements_df, sample_location, output_folder, outlier_detection_methods, interactive=False):
    create_output_folder(output_folder)

    plot_frame = pd.DataFrame()
    for qual_type, outlier_flag, not_enough_flag in zip([Columns.AMMONIUM.value, Columns.CONDUCTIVITY.value],
                                                        [SewageFlag.AMMONIUM_OUTLIER.value, SewageFlag.CONDUCTIVITY_OUTLIER.value],
                                                        [SewageFlag.NOT_ENOUGH_AMMONIUM_VALUES.value, SewageFlag.NOT_ENOUGH_CONDUCTIVITY_VALUES.value]):
        dat = pd.DataFrame()
        dat['date'] = measurements_df[Columns.DATE.value]
        dat['value'] = measurements_df[qual_type]
        dat['type'] = qual_type
        dat['outlier'] = measurements_df[CalculatedColumns.FLAG.value]
        dat['outlier'] = np.where((measurements_df[CalculatedColumns.FLAG.value] & not_enough_flag > 0), 'not tested',
                                     np.where((measurements_df[CalculatedColumns.FLAG.value] & outlier_flag > 0), 'outlier', 'inlier'))
        plot_frame = pd.concat([plot_frame, dat])

    g = sns.FacetGrid(plot_frame, col="type", col_wrap=1, margin_titles=True, height=3.5, aspect=6, sharey=False, legend_out=True)
    min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
    max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
    g.set(xlim=(min_date, max_date))
    g.map_dataframe(sns.scatterplot, x="date", y="value", hue="outlier", palette=get_label_colors())
    g.add_legend(frameon=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Water quality control for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods))
    if interactive:
        plt.show()
    g.savefig(os.path.join(output_folder, "Water_quality_control_{}.png".format(sample_location)), dpi=300)
    plt.cla()
    plt.close()
    return g.fig


def plot_sewage_flow(measurements_df, sample_location, output_folder, is_mean_dry_weather_flow_available, dry_weather_flow, interactive=False):
    create_output_folder(output_folder)
    plot_frame = pd.DataFrame()
    plot_frame['date'] = measurements_df[Columns.DATE.value]
    plot_frame['value'] = measurements_df[Columns.MEAN_SEWAGE_FLOW.value]
    plot_frame['outlier'] = np.where(
        (SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value], SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION)) |
        (SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value], SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO)),
        'outlier', np.where((SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value], SewageFlag.SEWAGE_FLOW_NOT_ENOUGH_PREVIOUS_VALUES)),
        'not tested', 'inlier'))
    plt.figure(figsize=(20, 8))
    g = sns.scatterplot(data=plot_frame, x="date", y="value", hue="outlier", palette=get_label_colors())
    if dry_weather_flow:
        g.axhline(dry_weather_flow, label="Dry weather flow", linestyle='dashed', c='black')
    plt.legend(bbox_to_anchor=(1.01, 0.5), loc='center left', borderaxespad=0)
    if dry_weather_flow:
        dry_weather_flow_type = 'estimated' if is_mean_dry_weather_flow_available else 'specified'
        plt.title("Mean sewage flow for '{}' - Dry weather flow: '{}' - Type: '{}'".format(sample_location, round(dry_weather_flow,1), dry_weather_flow_type))
    else:
        plt.title("Mean sewage flow for '{}'".format(sample_location))
    plt.ylabel('Mean sewage flow')
    plt.tight_layout()
    if interactive:
        plt.show()
    g.get_figure().savefig(os.path.join(output_folder, "Sewage_flow_quality_control_{}.png".format(sample_location)), dpi=300)
    plt.cla()
    plt.close()
    return g.figure


def plot_biomarker_normalization(measurements_df, sample_location, output_folder, interactive=False):
    create_output_folder(output_folder)
    plot_frame = pd.DataFrame()
    for column_type in [CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value, CalculatedColumns.BASE_REPRODUCTION_FACTOR.value]:
        dat = pd.DataFrame()
        dat['date'] = measurements_df[Columns.DATE.value]
        dat['value'] = measurements_df[column_type]
        dat['type'] = column_type
        dat['outlier'] = np.where(
            (SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value], SewageFlag.REPRODUCTION_NUMBER_OUTLIER)),
            'outlier', np.where((SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value], SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)),
                                'not tested', 'inlier'))
        plot_frame = pd.concat([plot_frame, dat])
    g = sns.FacetGrid(plot_frame, col="type", col_wrap=1, margin_titles=True, height=3.5, aspect=6, sharey=False, legend_out=True)
    min_date = plot_frame['date'].min() + dateutil.relativedelta.relativedelta(days=-2)
    max_date = plot_frame['date'].max() + dateutil.relativedelta.relativedelta(days=2)
    g.set(xlim=(min_date, max_date))
    g.map_dataframe(sns.scatterplot, x="date", y="value", hue="outlier", palette=get_label_colors())
    g.add_legend(frameon=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Biomarker normalization for '{}'".format(sample_location))
    g.fig.axes[1].set_yscale("symlog", base=2)
    if interactive:
        plt.show()
    g.savefig(os.path.join(output_folder, "Biomarker_normalization_{}.png".format(sample_location)), dpi=300)
    plt.cla()
    plt.close()
    return g.fig


def plot_general_outliers(measurements_df, sample_location, output_folder, interactive=False):
    create_output_folder(output_folder)
    plot_frame = pd.DataFrame()
    plot_frame['date'] = measurements_df[Columns.DATE.value]
    plot_frame['value'] = measurements_df[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value]
    plot_frame['outlier'] = measurements_df[CalculatedColumns.OUTLIER_REASON.value]
    plt.figure(figsize=(20, 10))
    g = sns.scatterplot(data=plot_frame, x="date", y="value", hue="outlier")
    sns.move_legend(
        g, loc="lower center",
        bbox_to_anchor=(.5, 1.2), ncol=2, title=None, frameon=True
    )
    plt.title("Outliers - Normalized mean biomarkers for '{}'".format(sample_location), pad=10)
    plt.subplots_adjust(top=0.72)
    if interactive:
        plt.show()
    g.get_figure().savefig(os.path.join(output_folder, "Outlier_{}.png".format(sample_location)), dpi=300)
    plt.cla()
    plt.close()
    return g.figure


def plot_all(plots_arr: [], output_file: str):
    with PdfPages('count.pdf') as pdf_pages:
        for i, fig in enumerate(plots_arr):
            figure = plt.figure(fig)
            pdf_pages.savefig(figure)
        pdf_pages.close()

