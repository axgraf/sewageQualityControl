# Created by alex at 28.06.23
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def plot_biomarker_outlier_summary(measurements_df, sample_location, output_folder, outlier_detection_methods):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    is_outliers_calculated = measurements_df.loc[:, measurements_df.columns.str.contains("outlier")].empty
    biomarker_ratio_columns = measurements_df.filter(regex="^biomarker_[a-zA-Z0-9]+/biomarker_").dropna(axis=1, how='all')
    sns.set_style('darkgrid')
    sns.set_palette(sns.color_palette(['black', '#1f77b4', '#ff7f0e']))
    plot_frame = pd.DataFrame()

    for idx, biomarker_ratio in enumerate(biomarker_ratio_columns.columns):
        dat = pd.DataFrame()
        outlier_method = "outlier_" + biomarker_ratio
        dat['date'] = measurements_df['collectionDate']
        dat['ratio'] = measurements_df[biomarker_ratio]
        dat['biomarker_ratio'] = biomarker_ratio
        dat['outlier'] = measurements_df[outlier_method] if not is_outliers_calculated else np.nan
        plot_frame = plot_frame.append(dat)
    plot_frame['outlier'] = np.where((plot_frame['outlier'].isna()) , 'not tested', np.where((plot_frame['outlier'] == False), "inlier",
                                     "outlier"))

    g = sns.FacetGrid(plot_frame, col="biomarker_ratio", col_wrap=1, margin_titles=True, height=2.5, aspect=6, sharey=False)
    g.map_dataframe(sns.scatterplot, x="date", y="ratio", hue="outlier")
    g.add_legend()
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Biomarker ratios for '{}' -  Outlier detection methods: {}".format(sample_location, outlier_detection_methods), fontsize=16)
    plt.show()
    g.savefig(os.path.join(output_folder, "Biomarker_qc_{}.png".format(sample_location)), dpi=300)





