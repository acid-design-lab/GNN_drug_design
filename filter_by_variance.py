# -*- coding: utf-8 -*-
"""filter_by_variance

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1S8juEkdJbdUU8jfuKUxpb-gONUcUAB5X
"""

def filter_by_variance(dataframe,
                        var_threshold,
                        features_to_keep,
                        target_feature):
        """
        :param dataframe: table in the format of DataFrame object.
        :param var_threshold: the minimal level of tolerable feature variance.
        :param features_to_keep: list of column names of the features to keep.
        :param target_feature: column name of a target feature to predict further.
        :return: a DataFrame object M columns (where M <= N) with chemical descriptors from Rdkit
        + 1 column with disease categories.
        """

        features_to_keep.append(target_feature)

        descriptors_to_filt = dataframe.drop(['Disease'] + features_to_keep, axis = 1)
        var_filtered_dataset_with_descriptors = pd.DataFrame()

        to_save = dataframe[features_to_keep]

        scaler = preprocessing.MinMaxScaler()
        scaled_descriptors = pd.DataFrame(scaler.fit_transform(descriptors_to_filt), columns = descriptors_to_filt.columns)

        filter = scaled_descriptors.loc[:, scaled_descriptors.var() < var_threshold]
        filtered_data = descriptors_to_filt.drop(columns = filter.columns , axis = 1)


        var_filtered_dataset_with_descriptors = filtered_data.join(to_save, how = 'outer')


        return var_filtered_dataset_with_descriptors