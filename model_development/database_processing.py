def main():
    def generate_descriptors(file_path,
                             smiles_column_name):
        """

        :param file_path: path to the .csv table.
        :param smiles_column_name: string reflecting the SMILES column name.

        :return: a DataFrame object containing N columns with chemical descriptors from Rdkit
        + 1 column with disease categories.

        """

        dataset_with_descriptors = None

        return dataset_with_descriptors

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

        var_filtered_dataset_with_descriptors = None

        return var_filtered_dataset_with_descriptors

    def drop_correlations(dataframe,
                          corr_threshold,
                          features_to_keep,
                          target_feature):
        """

        :param dataframe: table in the format of DataFrame object.
        :param corr_threshold: the maximal level of tolerable correlation.
        :param features_to_keep: list of column names of the features to keep.
        :param target_feature: column name of a target feature to predict further.

        :return: a DataFrame object M columns (where M <= N) with chemical descriptors from Rdkit
        + 1 column with disease categories.

        """

        corr_filtered_dataset_with_descriptors = None

        return corr_filtered_dataset_with_descriptors


if __name__ == '__main__':
    main()
