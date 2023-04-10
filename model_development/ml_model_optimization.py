def main():
    def select_features(dataframe,
                        classification_model,
                        feature_selection_method,
                        num_of_final_features,
                        features_to_keep):
        """

        :param dataframe: table in the format of DataFrame object.
        :param classification_model: function calling a model to be used for feature set performance
        evaluation.
        :param feature_selection_method: function calling a feature selection method; can be forward
        or backward selection, recursive feature elimination, genetic algorithm.
        :param num_of_final_features: maximal number of features to select.
        :param features_to_keep: list of column names of the features to keep.

        :return:
        1) a DataFrame object M columns (where M < N) with selected chemical descriptors from Rdkit
        + 1 column with disease categories;
        2) model object trained on selected features

        """

        dataset_after_feature_selection = None
        trained_model = None

        return dataset_after_feature_selection, trained_model

    def optimize_model_parameters(model,
                                  grid_search_parameters):

        """

        :param model: model object to optimize.
        :param grid_search_parameters: dictionary containing model parameters names as the keys and
        3-element long tuples reflecting start, end, and step of parameter iteration.

        :return: optimized model object.

        """

        optimized_model = None

        return optimized_model


if __name__ == '__main__':
    main()
