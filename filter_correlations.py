def filter_correlations(dataframe,
                            max_corr_threshold,
                            mean_corr_threshold,
                            features_to_keep=''):
        """
        For every feature
        1) skips all the next steps if in features_to_keep
        2) counts mean and max correlations for Q3 quantile
        3) compares max and mean values with max_corr_threshold and mean_corr_threshold, respectively
        4) deletes if at least one number is bigger than the threshold
        :param dataframe: table in the format of DataFrame object.
        :param max_corr_threshold: the minimal level of correlation considered as "strong".
        :param mean_corr_threshold: the minimal level of mean correlation considered as "strong".
        :param features_to_keep: set of important chemical features to keep in any case (includes a target feature).
        :return: a DataFrame object M columns (where M <= N) with chemical descriptors from Rdkit
        + 1 column (target feature) with disease categories.
        """
        corr_filtered_dataset_with_descriptors = pd.DataFrame()
        Descriptors = dataframe.drop('Disease', axis=1)

        #Save columns from features_to_keep
        for i in range(len(features_to_keep)):
          corr_filtered_dataset_with_descriptors[features_to_keep[i]] = dataframe[features_to_keep[i]]

        #Save disease column
        Disease = dataframe['Disease']
        
        #Update initial columns set with features to filt
        features_to_filt = list(set(Descriptors.columns) - set(features_to_keep))
        DescriptorsToFilt = Descriptors[features_to_filt]

        #Create "feature-feature" matrix of Pearson's criterias
        PearsCoefMatrix = pd.DataFrame(index=range(len(features_to_filt)),columns=range(len(features_to_filt)))
        PearsCoefMatrix.columns, PearsCoefMatrix.index = features_to_filt, features_to_filt
        
        #Fill "feature-feature" matrix by Pearson's criterias values
        for i in range(len(features_to_filt) - 1):
          for j in range(i, len(features_to_filt)):
            PearsCoefMatrix.iloc[i][j] = abs(stats.pearsonr(DescriptorsToFilt[DescriptorsToFilt.columns[i]],
                                                            DescriptorsToFilt[DescriptorsToFilt.columns[j]])[0])
        PearsCoefMatrix = PearsCoefMatrix.replace(np.NaN, 0).transpose()
        
        #Find failed descriptors by max_corr_threshold and mean_corr_threshold
        Fail_descriptors = []
        for i in range(1,len(PearsCoefMatrix.columns)):
          if np.mean(PearsCoefMatrix[PearsCoefMatrix.columns[i]][i+1:]) > mean_corr_threshold or np.max(PearsCoefMatrix[PearsCoefMatrix.columns[i]][i+1:]) > max_corr_threshold:
            Fail_descriptors.append(PearsCoefMatrix.columns[i])

        #Saving final dataset
        final_descriptors = list(set(features_to_filt) - set(Fail_descriptors) | set(features_to_keep))
        corr_filtered_dataset_with_descriptors = Descriptors[final_descriptors]
        corr_filtered_dataset_with_descriptors["Disease"] = Disease

        return corr_filtered_dataset_with_descriptors