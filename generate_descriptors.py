{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "provenance": []
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "source": [
        "Импортируем библиотеки"
      ],
      "metadata": {
        "id": "sLgVkF6YJCJ8"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "pip install rdkit"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "F-mTKJArJiqC",
        "outputId": "3cd88018-66e2-49a5-8068-8423b27a297d"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Looking in indexes: https://pypi.org/simple, https://us-python.pkg.dev/colab-wheels/public/simple/\n",
            "Collecting rdkit\n",
            "  Downloading rdkit-2022.9.5-cp39-cp39-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (29.4 MB)\n",
            "\u001b[2K     \u001b[90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\u001b[0m \u001b[32m29.4/29.4 MB\u001b[0m \u001b[31m32.1 MB/s\u001b[0m eta \u001b[36m0:00:00\u001b[0m\n",
            "\u001b[?25hRequirement already satisfied: numpy in /usr/local/lib/python3.9/dist-packages (from rdkit) (1.22.4)\n",
            "Requirement already satisfied: Pillow in /usr/local/lib/python3.9/dist-packages (from rdkit) (8.4.0)\n",
            "Installing collected packages: rdkit\n",
            "Successfully installed rdkit-2022.9.5\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "7N3MIdjJXpja"
      },
      "outputs": [],
      "source": [
        "import pandas as pd\n",
        "import numpy as np\n",
        "from rdkit import Chem\n",
        "from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator"
      ]
    },
    {
      "cell_type": "markdown",
      "source": [
        "Создаём функцию generate_descriptors()"
      ],
      "metadata": {
        "id": "c7ATyp4RJFOg"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "def generate_descriptors(file_path,\n",
        "                             smiles_column_name, disease_column_name):\n",
        "        \"\"\"\n",
        "        :param file_path: path to the .csv table.\n",
        "        :param smiles_column_name: string reflecting the SMILES column name.\n",
        "        :return: a DataFrame object containing N columns with chemical descriptors from Rdkit\n",
        "        + 1 column with disease categories.\n",
        "        \"\"\"\n",
        "        Smiles_and_Disease = pd.read_csv(file_path)[[smiles_column_name, disease_column_name]]\n",
        "        Smiles_and_Disease = Smiles_and_Disease.drop_duplicates().dropna().reset_index(drop=True)\n",
        "        \n",
        "        descriptors = [x[0] for x in Chem.Descriptors._descList]\n",
        "        calculator = MolecularDescriptorCalculator(descriptors)\n",
        "\n",
        "        dataset_with_descriptors = []\n",
        "        \n",
        "        for smile in Smiles_and_Disease[smiles_column_name]:\n",
        "          mol = Chem.MolFromSmiles(smile)                             \n",
        "          dataset_with_descriptors.append(calculator.CalcDescriptors(mol))\n",
        "          \n",
        "        dataset_with_descriptors = pd.DataFrame(dataset_with_descriptors)\n",
        "        dataset_with_descriptors.columns = descriptors\n",
        "        \n",
        "        dataset_with_descriptors[\"Disease\"] = Smiles_and_Disease[disease_column_name]\n",
        "        dataset_with_descriptors = dataset_with_descriptors.dropna().reset_index(drop=True)                                              \n",
        "\n",
        "        #dataset_with_descriptors = None\n",
        "        \n",
        "        return dataset_with_descriptors"
      ],
      "metadata": {
        "id": "p6n5m4U0apIF"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "A = generate_descriptors(\"/content/full_drug-disease_db.csv\", \"iso_smiles\", \"disease category\")"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "KLVBJFO5PAak",
        "outputId": "7226f847-21ca-4ac4-b51e-29f4cd17d00d"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "<ipython-input-3-37817f98316a>:9: DtypeWarning: Columns (26,27,28,31,34) have mixed types. Specify dtype option on import or set low_memory=False.\n",
            "  Smiles_and_Disease = pd.read_csv(file_path)[[smiles_column_name, disease_column_name]]\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "A.isna().sum()"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "73AQ1qg6jPuV",
        "outputId": "a777271e-60fc-41d2-d639-bd1f993feea8"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "execute_result",
          "data": {
            "text/plain": [
              "MaxEStateIndex       0\n",
              "MinEStateIndex       0\n",
              "MaxAbsEStateIndex    0\n",
              "MinAbsEStateIndex    0\n",
              "qed                  0\n",
              "                    ..\n",
              "fr_thiocyan          0\n",
              "fr_thiophene         0\n",
              "fr_unbrch_alkane     0\n",
              "fr_urea              0\n",
              "Disease              0\n",
              "Length: 209, dtype: int64"
            ]
          },
          "metadata": {},
          "execution_count": 9
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "len(A)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "XF-j0Lpg2HQM",
        "outputId": "6d5c3675-4c64-4565-df48-aac0cd4f908c"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "execute_result",
          "data": {
            "text/plain": [
              "710"
            ]
          },
          "metadata": {},
          "execution_count": 10
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "print(A)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "_Ll4zOBoVTer",
        "outputId": "18966863-54b1-4109-d1e4-3d15e6f8666a"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "     MaxEStateIndex  MinEStateIndex  MaxAbsEStateIndex  MinAbsEStateIndex  \\\n",
            "0         10.309105       -1.064630          10.309105           0.106065   \n",
            "1         10.884259       -0.167824          10.884259           0.167824   \n",
            "2         12.578431       -1.755831          12.578431           0.010081   \n",
            "3         10.327917       -0.961296          10.327917           0.047454   \n",
            "4         12.881264       -0.751920          12.881264           0.109219   \n",
            "..              ...             ...                ...                ...   \n",
            "705       10.944576       -5.406512          10.944576           0.956528   \n",
            "706       10.944576       -5.406512          10.944576           0.956528   \n",
            "707       10.838410       -3.561044          10.838410           0.301574   \n",
            "708       10.838410       -3.561044          10.838410           0.301574   \n",
            "709       10.838410       -3.561044          10.838410           0.301574   \n",
            "\n",
            "          qed    MolWt  HeavyAtomMolWt  ExactMolWt  NumValenceElectrons  \\\n",
            "0    0.491543  163.198         154.126  163.030314                   58   \n",
            "1    0.519237  136.114         132.082  136.038511                   50   \n",
            "2    0.113064  585.608         542.264  585.285736                  234   \n",
            "3    0.532393  131.131         122.059  131.058243                   52   \n",
            "4    0.502326  408.882         383.682  408.145200                  152   \n",
            "..        ...      ...             ...         ...                  ...   \n",
            "705  0.431964  272.090         262.010  271.996324                   92   \n",
            "706  0.431964  272.090         262.010  271.996324                   92   \n",
            "707  0.789958  212.230         204.166  212.025563                   74   \n",
            "708  0.789958  212.230         204.166  212.025563                   74   \n",
            "709  0.789958  212.230         204.166  212.025563                   74   \n",
            "\n",
            "     NumRadicalElectrons  ...  fr_sulfonamd  fr_sulfone  fr_term_acetylene  \\\n",
            "0                      0  ...             0           0                  0   \n",
            "1                      0  ...             0           0                  0   \n",
            "2                      0  ...             0           0                  0   \n",
            "3                      0  ...             0           0                  0   \n",
            "4                      0  ...             0           0                  0   \n",
            "..                   ...  ...           ...         ...                ...   \n",
            "705                    0  ...             0           0                  0   \n",
            "706                    0  ...             0           0                  0   \n",
            "707                    0  ...             1           0                  0   \n",
            "708                    0  ...             1           0                  0   \n",
            "709                    0  ...             1           0                  0   \n",
            "\n",
            "     fr_tetrazole  fr_thiazole  fr_thiocyan  fr_thiophene  fr_unbrch_alkane  \\\n",
            "0               0            0            0             0                 0   \n",
            "1               0            0            0             0                 0   \n",
            "2               0            0            0             0                 0   \n",
            "3               0            0            0             0                 0   \n",
            "4               0            0            0             0                 1   \n",
            "..            ...          ...          ...           ...               ...   \n",
            "705             0            0            0             0                 0   \n",
            "706             0            0            0             0                 0   \n",
            "707             0            0            0             0                 0   \n",
            "708             0            0            0             0                 0   \n",
            "709             0            0            0             0                 0   \n",
            "\n",
            "     fr_urea                            Disease  \n",
            "0          0    Mental and behavioural disorder  \n",
            "1          0                Hematologic disease  \n",
            "2          0  Congenital disorder of metabolism  \n",
            "3          0                             Cancer  \n",
            "4          0             Cardiovascular disease  \n",
            "..       ...                                ...  \n",
            "705        0             Nervous system disease  \n",
            "706        0                             Cancer  \n",
            "707        0                              Other  \n",
            "708        0             Nervous system disease  \n",
            "709        0    Mental and behavioural disorder  \n",
            "\n",
            "[710 rows x 209 columns]\n"
          ]
        }
      ]
    }
  ]
}