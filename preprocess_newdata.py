import pandas as pd

def preprocess_input_data(
    file_path,
    log_transformer,
    scaler
):
    """
    Preprocess input data for model inference.
    
    Parameters
    ----------
    file_path : str
        Path to Excel file containing the raw input data.
    log_transformer : transformer object
        Fitted transformer for log transformation.
    scaler : transformer object
        Fitted scaler for numerical features.

    Returns
    -------
    pd.DataFrame
        Processed and scaled dataframe ready for inference.
    """

    # === 1. Load raw data ===
    df = pd.read_excel(file_path).drop(columns=["Unnamed: 0"], errors="ignore")

    # === 2. Columns after log transformation ===
    columns_processed = [
        "KMF", "Class I TBs", "Class II TBs", "Number of O", "Number of C", "H Bond Ratio",
        "SMRVSA5", "Number of N3", "MW", "Number of N", "Oxygen Balance", "NO2 adj to NO2",
        "NO2 adj to NH2", "NO2 adj to CO", "NO2 adj to CH3", "NO2 adj to OH", "NO2 adj to NH",
        "Rot Bonds", "Num Heteroatoms", "Total Rings", "Aromatic Rings", "Aliphatic Rings",
        "VSAEState8", "TPSA"
    ]

    df_processed = pd.DataFrame(
        log_transformer.transform(df),
        columns=columns_processed
    )

    # === 3. Feature selection ===
    to_drop = ["MW", "Number of O", "Rot Bonds", "Num Heteroatoms", "Total Rings"]

    columns_to_scale = [
        "KMF", "Class I TBs", "Class II TBs", "Number of N3", "H Bond Ratio",
        "Oxygen Balance", "Aromatic Rings", "Aliphatic Rings", "SMRVSA5", "TPSA",
        "VSAEState8", "Number of N"
    ]

    other_features = [
        col for col in df_processed.columns
        if col not in columns_to_scale and col not in to_drop
    ]

    all_feature_names = columns_to_scale + other_features

    # === 4. Scale features ===
    df_scaled = pd.DataFrame(
        scaler.transform(df_processed),
        columns=all_feature_names
    )

    # === 5. Rename columns (subscripts & formatting) ===
    rename_columns = {
        "SMRVSA5": "SMR_VSA5",
        "VSAEState8": "VSA_EState8",
        "Number of N3": "Number of N\u2083",
        "NO2 adj to NO2": "NO\u2082 adj to NO\u2082",
        "NO2 adj to NH2": "NO\u2082 adj. to NH\u2082",
        "NO2 adj to CO": "NO\u2082 adj. to CO",
        "NO2 adj to CH3": "NO\u2082 adj. to CH\u2083",
        "NO2 adj to OH": "NO\u2082 adj to OH",
        "NO2 adj to NH": "NO\u2082 adj to NH",
        "Aromatic Rings": "No. of Aromatic Rings",
        "Aliphatic Rings": "No. of Aliphatic Rings"
    }

    df_scaled.rename(columns=rename_columns, inplace=True)

    return df_scaled